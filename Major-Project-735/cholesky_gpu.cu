// Compute triangular factors of an SPD matrix using GPU
// - given an SPD matrix A, compute upper triangular matrix R such that
//   A = R'*R, where R' is the transpose of R
//
// SPD = Symmetric Positive Definite, doesn't require pivoting during factorization
//
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <new>
#include <cuda_runtime.h>
#include <cublas_v2.h>

#define MAX_MATRIX_SIZE 4096
#define TOL 1.0e-8

#define ERR_MALLOC 1
#define ERR_MEMCPY 2
#define ERR_KERNEL 3

#define DEBUG 1

// Define Matrix
typedef struct {
    int    n;          // order of matrix (number of rows and columns)
    double **elements; // host view: 2D pointer for convenience
    double *array;     // linear array, stores matrix row-by-row
} Matrix;

// ============================================================================
// Device routines
// ============================================================================

__global__ void device_create_matrix_on_device(Matrix A);

// New tiled kernels
// Tile size (adjustable)
#ifndef TILE
#define TILE 16
#endif

__global__ void potrf_tile_kernel(Matrix A, int k, int B);
__global__ void trsm_tile_kernel(Matrix A, int k, int B);
__global__ void syrk_gemm_update_kernel(Matrix A, int k, int B);

// Initializes A.elements of matrix A to point to the start of each row in A.array
// - allows access to matrix entries in the standard way:
//   A.elements[i][j] points to Aij element that is stored in A.array
// - both A.elements and A.array are allocated on the device, therefore not
//   accessible to host
//
__global__ void device_create_matrix_on_device(Matrix A) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < A.n) {
        A.elements[i] = &(A.array[i * A.n]);
    }
}

// potrf on a diagonal tile stored row-major in A.array
// Each diagonal tile is BxB at block index k. This kernel is launched with
// one block per diagonal tile. It uses shared memory and performs a small
// serial Cholesky on the tile (thread 0 does the work) then writes back.
__global__ void potrf_tile_kernel(Matrix A, int k, int B) {
    int n = A.n;
    int bi = k * B; // start row/col of tile

    extern __shared__ double s[]; // shared tile B x B
    double *sA = s; // size B*B (may contain unused entries if partial tile)

    // Initialize shared tile with global contents (handle partial tiles)
    for (int i = threadIdx.x; i < B * B; i += blockDim.x) {
        int r = i / B;
        int c = i % B;
        int gr = bi + r;
        int gc = bi + c;
        double val = 0.0;
        if (gr < n && gc < n) val = A.array[(size_t)gr * n + gc];
        sA[r * B + c] = val;
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        // serial Cholesky on shared tile (upper-triangular output)
        int mb = (B < (n - bi)) ? B : (n - bi);
        for (int t = 0; t < mb; ++t) {
            double pivot = sA[t * B + t];
            pivot = sqrt(pivot);
            sA[t * B + t] = pivot;
            for (int j = t + 1; j < mb; ++j) {
                sA[t * B + j] = sA[t * B + j] / pivot;
            }
            for (int i = t + 1; i < mb; ++i) {
                for (int j = t + 1; j < mb; ++j) {
                    sA[i * B + j] -= sA[t * B + j] * sA[t * B + i];
                }
            }
            for (int i = t + 1; i < mb; ++i) {
                sA[i * B + t] = 0.0;
            }
        }
    }
    __syncthreads();

    // Write back shared tile into global A
    for (int i = threadIdx.x; i < B * B; i += blockDim.x) {
        int r = i / B;
        int c = i % B;
        int gr = bi + r;
        int gc = bi + c;
        if (gr < n && gc < n) {
            A.array[(size_t)gr * n + gc] = sA[r * B + c];
        }
    }
}

// TRSM for tiles: solve Rkk^T * X = Btile  (Rkk is upper triangular)
// We treat Rkk^T as lower triangular and perform forward substitution.
// Each block handles one tile in the k-th row (tile columns j>k).
__global__ void trsm_tile_kernel(Matrix A, int k, int B) {
    int n = A.n;
    int bi = k * B;
    int tile_j = k + 1 + blockIdx.x; // column block index
    if (tile_j * B >= n) return;

    int br = threadIdx.x; // column within tile (we parallelize over columns)

    // Load Rkk (upper) into shared memory
    __shared__ double sR[TILE * TILE];
    __shared__ double sB[TILE * TILE];

    int mb = (B < (n - bi)) ? B : (n - bi);
    int nb = (B < (n - tile_j * B)) ? B : (n - tile_j * B);

    // load Rkk
    for (int i = threadIdx.x; i < mb * mb; i += blockDim.x) {
        int r = i / mb;
        int c = i % mb;
        sR[r * mb + c] = A.array[(size_t)(bi + r) * n + (bi + c)];
    }

    // load Btile (row k, col tile_j)
    for (int i = threadIdx.x; i < mb * nb; i += blockDim.x) {
        int r = i / nb;
        int c = i % nb;
        sB[r * nb + c] = A.array[(size_t)(bi + r) * n + (tile_j * B + c)];
    }
    __syncthreads();

    // Forward substitution per column (each thread handles a column index c)
    for (int col = br; col < nb; col += blockDim.x) {
        for (int i = 0; i < mb; ++i) {
            double sum = sB[i * nb + col];
            for (int t = 0; t < i; ++t) {
                // Rkk^T has element L[i,t] = sR[t*mb + i]
                sum -= sR[t * mb + i] * sB[t * nb + col];
            }
            // diagonal of L is sR[i*mb + i]
            sB[i * nb + col] = sum / sR[i * mb + i];
        }
    }
    __syncthreads();

    // Write back computed R[k, tile_j]
    for (int i = threadIdx.x; i < mb * nb; i += blockDim.x) {
        int r = i / nb;
        int c = i % nb;
        A.array[(size_t)(bi + r) * n + (tile_j * B + c)] = sB[r * nb + c];
    }
}

// Update trailing blocks C(bi,bj) -= R(k,bi)^T * R(k,bj)
// Launch grid over trailing block rows (y) and cols (x)
__global__ void syrk_gemm_update_kernel(Matrix A, int k, int B) {
    int n = A.n;
    int bi_idx = k + 1 + blockIdx.y; // block-row index
    int bj_idx = k + 1 + blockIdx.x; // block-col index
    if (bi_idx * B >= n || bj_idx * B >= n) return;

    int bi = bi_idx * B;
    int bj = bj_idx * B;
    int bk = k * B;

    int mb = (B < (n - bi)) ? B : (n - bi);
    int nb = (B < (n - bj)) ? B : (n - bj);
    int kb = (B < (n - bk)) ? B : (n - bk);

    // Each thread computes one element within the mb x nb tile
    int local_i = threadIdx.y;
    int local_j = threadIdx.x;

    // Bounds check so partial tiles don't read/write out of range
    if (local_i >= mb || local_j >= nb) {
        return;
    }

    double sum = 0.0;
    for (int t = 0; t < kb; ++t) {
        double a = A.array[(size_t)(bk + t) * n + (bi + local_i)];
        double b = A.array[(size_t)(bk + t) * n + (bj + local_j)];
        sum += a * b;
    }

    int gr = bi + local_i;
    int gc = bj + local_j;
    if (gr < n && gc < n) {
        A.array[(size_t)gr * n + gc] -= sum;
    }
}

// ============================================================================
// Host routines
// ============================================================================

Matrix cholesky_factorization(Matrix&);
Matrix product_with_transpose(Matrix& R);
int compare_matrix(Matrix&, Matrix&);
Matrix clone_matrix(Matrix& A);
void initialize_spd_matrix(Matrix&, double);
Matrix create_matrix(int, int);
void free_matrix_memory(Matrix&);
void print_matrix(Matrix&);
void check_error(cudaError_t, int);
void print_device_properties();

// Cholesky factorization on the host (provided for reference only)
// - return R where R is an upper triangular matrix
//   such that A = R'*R (R' = transpose of R)
//
Matrix cholesky_factorization(Matrix& A) {
    double sqrt_pivot;
    Matrix R = clone_matrix(A);
    for (int k = 0; k < R.n; k++) {
        sqrt_pivot = sqrt(R.elements[k][k]);
        for (int j = k; j < R.n; j++) {
            R.elements[k][j] = R.elements[k][j] / sqrt_pivot;
        }
        for (int i = k + 1; i < R.n; i++) {
            for (int j = k + 1; j < R.n; j++) {
                R.elements[i][j] -= R.elements[k][j] * R.elements[k][i];
            }
        }
        for (int j = k + 1; j < R.n; j++) {
            R.elements[j][k] = 0.0;
        }
    }
    return R;
}

// Product with transpose
// - return C = R' * R
Matrix product_with_transpose(Matrix& R) {
    Matrix C = create_matrix(R.n, R.n);
    for (int i = 0; i < R.n; i++) {
        for (int j = 0; j < R.n; j++) {
            C.elements[i][j] = 0.0;
            for (int k = 0; k < R.n; k++)
                C.elements[i][j] += R.elements[k][i] * R.elements[k][j];
        }
    }
    return C;
}

// Compare if matrix is identical to another by checking if
// their elements are identical within specified tolerance
int compare_matrix(Matrix& A, Matrix& B) {
    int error = 0;
    if (A.n != B.n) return 1;
    for (int i = 0; i < A.n; i++) {
        for (int j = 0; j < A.n; j++) {
            if (fabs(A.elements[i][j] - B.elements[i][j]) > TOL) {
                error = 1;
                break;
            }
        }
        if (error) break;
    }
    return error;
}

// Clone matrix
Matrix clone_matrix(Matrix& A) {
    Matrix C = create_matrix(A.n, A.n);
    for (int i = 0; i < C.n; i++) {
        for (int j = 0; j < C.n; j++) {
            C.elements[i][j] = A.elements[i][j];
        }
    }
    return C;
}

// Initialize an SPD matrix (for testing factorization routine)
void initialize_spd_matrix(Matrix& A, double delta) {
    double value;
    for (int i = 0; i < A.n; i++) {
        A.elements[i][i] = delta;
    }
    for (int i = 0; i < A.n; i++) {
        for (int j = i + 1; j < A.n; j++) {
            value = (double)(rand()) / (double)(RAND_MAX);
            A.elements[i][j] = value;
            A.elements[j][i] = value;
            A.elements[i][i] += fabs(A.elements[i][j]);
            A.elements[j][j] += fabs(A.elements[i][j]);
        }
    }
}

// Create new matrix (square)
// - matrix entries are uninitialized
Matrix create_matrix(int num_rows, int num_cols) {
    Matrix A;
    A.n = num_rows;  // square; assume num_rows == num_cols
    A.elements = new double *[A.n];
    A.array    = new double[A.n * A.n];
    for (int i = 0; i < A.n; i++)
        A.elements[i] = &(A.array[i * A.n]);
    return A;
}

// Delete matrix arrays
void free_matrix_memory(Matrix& A) {
    delete [] A.elements;
    delete [] A.array;
}

// Print matrix (for debugging)
void print_matrix(Matrix& A) {
    printf("\n... Printing matrix ... \n");
    for (int i = 0; i < A.n; i++) {
        for (int j = 0; j < A.n; j++) {
            printf(" %8.4f", A.elements[i][j]);
        }
        printf("\n");
    }
}

// Generic error
void check_error(cudaError_t err, int type) {
    if (err != cudaSuccess) {
        switch (type) {
            case ERR_MALLOC:
                fprintf(stderr, "Failed cudaMalloc (error code %s)!\n",
                        cudaGetErrorString(err));
                break;
            case ERR_MEMCPY:
                fprintf(stderr, "Failed cudaMemcpy (error code %s)!\n",
                        cudaGetErrorString(err));
                break;
            case ERR_KERNEL:
                fprintf(stderr, "Failed kernel launch (error code %s)!\n",
                        cudaGetErrorString(err));
                break;
        }
        exit(0);
    }
}

// Print device properties
void print_device_properties() {
    int i, deviceCount;
    cudaDeviceProp deviceProp;
    cudaGetDeviceCount(&deviceCount);
    printf("------------------------------------------------------------\n");
    printf("Number of GPU devices found = %d\n", deviceCount);
    for (i = 0; i < deviceCount; ++i) {
        cudaGetDeviceProperties(&deviceProp, i);
        printf("[Device: %1d] Compute Capability %d.%d.\n",
               i, deviceProp.major, deviceProp.minor);
        printf(" ... multiprocessor count  = %d\n", deviceProp.multiProcessorCount);
        printf(" ... max threads per multiprocessor = %d\n",
               deviceProp.maxThreadsPerMultiProcessor);
        printf(" ... max threads per block = %d\n",
               deviceProp.maxThreadsPerBlock);
        printf(" ... max block dimension   = %d, %d, %d (along x, y, z)\n",
               deviceProp.maxThreadsDim[0],
               deviceProp.maxThreadsDim[1],
               deviceProp.maxThreadsDim[2]);
        printf(" ... max grid size         = %d, %d, %d (along x, y, z)\n",
               deviceProp.maxGridSize[0],
               deviceProp.maxGridSize[1],
               deviceProp.maxGridSize[2]);
        printf(" ... warp size             = %d\n", deviceProp.warpSize);
        printf(" ... clock rate            = %d MHz\n", deviceProp.clockRate / 1000);
    }
    printf("------------------------------------------------------------\n");
}

// ============================================================================
// Main Program
// ============================================================================

int main(int argc, char *argv[]) {
    cudaError_t err = cudaSuccess;

    // Timing variables
    cudaEvent_t start, stop;    // GPU timing variables
    float time_array[5];

    // Print device properties
    print_device_properties();

    // Get device information and set device to use
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount > 0) {
        cudaSetDevice(0);
    } else {
        printf("Warning: No GPU device found ... results may be incorrect\n");
    }

    // Timing initializations
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Read input, validate
    if (argc != 2) {
        printf("Need one integer as input \n");
        printf("Use: <executable_name> <matrix_size>\n");
        exit(0);
    }
    int matrix_size = atoi(argv[argc - 1]);
    if (matrix_size > MAX_MATRIX_SIZE) {
        printf("Maximum matrix size allowed: %d.\n", MAX_MATRIX_SIZE);
        exit(0);
    }

    // Initialize matrix A
    Matrix A = create_matrix(matrix_size, matrix_size);
    initialize_spd_matrix(A, 1.0);

    // Create a copy of A on the device:
    // - initialize number of rows of the copy
    // - allocate dA.elements and dA.array on device
    // - copy A.array that has matrix entries to the device
    // - initialize dA.elements on the device
    Matrix dA;
    dA.n = A.n;

    // Allocate linear arrays on device
    size_t size_elements = dA.n * sizeof(double*);
    size_t size_array    = (size_t)dA.n * dA.n * sizeof(double);
    err = cudaMalloc(&dA.elements, size_elements); check_error(err, ERR_MALLOC);
    err = cudaMalloc(&dA.array, size_array);      check_error(err, ERR_MALLOC);

    // Copy matrix elements to device
    err = cudaMemcpy(dA.array, A.array, size_array, cudaMemcpyHostToDevice);
    check_error(err, ERR_MEMCPY);

    // Initialize row pointer array dA.elements on device
    {
        int blockSize = 256;
        int gridSize  = (dA.n + blockSize - 1) / blockSize;
        device_create_matrix_on_device<<<gridSize, blockSize>>>(dA);
        err = cudaGetLastError(); check_error(err, ERR_KERNEL);
        cudaDeviceSynchronize();
    }

    // ------------------------------------------------------------------------
    // Parallel Cholesky factorization on device:
    // Host loop over k, multiple parallel kernels per iteration
    // ------------------------------------------------------------------------

    cudaEventRecord(start, 0);

    int n = dA.n;
    int B = TILE; // tile size
    int nb = (n + B - 1) / B; // number of tiles

    // Blocked tiled Cholesky
    for (int kb = 0; kb < nb; ++kb) {
        // 1. Factor diagonal tile (kb,kb)
        size_t shared_bytes = (size_t)B * B * sizeof(double);
        potrf_tile_kernel<<<1, 128, shared_bytes>>>(dA, kb, B);
        err = cudaGetLastError(); check_error(err, ERR_KERNEL);
        cudaDeviceSynchronize();

        int trailing = nb - kb - 1;
        if (trailing > 0) {
            // 2. Solve for row tiles R[kb, j] for j = kb+1..nb-1
            int threads_trsm = 128;
            trsm_tile_kernel<<<trailing, threads_trsm>>>(dA, kb, B);
            err = cudaGetLastError(); check_error(err, ERR_KERNEL);
            cudaDeviceSynchronize();

            // 3. Update trailing block submatrix: for bi,bj in trailing tiles
            dim3 block2D(TILE, TILE);
            dim3 grid2D(trailing, trailing);
            syrk_gemm_update_kernel<<<grid2D, block2D>>>(dA, kb, B);
            err = cudaGetLastError(); check_error(err, ERR_KERNEL);
            cudaDeviceSynchronize();
        }
    }

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&(time_array[0]), start, stop);

    // Copy result matrix from device to host
    Matrix R = create_matrix(dA.n, dA.n);
    size_array = (size_t)R.n * R.n * sizeof(double);

    err = cudaMemcpy(R.array, dA.array, size_array, cudaMemcpyDeviceToHost);
    check_error(err, ERR_MEMCPY);

    // enforce upper-triangular R
    for (int i = 0; i < R.n; ++i) {
        for (int j = 0; j < i; ++j) {   // j < i  => below diagonal
            R.elements[i][j] = 0.0;
        }
    }
    // Compute C = R'*R
    Matrix RtR = product_with_transpose(R);

    // Compare A with C = R'*R
    int error = compare_matrix(A, RtR);

    if (error != 0) {
        printf("+++  Houston, we have a problem!\n");
    } else {
        printf("+++  Matrix successfully factored\n");
        printf("Matrix size: %d, GPU execution time (parallel): %8.4f ms\n",
               A.n, time_array[0]);
    }

    // Free allocated arrays for A, R, RtR on host
    free_matrix_memory(A);
    free_matrix_memory(R);
    free_matrix_memory(RtR);

    // Free allocated arrays for dA on device
    cudaFree(dA.elements);
    cudaFree(dA.array);

    return 0;
}
