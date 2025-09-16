//
// Sorts a list using multiple threads
//
#define _XOPEN_SOURCE 600   // enables POSIX barrier and clock_gettime
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#define MAX_THREADS     65536
#define MAX_LIST_SIZE   1000000000 // increased to handle case of 2^28 elements 

#define DEBUG 0

// Thread variables
//
// VS: ... declare thread variables, mutexes, condition varables, etc.,
// VS: ... as needed for this assignment 
//
static pthread_t *threads = NULL; // thread handles
static int       *tids = NULL; // thread IDs
static pthread_barrier_t barrier; // barrier for thread sync
static int *ptr_idx = NULL;   // start index of each original sublist

// Global variables
int num_threads;		// Number of threads to create - user input 
int list_size;			// List size
int *list;			// List of values
int *work;			// Work array
int *list_orig;			// Original list of values, used for error checking

// Print list - for debugging
void print_list(int *list, int list_size) {
    int i;
    for (i = 0; i < list_size; i++) {
        printf("[%d] \t %16d\n", i, list[i]); 
    }
    printf("--------------------------------------------------------------------\n"); 
}

// Comparison routine for qsort (stdlib.h) which is used to 
// a thread's sub-list at the start of the algorithm
int compare_int(const void *a0, const void *b0) {
    int a = *(int *)a0;
    int b = *(int *)b0;
    if (a < b) {
        return -1;
    } else if (a > b) {
        return 1;
    } else {
        return 0;
    }
}

// Return index of first element larger than or equal to v in sorted list
// ... return last if all elements are smaller than v
// ... elements in list[first], list[first+1], ... list[last-1]
//
//   int idx = first; while ((v > list[idx]) && (idx < last)) idx++;
//
int binary_search_lt(int v, int *list, int first, int last) {
   
    // Linear search code
    // int idx = first; while ((v > list[idx]) && (idx < last)) idx++; return idx;

    int left = first; 
    int right = last-1; 

    if (list[left] >= v) return left;
    if (list[right] < v) return right+1;
    int mid = (left+right)/2; 
    while (mid > left) {
        if (list[mid] < v) {
	    left = mid; 
	} else {
	    right = mid;
	}
	mid = (left+right)/2;
    }
    return right;
}
// Return index of first element larger than v in sorted list
// ... return last if all elements are smaller than or equal to v
// ... elements in list[first], list[first+1], ... list[last-1]
//
//   int idx = first; while ((v >= list[idx]) && (idx < last)) idx++;
//
int binary_search_le(int v, int *list, int first, int last) {

    // Linear search code
    // int idx = first; while ((v >= list[idx]) && (idx < last)) idx++; return idx;
 
    int left = first; 
    int right = last-1; 

    if (list[left] > v) return left; 
    if (list[right] <= v) return right+1;
    int mid = (left+right)/2; 
    while (mid > left) {
        if (list[mid] <= v) {
	    left = mid; 
	} else {
	    right = mid;
	}
	mid = (left+right)/2;
    }
    return right;
}
// --------------------- Thread worker ---------------------
static void *worker(void *arg) { 
    /*
    Purpose: Each thread sorts its own sublist and does only its own share of work. 
                Barriers were added so all threads stay in lockstep between phases,
                then repeatedly doubles size of sorted sublists. 
                Each thread first sorts its own sublist using qsort, then merges 
                sublists in parallel with other threads.
    Input:   arg = thread ID (0, 1, ..., num_threads-1)
    Output:  None
    Returns: NULL
    */
    int my_id = *(int *)arg; // my thread ID
    int np = list_size / num_threads; // sublist size
    int my_first = ptr_idx[my_id]; // first index of sublist
    int my_last = ptr_idx[my_id+1]; // last index of sublist

    // local qsort of this thread's leaf sublist
    qsort(&list[my_first], my_last - my_first, sizeof(int), compare_int);

    // all threads must finish qsort before merging
    pthread_barrier_wait(&barrier);

    // performs exactly q levels as list_size == (num_threads * np)
    // leaf size == np, doubling each level until reaching n
    int n = list_size;
    int q = 0; for (int s = np; s < n; s <<= 1) q++;

    // each thread merges 2 sublists of size (np * 2^level) at each level
    for (int level = 0; level < q; level++) {
        
        // size of each block being merged
        int my_blk_size = np * (1 << level); 
        // my block number
        int my_own_blk = ((my_id >> level) << level);
        // start index of my block
        int my_own_idx = ptr_idx[my_own_blk]; 
        // block to search
        int my_search_blk = ((my_id >> level) << level) ^ (1 << level); 
        // start index of block to search
        int my_search_idx = ptr_idx[my_search_blk];  
        // end index of block to search 
        int my_search_max = my_search_idx + my_blk_size;
        // limit search to list size
        int my_write_blk = ((my_id >> (level+1)) << (level+1));
        // start index of block to write into
        int my_write_idx = ptr_idx[my_write_blk]; //
        // current index into search block
        int idx = my_search_idx; 
        // number of elements found in search block
        int my_search_count = 0; 

        // first element through binary search
        if (my_search_blk > my_own_blk) {
            idx = binary_search_lt(list[my_first], list, my_search_idx, my_search_max);
        } else {
            idx = binary_search_le(list[my_first], list, my_search_idx, my_search_max);
        }
        my_search_count = idx - my_search_idx;
        int i_write = my_write_idx + my_search_count + (my_first - my_own_idx);
        work[i_write] = list[my_first];

        // remaining element - checks bound before reading list[idx]
        for (int i = my_first + 1; i < my_last; i++) {
            if (my_search_blk > my_own_blk) {
                while ((idx < my_search_max) && (list[i] > list[idx])) {
                    idx++; my_search_count++;
                }
            } else {
                while ((idx < my_search_max) && (list[i] >= list[idx])) {
                    idx++; my_search_count++;
                }
            }
            i_write = my_write_idx + my_search_count + (i - my_own_idx);
            work[i_write] = list[i];
        }

        // all threads done scattering into work[]
        pthread_barrier_wait(&barrier);

        // swap buffers once per level 
        if (my_id == 0) {
            int *tmp = list; 
            list = work; 
            work = tmp;
        }
        // ensure swapped pointers before next level
        pthread_barrier_wait(&barrier);

#if DEBUG
        if (my_id == 0) print_list(list, list_size);
        pthread_barrier_wait(&barrier);
#endif
    }

    return NULL;
}
// Sort list via parallel merge sort
//
// VS: ... to be parallelized using threads ...
//
void sort_list(int q) {

    int np = list_size/num_threads; // sublist size 

    // build ptr[] at start of each sublist
    ptr_idx = (int *)malloc((num_threads + 1) * sizeof(int));
    
    // Initialize starting position for each sublist
    for (int my_id = 0; my_id < num_threads; my_id++) {
        ptr_idx[my_id] = my_id * np;
    }
    ptr_idx[num_threads] = list_size;

    // threads and barrier
    threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t)); // thread handles
    tids    = (int *)malloc(num_threads * sizeof(int)); // thread IDs
    pthread_barrier_init(&barrier, NULL, num_threads); // barrier for thread sync

    // launch workers
    for (int t = 0; t < num_threads; t++) {
        tids[t] = t;
        pthread_create(&threads[t], NULL, worker, &tids[t]);
    }

    // join workers
    for (int t = 0; t < num_threads; t++) {
        pthread_join(threads[t], NULL);
    }

    // thread cleanup
    pthread_barrier_destroy(&barrier);
    free(threads); free(tids);
    free(ptr_idx);
}

// Main program - set up list of random integers and use threads to sort the list
//
// Input: 
//	k = log_2(list size), therefore list_size = 2^k
//	q = log_2(num_threads), therefore num_threads = 2^q
//
int main(int argc, char *argv[]) {

    struct timespec start, stop, stop_qsort;
    double total_time, time_res, total_time_qsort;
    int k, q, j, error; 

    // Read input, validate
    if (argc != 3) {
	printf("Need two integers as input \n"); 
	printf("Use: <executable_name> <log_2(list_size)> <log_2(num_threads)>\n"); 
	exit(0);
    }
    k = atoi(argv[argc-2]);
    if ((list_size = (1 << k)) > MAX_LIST_SIZE) {
	printf("Maximum list size allowed: %d.\n", MAX_LIST_SIZE);
	exit(0);
    }; 
    q = atoi(argv[argc-1]);
    if ((num_threads = (1 << q)) > MAX_THREADS) {
	printf("Maximum number of threads allowed: %d.\n", MAX_THREADS);
	exit(0);
    }; 
    if (num_threads > list_size) {
	printf("Number of threads (%d) < list_size (%d) not allowed.\n", 
	   num_threads, list_size);
	exit(0);
    }; 

    // Allocate list, list_orig, and work

    list = (int *) malloc(list_size * sizeof(int));
    list_orig = (int *) malloc(list_size * sizeof(int));
    work = (int *) malloc(list_size * sizeof(int));

//
// VS: ... May need to initialize mutexes, condition variables, 
// VS: ... and their attributes
//

    // Initialize list of random integers; list will be sorted by 
    // multi-threaded parallel merge sort
    // Copy list to list_orig; list_orig will be sorted by qsort and used
    // to check correctness of multi-threaded parallel merge sort    
    srand48(0); 	// seed the random number generator
    for (j = 0; j < list_size; j++) {
	list[j] = (int) lrand48();
	list_orig[j] = list[j];
    }
    // duplicate first value at last location to test for repeated values
    list[list_size-1] = list[0]; list_orig[list_size-1] = list_orig[0];

    // Create threads; each thread executes find_minimum
    clock_gettime(CLOCK_REALTIME, &start);

//
// VS: ... may need to initialize mutexes, condition variables, and their attributes
//

// Serial merge sort 
// VS: ... replace this call with multi-threaded parallel routine for merge sort
// VS: ... need to create threads and execute thread routine that implements 
// VS: ... parallel merge sort

    sort_list(q);

    // Compute time taken
    clock_gettime(CLOCK_REALTIME, &stop);
    total_time = (stop.tv_sec-start.tv_sec)
	+0.000000001*(stop.tv_nsec-start.tv_nsec);

    // Check answer
    qsort(list_orig, list_size, sizeof(int), compare_int);
    clock_gettime(CLOCK_REALTIME, &stop_qsort);
    total_time_qsort = (stop_qsort.tv_sec-stop.tv_sec)
	+0.000000001*(stop_qsort.tv_nsec-stop.tv_nsec);

    error = 0; 
    for (j = 1; j < list_size; j++) {
	if (list[j] != list_orig[j]) error = 1; 
    }

    if (error != 0) {
	printf("Houston, we have a problem!\n"); 
    }

    // Print time taken
    printf("List Size = %d, Threads = %d, error = %d, time (sec) = %8.4f, qsort_time = %8.4f\n", 
	    list_size, num_threads, error, total_time, total_time_qsort);

// VS: ... destroy mutex, condition variables, etc.

    free(list); free(work); free(list_orig); 

}

