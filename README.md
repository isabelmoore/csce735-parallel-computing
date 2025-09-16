# csce735-parallel-computing

To log into HPRC:
```bash
ssh -p 22 isabelmoore717@grace.hprc.tamu.edu
```

To copy files over:

From the local folder:

```bash
scp -rv /home/wizard/CSCE735/HW2-735.zip isabelmoore717@grace.hprc.tamu.edu:/home/isabelmoore717/CSCE735/homework/

```

To check the status of the workers:

```bash

squeue -u $USER
```

To cancel a certain worker:
```bash
scancel <batch_job_number>
```
