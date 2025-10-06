# csce735-parallel-computing

## SSH into Grace HRPC:
```bash
ssh -p 22 isabelmoore717@grace.hprc.tamu.edu
```

## Transfer files from local to remote:
```bash
scp -rv /home/wizard/CSCE735/HW2-735.zip isabelmoore717@grace.hprc.tamu.edu:/home/isabelmoore717/CSCE735/homework/
```

## Check the number of tasks running:
```bash
squeue -u $USER
```

## Cancel a particular task:
```bash
scancel 16643980
```
