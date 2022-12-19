#!/bn/bash
rm *_out*
rm *_err*
for FILE in ./*.job
do
  sbatch $FILE
done
