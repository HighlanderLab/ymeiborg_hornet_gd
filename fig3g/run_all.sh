#!/bn/bash
for FILE in ./*.job
do
  sbatch $FILE
done
