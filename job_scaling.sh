#BSUB -J Scaling
#BSUB -W 1:30
#BSUB -n 40
#BSUB -o Scaling.%J
bash scaling.sh