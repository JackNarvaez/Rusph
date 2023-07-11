#BSUB -J Scaling
#BSUB -W 3:00
#BSUB -n 42
#BSUB -o Scaling.%J
bash '../speed_up/scaling.sh'