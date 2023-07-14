#BSUB -J Scaling
#BSUB -W 0:40
#BSUB -n 40
#BSUB -q normal
#BSUB -o Scaling.%J
bash './speed_up/scaling.sh'
