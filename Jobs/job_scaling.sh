#BSUB -J Scaling
#BSUB -W 2:30
#BSUB -n 40
#BSUB -q gpu
#BSUB -o Scaling.%J
bash './speed_up/scaling.sh'
