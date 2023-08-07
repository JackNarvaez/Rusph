#BSUB -J Scaling
#BSUB -W 40:00
#BSUB -n 70
#BSUB -q normal
#BSUB -o Scaling.%J
bash './speed_up/scaling.sh'
