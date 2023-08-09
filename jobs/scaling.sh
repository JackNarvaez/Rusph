#BSUB -J Speedup
#BSUB -W 40:00
#BSUB -n 40
#BSUB -q normal
#BSUB -o Speedup.%J
bash './speed_up/speedup.sh'
