#BSUB -J Sedov
#BSUB -W 70:00
#BSUB -n 40
#BSUB -q gpu
#BSUB -o Sedov.%J
make Sedov
