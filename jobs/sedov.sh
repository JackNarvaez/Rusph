#BSUB -J Sedov
#BSUB -W 70:00
#BSUB -n 40
#BSUB -q normal
#BSUB -o Sedov.%J
make Sedov
