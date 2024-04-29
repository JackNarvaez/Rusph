#BSUB -J AccDisc
#BSUB -W 70:00
#BSUB -n 40
#BSUB -q normal
#BSUB -o AccDisc.%J
make Accretiondisc
