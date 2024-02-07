#BSUB -J ToyStar
#BSUB -W 0:10
#BSUB -n 20
#BSUB -q normal
#BSUB -o ToyStar.%J
make Toystar
