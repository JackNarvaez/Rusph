#BSUB -J Kelvin_Helmholtz
#BSUB -W 70:00
#BSUB -n 40
#BSUB -q normal
#BSUB -o KelvinHelmholtz.%J

make Kelvinhelmholtz
