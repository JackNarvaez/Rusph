#BSUB -J Turbulence
#BSUB -W 70:00
#BSUB -n 40
#BSUB -q normal
#BSUB -o Turbulence.%J
bash './turbulent_gas/turbulence.sh'