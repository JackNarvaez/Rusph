#BSUB -J Kelvin_Helmholtz
#BSUB -W 70:00
#BSUB -n 40
#BSUB -q gpu
#BSUB -o KelvinHelmholtz.%J
bash './tests/kelvin_helmholtz/kh.sh'
