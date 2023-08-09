Np=(40 32 16 8 4 2 1)
for ii in ${!Np[@]}; do
     s=${Np[$ii]}
     d=$(RAYON_NUM_THREADS=${s} cargo run -p speed_up);
     echo ${s} $d
done > "./Data/speed_up/turbulence/scaling32K.txt"