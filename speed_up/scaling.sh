Np=6
for ii in $(seq 0 ${Np}); do
     s=$((2**(6-ii)))
     d=$(RAYON_NUM_THREADS=${s} cargo run -p speed_up);
     echo ${s} $d
done > "./Data/speed_up/turbulence/scaling64K.txt"