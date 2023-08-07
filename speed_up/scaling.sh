Np=6
for ii in $(seq 1 ${Np}); do
     s=$((2**ii))
     d=$(RAYON_NUM_THREADS=${s} cargo run -p speed_up);
     echo ${s} $d
done > "./Data/speed_up/scaling.txt"
