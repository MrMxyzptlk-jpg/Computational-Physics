#!/bin/bash

start=1.0
end=5.0
increment=0.1

input_file="input.nml"
template_file="input_template.nml"
stats_file="datos/energy_vs_rcut.dat"

cp "$input_file" "$template_file"

echo "# rcut  E_pot_avg  E_pot_std  E_kin_avg  E_kin_std  E_tot_avg  E_tot_std" > "$stats_file"

current=$start

while (( $(echo "$current <= $end" | bc -l) )); do
    echo "Running with radius_cutoff = $current"

    # Update input.nml
    awk -v rc="$current" '
    /^\s*radius_cutoff/ {
        print "    radius_cutoff      = " rc "              ! In Bohr"
        next
    }
    { print }
    ' "$template_file" > "$input_file"

    grep radius_cutoff "$input_file"

    # Run simulation
    ./run.exe > temp_output.txt

    # Extract energy stats
    grep 'Energy' datos/INFO.out > energy_block.txt

    Epot_avg=$(awk '/Potential Energy:/ {print $5}' energy_block.txt)
    Epot_std=$(awk '/Potential Energy:/ {print $9}' energy_block.txt)

    Ekin_avg=$(awk '/Kinetic Energy:/ {print $5}' energy_block.txt)
    Ekin_std=$(awk '/Kinetic Energy:/ {print $9}' energy_block.txt)

    Etot_avg=$(awk '/Total Energy:/ {print $5}' energy_block.txt)
    Etot_std=$(awk '/Total Energy:/ {print $9}' energy_block.txt)

    echo "$current  $Epot_avg  $Epot_std  $Ekin_avg  $Ekin_std  $Etot_avg  $Etot_std" >> "$stats_file"

    gnuplot plot_errors.gp

    current=$(echo "$current + $increment" | bc -l)
done


rm energy_block.txt
rm temp_output.txt