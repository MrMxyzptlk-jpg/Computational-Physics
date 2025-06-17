#!/bin/bash

start=0.1
end=0.88
increment=0.02

input_file="input.nml"
template_file="input_template.nml"
stats_file="datos/pressure_vs_density.dat"

cp "$input_file" "$template_file"

echo "# density  | Pressure_avg  |  Pressure_std  | Sfactor_avg  |  Sfactor_std | Diffusion*6 | Diffusion_stddev*6" > "$stats_file"

current=$start

while (( $(echo "$current <= $end" | bc -l) )); do
    echo "Running with density = $current"
    # Update input.nml
    awk -v rc="$current" '
    /^\s*density/ {
        print "    density             = " rc "           ! Number of atoms/volume"
        next
    }
    { print }
    ' "$template_file" > "$input_file"

    grep density "$input_file"

    # Run simulation
    ./run.exe > temp_output.txt

    # Extract pressure stats
    grep 'Pressure' datos/INFO.out > Density_Block.txt
    Pressure_avg=$(awk '/Pressure:/ {print $4}' Density_Block.txt)
    Pressure_std=$(awk '/Pressure:/ {print $8}' Density_Block.txt)
    # Extract Structure Factors stats
    grep 'Structure Factor' datos/INFO.out > Density_Block.txt
    Sfactor_avg=$(awk '/Structure Factor:/ {print $5}' Density_Block.txt)
    Sfactor_std=$(awk '/Structure Factor:/ {print $9}' Density_Block.txt)

    # For Mean Square Displacement
    gnuplot plot_msd.gp
    grep 'Final set of parameters ' -A 3 fit.log > Density_Block.txt
    Difusion_avg=$(awk '/a               =/ {print $3}' Density_Block.txt)
    Difusion_std=$(awk '/a               =/ {print $5}' Density_Block.txt)
    rm fit.log

    echo "$current  $Pressure_avg  $Pressure_std  $Sfactor_avg  $Sfactor_std  $Difusion_avg  $Difusion_std" >> "$stats_file"


    cp "datos/mean_sqr_displacement.out" "datos/msd_density_$current.dat"

    gnuplot plot_observables.gp
    cp PLOT.png "figs/observables_density_$current.png"
    rm fit.log

    current=$(echo "$current + $increment" | bc -l)
done


rm Density_Block.txt
rm temp_output.txt