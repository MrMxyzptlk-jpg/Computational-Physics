#!/bin/bash

start=0.2
end=1.05
increment=0.05

input_file="input.nml"
template_file="input_template.nml"
stats_file="datos/pressure_vs_density.dat"

cp "$input_file" "$template_file"

echo "# density  |  D*D0 | Fit error" > "$stats_file"

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

    # For Mean Square Displacement
    gnuplot plot_msd.gp
    grep 'Final set of parameters ' -A 3 fit.log > Density_Block.txt
    Difusion_avg=$(awk '/a               =/ {print $3}' Density_Block.txt)
    Difusion_std=$(awk '/a               =/ {print $5}' Density_Block.txt)
    rm fit.log

    echo "$current  $Difusion_avg  $Difusion_std" >> "$stats_file"


    cp "datos/mean_sqr_displacement.out" "datos/msd_density_$current.dat"

    cp PLOT.png "figs/msd_density_$current.png"
    rm fit.log

    current=$(echo "$current + $increment" | bc -l)
done


rm Density_Block.txt
rm temp_output.txt