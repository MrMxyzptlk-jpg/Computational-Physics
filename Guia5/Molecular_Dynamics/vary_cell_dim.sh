#!/bin/bash

start=5
end=8
increment=1

input_file="input.nml"
template_file="input_template.nml"
stats_file="datos/time_vs_cell_dim.dat"

cp "$input_file" "$template_file"

echo "# cell_dim | Elapsed time" > "$stats_file"

current=$start

while (( $(echo "$current <= $end" | bc -l) )); do
    echo "Running with cell_dim = $current"

    # Update input.nml
    awk -v rc="$current,$current,$current" '
    /^\s*cell_dim/ {
        print "    cell_dim             = " rc " "
        next
    }
    { print }
    ' "$template_file" > "$input_file"

    grep cell_dim "$input_file"

    # Run simulation
    ./run.exe > temp_output.txt

    # Extract energy stats
    grep 'Elapsed time' datos/INFO.out  > time_block.txt

    time=$(awk '/Elapsed time:/ {print $3}' time_block.txt)

    echo "$current  $time" >> "$stats_file"

    current=$(echo "$current + $increment" | bc -l)
done


rm time_block.txt
rm temp_output.txt