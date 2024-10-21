#!/bin/bash

# run_scaling_study.sh

# Ensure the script is executable
# chmod +x run_scaling_study.sh

# Clean previous builds
make realclean

# Build the executable
make sph.x

# Check if the build was successful
if [ ! -f sph.x ]; then
    echo "Error: Build failed."
    exit 1
fi

# Particle counts to test (adjust as needed)
particle_counts=(500)

# Simulation parameters (adjust as needed)
nframes=1
npframe=100

# Output directory for timing data
mkdir -p timing_results

# Loop over different particle counts
for n_particles in "${particle_counts[@]}"; do
    echo "Running simulation with approximately $n_particles particles..."

    # Calculate the h value to achieve the desired number of particles
    h_value=$(awk -v np="$n_particles" 'BEGIN {print 1.0 / (np^(1/3)) * 1.3}')

    # Output file name
    output_file="run_${n_particles}.out"

    # Run the simulation
    ./sph.x -s $h_value -o $output_file -F $nframes -f $npframe

    # Check if the simulation ran successfully
    if [ $? -ne 0 ]; then
        echo "Error: Simulation failed for n_particles = $n_particles."
        continue
    fi

    # Move the timing data file to the results directory
    timing_data_file="timing_data.csv"
    if [ -f $timing_data_file ]; then
        mv $timing_data_file "timing_results/timing_data_${n_particles}.csv"
    else
        echo "Warning: Timing data file not found for n_particles = $n_particles."
    fi

    # Optionally, save the output file
    mv $output_file "timing_results/$output_file"

done

echo "Timing experiments completed. Timing data files are in the timing_results/ directory."

