#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def plot_error_vs_time(differences_file, time_step_size, output_plot):
    """
    Reads the differences file, computes error norms at each time step,
    and plots the error versus time on a log-log plot.

    Parameters:
    - differences_file: Path to the file containing coordinate differences.
    - time_step_size: The time step size (h) used in the simulation.
    - output_plot: Filename for saving the output plot.
    """
    num_particles = 17576  # Fixed number of particles
    error_norms = []
    times = []

    try:
        # Read the differences file
        with open(differences_file, 'r') as f:
            lines = f.readlines()

        # Filter out only data lines that contain particle differences (lines that start with 'Line')
        data_lines = [line.strip() for line in lines if line.startswith('Line')]

        # Total number of data lines
        total_lines = len(data_lines)

        # Calculate the number of time steps
        num_time_steps = total_lines // num_particles

        if num_time_steps == 0:
            raise ValueError("No time steps detected in the file.")

        # Process data for each time step
        for step in range(num_time_steps):
            # Extract the relevant lines for the current time step
            start_idx = step * num_particles
            end_idx = start_idx + num_particles
            step_lines = data_lines[start_idx:end_idx]

            # Extract the coordinate differences and compute the error norm
            diffs = []
            for line in step_lines:
                # Split the line to extract the dx, dy, dz values
                tokens = line.split(':')[1].strip().split()
                coords = [float(token) for token in tokens]
                diffs.append(coords)

            # Convert to a numpy array for vectorized operations
            diffs_array = np.array(diffs)

            # Compute the RMS error norm for the time step
            error_norm = np.sqrt(np.mean(np.sum(diffs_array**2, axis=1)))
            error_norms.append(error_norm)

            # Record the time for the current time step
            time = (step + 1) * time_step_size
            times.append(time)

        # Convert lists to numpy arrays for plotting
        times = np.array(times)
        error_norms = np.array(error_norms)

        # Create the log-log plot of error vs. time
        plt.figure(figsize=(10, 6))
        plt.loglog(times, error_norms, marker='o', label='Error vs. Time')

        # Plot a reference line with slope 1 (for comparison)
        reference_error = error_norms[0] * (times / times[0])
        plt.loglog(times, reference_error, linestyle='--', label='Reference Slope 1')

        # Customize plot labels, title, and grid
        plt.xlabel('Time (s)')
        plt.ylabel('Error Norm')
        plt.title('Error vs. Time (Log-Log Plot)')
        plt.legend()
        plt.grid(True, which='both', linestyle='--')
        plt.tight_layout()

        # Save and show the plot
        plt.savefig(output_plot)
        plt.show()

        print(f"Plot saved to {output_plot}.")

    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage when running the script directly
if __name__ == "__main__":
    import sys
    # Check for correct number of arguments
    if len(sys.argv) != 4:
        print("Usage: python3 plot_error_vs_time.py <differences_file> <time_step_size> <output_plot>")
        sys.exit(1)

    # Get parameters from command-line arguments
    differences_file = sys.argv[1]
    time_step_size = float(sys.argv[2])
    output_plot = sys.argv[3]

    # Call the function to plot error vs. time
    plot_error_vs_time(differences_file, time_step_size, output_plot)

