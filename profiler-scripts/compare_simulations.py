#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys

def load_energy_data(file_path):
    """
    Loads energy data from a file and returns time and kinetic energy arrays.
    
    Parameters:
    - file_path: Path to the energy data file
    
    Returns:
    - time: numpy array of time values
    - kinetic_energy: numpy array of kinetic energy values
    """
    data = np.loadtxt(file_path, skiprows=1)  # Skip header row
    time = data[:, 0]  # First column is time
    kinetic_energy = data[:, 1]  # Third column is kinetic energy
    return time, kinetic_energy

def plot_kinetic_energy_difference(file1, file2, output_file):
    """
    Plots the difference in kinetic energy between two simulations over time.
    
    Parameters:
    - file1: Path to the first simulation's output file
    - file2: Path to the second simulation's output file
    - output_file: Path where the plot will be saved
    """
    # Load data from both files
    time1, ke1 = load_energy_data(file1)
    time2, ke2 = load_energy_data(file2)
    
    # Check if both files have the same number of time steps
    if len(time1) != len(time2):
        print("Error: The two files have different numbers of time steps.")
        return

    # Calculate the difference in kinetic energy
    ke_difference = ke1 - ke2

    # Plot the kinetic energy difference
    plt.figure(figsize=(10, 6))
    plt.plot(time1, ke_difference, label='Kinetic Energy Difference', marker='o')
    
    # Add plot labels and title
    plt.xlabel('Time')
    plt.ylabel('Kinetic Energy Difference')
    plt.title('Difference in Kinetic Energy between Simulations')
    plt.legend()
    plt.grid(True)

    # Save the plot to the specified file
    plt.savefig(output_file)
    plt.show()

    print(f"Kinetic energy difference plot saved to {output_file}")

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python3 compare_simulations.py <file1> <file2> <output_file>")
        sys.exit(1)

    # Get file paths and output file from command-line arguments
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3]

    # Generate the plot
    plot_kinetic_energy_difference(file1, file2, output_file)

