#!/usr/bin/env python3

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    # Directory containing the timing data files
    timing_results_dir = 'timing_results'

    # Find all timing_data_*.csv files in the directory
    timing_files = glob.glob(os.path.join(timing_results_dir, 'timing_data_*.csv'))

    # Initialize a dictionary to hold the data
    data = {}

    # Process each timing data file
    for filepath in timing_files:
        # Extract the number of particles from the filename
        filename = os.path.basename(filepath)
        n_particles_str = filename.replace('timing_data_', '').replace('.csv', '')
        n_particles = int(n_particles_str)

        # Read the CSV file into a DataFrame
        df = pd.read_csv(filepath)

        # Add the number of particles to the DataFrame
        df['n_particles'] = n_particles

        # Append the DataFrame to the data dictionary
        data[n_particles] = df

    # Combine all DataFrames into one
    combined_df = pd.concat(data.values(), ignore_index=True)

    # Pivot the DataFrame to have functions as columns
    pivot_df = combined_df.pivot(index='n_particles', columns='Function', values='Time(s)')

    # Sort the DataFrame by number of particles
    pivot_df = pivot_df.sort_index()

    # Plotting
    plt.figure(figsize=(10, 6))

    # Plot execution time vs. number of particles for each function
    for function in pivot_df.columns:
        plt.plot(pivot_df.index, pivot_df[function], marker='o', label=function)

    plt.xlabel('Number of Particles')
    plt.ylabel('Execution Time (s)')
    plt.title('Execution Time vs. Number of Particles')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # Save the plot to a file
    plt.savefig('execution_time_vs_particles.png')

    # Show the plot
    plt.show()

    # Optional: Plot total execution time if available
    if 'total_t' in pivot_df.columns:
        plt.figure(figsize=(10, 6))
        plt.plot(pivot_df.index, pivot_df['total_t'], marker='o', label='Total Time')
        plt.xlabel('Number of Particles')
        plt.ylabel('Total Execution Time (s)')
        plt.title('Total Execution Time vs. Number of Particles')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig('total_execution_time_vs_particles.png')
        plt.show()

    # Optionally, you can plot the execution time on a log-log scale
    plt.figure(figsize=(10, 6))
    for function in pivot_df.columns:
        plt.loglog(pivot_df.index, pivot_df[function], marker='o', label=function)

    plt.xlabel('Number of Particles (log scale)')
    plt.ylabel('Execution Time (s, log scale)')
    plt.title('Execution Time vs. Number of Particles (Log-Log Scale)')
    plt.legend()
    plt.grid(True, which='both', linestyle='--')
    plt.tight_layout()
    plt.savefig('execution_time_vs_particles_loglog.png')
    plt.show()

if __name__ == '__main__':
    main()

