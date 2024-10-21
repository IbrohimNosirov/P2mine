#!/usr/bin/env python3

import sys

def compare_simulations(file1, file2, output_file, tolerance=1e-8):
    """
    Compares the coordinates in two simulation output files line by line and computes the differences.

    Parameters:
    - file1: Path to the first simulation output file.
    - file2: Path to the second simulation output file.
    - output_file: Path to the output file where differences will be written.
    - tolerance: A float specifying the tolerance for considering differences as zero.
    """
    try:
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            # Read headers from both files
            header1 = f1.readline().strip()
            header2 = f2.readline().strip()

            # Optionally, compare headers to ensure they match
            if header1 != header2:
                print("Warning: Headers do not match.")
                print(f"File1 header: {header1}")
                print(f"File2 header: {header2}")
                # Decide whether to proceed or exit
                # sys.exit(1)

            # Initialize lists to store differences
            differences = []

            line_number = 2  # Start from line 2 since we've read the header
            while True:
                line1 = f1.readline()
                line2 = f2.readline()

                # Check for end of file
                if not line1 and not line2:
                    break  # Both files have been read completely
                elif not line1 or not line2:
                    print(f"Error: Files have different number of lines at line {line_number}.")
                    break

                # Strip whitespace and split into components
                tokens1 = line1.strip().split()
                tokens2 = line2.strip().split()

                # Skip empty lines (if any)
                if not tokens1 and not tokens2:
                    line_number += 1
                    continue
                elif not tokens1 or not tokens2:
                    print(f"Error: Mismatched lines at line {line_number}.")
                    line_number += 1
                    continue

                # Convert tokens to floats
                try:
                    coords1 = [float(token) for token in tokens1]
                    coords2 = [float(token) for token in tokens2]
                except ValueError as e:
                    print(f"Error parsing line {line_number}: {e}")
                    line_number += 1
                    continue

                # Ensure both lines have the same number of coordinates
                if len(coords1) != len(coords2):
                    print(f"Error: Different number of coordinates at line {line_number}.")
                    line_number += 1
                    continue

                # Compute differences
                diff = [c1 - c2 for c1, c2 in zip(coords1, coords2)]

                # Append differences to the list
                differences.append(diff)

                line_number += 1

        # Write differences to the output file
        with open(output_file, 'w') as fout:
            fout.write("Differences between simulations\n")
            fout.write("Format: dx dy dz (line numbers correspond to original files)\n\n")
            for idx, diff in enumerate(differences, start=2):  # Starting from line 2
                diff_line = ' '.join(f"{d:.6e}" for d in diff)
                # Optionally check if differences are within tolerance
                if all(abs(d) <= tolerance for d in diff):
                    pass  # Differences are within tolerance, can skip or write zeros
                    fout.write(f"Line {idx}: 0.000000e+00 0.000000e+00 0.000000e+00\n")
                else:
                    fout.write(f"Line {idx}: {diff_line}\n")

        print(f"Differences written to {output_file}.")

    except Exception as e:
        print(f"An error occurred: {e}")

# Example usage when running the script directly
if __name__ == "__main__":
    # Check for correct number of arguments
    if len(sys.argv) != 4:
        print("Usage: python3 compare_simulations.py <file1> <file2> <output_file>")
        sys.exit(1)

    # Get file paths from command-line arguments
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3]

    # Call the function to compare simulations
    compare_simulations(file1, file2, output_file)
