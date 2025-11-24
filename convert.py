# swap_values.py

# Input and output file names
input_file = "prbs7_b.txt"
output_file = "prbs7_b_inv.txt"

# Read the input file
with open(input_file, "r") as f:
    lines = f.readlines()

# Process each line
output_lines = []
for line in lines:
    parts = line.strip().split()
    if len(parts) == 2:
        time, value = parts
        if value == "0.75":
            value = "0"
        elif value == "0":
            value = "0.75"
        output_lines.append(f"{time} {value}\n")

# Write to output file
with open(output_file, "w") as f:
    f.writelines(output_lines)

print(f"Converted file saved as {output_file}")
