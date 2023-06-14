import utils
import subprocess
import os
import matplotlib.pyplot as plt
from operator import itemgetter

script_path = "run.sh"
versions = ["plain", "openmp", "cuda", "combined"]
instances_directory = "./instances"

instance_paths = []
for root, dirs, files in os.walk(instances_directory):
    for file in files:
        file_path = os.path.join(root, file)
        instance_paths.append(file_path)
instance_paths = sorted(instance_paths, key=utils.extract_number_from_filename)

arg1 = "plain"
arg2 = "./instances/5.tsp"

results = list()

for instance_path in instance_paths:
    for version in versions:
        # Execute the Bash script and capture the output
        console_result = subprocess.run(["bash", script_path, version, instance_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Get the output
        output_data = {"Version": version, "Instance": instance_path.split('/')[2]}
        output = console_result.stdout.decode("utf-8")
        for line in output.split('\n'):
            if len(line) > 3:
                k, v = line.split(': ')
                if k == "Time":
                    output_data[k] = float(v.rstrip())

        results.append(output_data)
        print(output_data)

print()
print(results)

# Prepare the data for plotting
versions = []
instances = []
times = []
colors = []

for entry in results:
    versions.append(entry['Version'])
    instances.append(entry['Instance'])
    times.append(entry['Time'])
    if entry['Version'] == "plain":
        colors.append("black")
    if entry['Version'] == "openmp":
        colors.append("blue")
    if entry['Version'] == "cuda":
        colors.append("green")
    if entry['Version'] == "combined":
        colors.append("red")

# Create a figure and axis for the bar graph
fig, ax = plt.subplots()

# Set the x-axis tick positions and labels
x_pos = range(len(results))
ax.set_xticks(x_pos)
ax.set_xticklabels(versions)

# Plot the bars
ax.bar(x_pos, times, align='center', color=colors)

# Set the labels and title
ax.set_xlabel('Version')
ax.set_ylabel('Time')
ax.set_title('Time Comparison')

# Add the instance names as annotations above the bars
for i, instance in enumerate(instances):
    ax.text(i, times[i], instance, ha='center', va='bottom')

# Save the graph
plt.savefig('bar_graph.png')

# Display the bar graph
plt.show()

## Logarithmic version

# Create a new figure and axis for the second plot
fig2, ax2 = plt.subplots()

# Set the x-axis tick positions and labels
ax2.set_xticks(x_pos)
ax2.set_xticklabels(versions)

# Plot the bars with logarithmic scale
ax2.bar(x_pos, times, align='center', color=colors)
ax2.set_yscale('log')  # Set y-axis scale to logarithmic

# Set the labels and title for the second plot
ax2.set_xlabel('Version')
ax2.set_ylabel('Time (log scale)')
ax2.set_title('Time Comparison (Logarithmic Scale)')

# Add the instance names as annotations above the bars
for i, instance in enumerate(instances):
    ax2.text(i, times[i], instance, ha='center', va='bottom')

# Save the second plot as a PNG image
plt.savefig('bar_graph_log_scale.png')

# Display the second plot
plt.show()
