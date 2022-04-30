import os
import re

path = "./experiments"
minimum = float('inf')
minimum_data = None
top_parameters = []
  
def get_file_result(file_path, file):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        result_line = lines[19]
        result = float(result_line.split(' ')[1])

        data = '\n'
        for i in range(12):
            data += lines[i]

        data += '\n'
        parameters_line = 21
        for i in range(parameters_line, len(lines)):
            
            data += lines[i]
        
        top_parameters.append((result, data))

        return result, data

for file in os.listdir(path):
    file_path = f"{path}/{file}"

    result, data = get_file_result(file_path, file)

    if result < minimum:
        minimum = result
        minimum_data = data 

top_parameters.sort(key=lambda x: x[0])

print("Minimum: ", minimum)
print("Parameters: ", minimum_data)
print("Top 10: ")
for i in range(10):
    print("Value: ", top_parameters[i][0], " , Parameters: ", top_parameters[i][1])