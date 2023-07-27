# wrapper take N
# run a loop for i from 3..N
# give i as input to generate.py
# run generate.py for j = 0..2
# for every run in j, calculate time taken by dcel1.cpp
# output that into a file
# add all of the three times, take avg
# write this time against i in separate file
import os
import sys
import matplotlib.pyplot as plt


num_test = int(input("Enter the number of max vertices: "))
open("time_output.txt", "w+").close()
for i in range(3, num_test+1):
    open("times.txt", "w+").close()
    for j in range(3):
        os.system("python3 generate.py " + str(i))
    with open('times.txt', 'r') as input_file, open('time_output.txt', 'a') as output_file:
        numbers = [float(number) for number in input_file.readline().split()]
        average = sum(numbers) / len(numbers)
        output_file.write(str(i) + " " + str(average) + "\n")

with open("time_output.txt", "r") as f:
    lines = f.readlines()

# Take input from the time taken vs number of vertices file (time_output.txt)
vertices = []
time_taken = []
for line in lines:
    v, t = line.strip().split(" ")
    vertices.append(int(v))
    time_taken.append(float(t)*1000)

# Plot time taken vs number of vertices 
plt.plot(vertices, time_taken, color = "black", linewidth = 1.5)
plt.title("avg. Time Taken vs Number of vertices", color = "Purple", style = "italic", size = 15)
plt.xlabel("Number of Vertices", size = 10, color = 'red')
plt.ylabel("Time Taken (ms)", size = 10, color = 'red')
plt.show()
