import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("You should call this script with just one argument. This contains the name of the main .gol file (the one without underscores).")

#Open main file and get metadata
main_file_name, _ = os.path.splitext(sys.argv[1])
main_file = open(main_file_name + ".gol")
data_str = main_file.readline()
(x_dim, y_dim, iter_step, iter_end, processes) = tuple(map(int, data_str.split(" ")))
main_file.close()

#Iterate through files and visualize
for iteration in range(0, iter_end+1, iter_step):
    board = np.zeros((x_dim, y_dim), dtype=np.bool)
    for process in range(processes):
        file = main_file_name + "_" + str(iteration) + "_" + str(process) + ".gol"
        print(file)
        part_file = open(file)
        #Get metadata
        x_range_str = part_file.readline()
        (x_min, x_max) = tuple(map(int, x_range_str.split(" ")))
        y_range_str = part_file.readline()
        (y_min, y_max) = tuple(map(int, y_range_str.split(" ")))
        if x_min < x_max and y_min < y_max:
            #Get data from file
            part_board = np.array([line.split() for line in part_file])
            #Write data to array
            board[x_min:x_max+1, y_min:y_max+1] = part_board
        part_file.close()
    #Visualize
    plt.pcolor(board)
    plt.title("Iteration=" + str(iteration))
    plt.show(block=False)
    plt.pause(0.5)