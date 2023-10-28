The project is an implementation and analysis of https://www.sciencedirect.com/science/article/pii/S0377221799000338#:~:text=They%20start%20with%20the%20whole,last%20polygon%20of%20the%20partition.https://github.com/arch1tek/Polygon-Decomposition/blob/main/readme.txt


To run the program with a randomly generated input with number of vertices  as user input : 

>> open the project directory, and then type in the command :
    python generate.py <number_of_vertices>

>> the output is stored in the directory Outputfiles as a png.

To run the program with a given already created input :

>> create a file named "generated_test2.txt" and place your vertices in Clockwise order
>> run the command in the correct directory as follows :
    g++ dcel1.cpp
    ./a.out <number_of_vertices>
>> the output is stored in the directory Outputfiles as a png.

The doxygen file can be created by using the following command : 
     doxygen -g <config_file>
>> necessary changes can be made in config_file
     doxygen 
>> run the index.html file from html directory in doxygen_docs created, to view the html page locally.




The wrapper.py is a python script used for calculating the average time vs the input size

The files name "time_same_n.txt", "plot_answer.txt", "generate_test1.txt", "generate_test2.txt" are all auxilary files used for plotting analysis

The files timeplotter.py and timeplotter2.py were also used to plot time analysis graphs
They can be ignored, or used as per the requirements and will of the user.

