This program solves inviscid Burgers equation using FVM for 3 different initial conditions that can be chosen as the program runs. 4 different flux functions are available. These must be chosen before executing the program by making changes in the source code itself.



Execute the program:


1. Using code-blocks

Open Project-1.cbp using code-blocks. Build and run the project.



2. Using terminal
Compile main.c, "gcc main.c -lm -o out"



During the execution:


1. Choose the number of divisions (any integer, say 100)

2. Choose the CFL number (less than 1, say 0.4)

3. Choose the initial condition - presently coded for jumps between 0 and 1 and a peak initial condition. Initial conditions can be changed to jumps between -1 and 1. The exact solutions will have to be changed. The rest of the code works fine. 



After execution:


Apart from object and exe files, the program generates 'burgersexact.dat'  and 'burgersnum.dat' .  These files contain the numerical and exact solutions to the equation solved. They can be plotted using the 'plotting-script' which generates the file 'burgers.pdf' or by manually ececuting the commands in that script.
 
The peak initial condition doesn't generate the burgersexact.dat file so the plot only contains the numerical solution. 
