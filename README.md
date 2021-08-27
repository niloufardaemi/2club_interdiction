### Code for interdicting the maximum s-club in a graph

This code accompanies the paper "Interdicting Low-Diameter Cohesive Subgroups in
Large-Scale Social Networks" and is written in C++. If you wish to use or cite this code, please cite the paper:


@article{Niloufar2020InterClubs,
	author = {Niloufar Daemi and Juan S. Borrero and Balabhaskar Balasundaram},
	journal = {INFORMS Journal on Optimization},
	month = {August},
	note = {Accepted for publication.},
	title = {Interdicting low-diameter cohesive subgroups in large-scale social networks},
	year = {2021}
	}
  
  
This repository includes two folders:

1. exact_separation: used to solve the separation problem (maximum s-club) to optimality.
2. heuristic_separation: used to first try heuristic approaches to solve the separation problem and solving it to optimality only if needed.



## Compiling the code

The following steps show how to compile and run the code to find r-robust s-clubs in a Linux environment using a makefile (you can also run the code in Mac or Windows environment by configuring your IDE appropriately). In the folder of r_robust_s_club, parameter.txt is used to configure parameters r and s, and the "data" folder includes InputFile.txt (used to determine which instances you want to test) and 10th DIMACS graph instances (you can also downlaod these instances from the website: https://www.cc.gatech.edu/dimacs10/archive/clustering.shtml). It is similar to run the code to find h-hereditary s-clubs.

## Steps to run the code to find r-robust s-clubs in Linux environment:

Download or clone the repository to your machine.
Go to the folder r_robust_s_club.
Open the "Makefile" and set GUROBI_HOME to the directory of your Gurobi installation, e.g.: /opt/gurobi/9.0.1/linux64.
From the terminal, go to the folder of r_robust_s_club.
Type "make" and hit enter to compile. After successful complilation, type "./main" to run the code.


### Terms and Use:

* Describe any prerequisites, libraries, OS version, etc., needed before installing program.
* ex. Windows 10



### Acknowledgments
