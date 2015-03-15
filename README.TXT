To compile this program, type make
To run this program. run ./TSP followed by the 4 command line arguments in the corresponding order
	first command line argument:  W = width
	second command line argument: H = height
	third command line argument:  N = number of vertices
	fourth command line argument: T = number of trials


The program will display unsupressed output from the provided PerfectMatching code while it runs the trials. At the end, it will display a chart that will look something like this: 

Number of Trials: 50
___________________________
MST Weight --> Mean: 104039  Std: 818.684
MST Time   --> Mean: 0.00599676  Std: 0.000798796
TSP-2 Weight --> Mean: 158925  Std: 1885.5
TSP-2 Time -->  Mean: 0.00123022  Std: 0.000187853
TSP-1.5 Weight --> Mean: 128460 Std: 1327.02
TSP-1.5 Time --> Mean: 0.0651483 Std: 0.00765194

This displays relevant information about the calculation time and the weight of the MST, TSP-2 Approximation and TSP-1.5 Approximation.


