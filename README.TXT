To compile this program, type make
To run this program. run ./TSP followed by the 5 command line arguments in the corresponding order
	first command line argument:  W = width
	second command line argument: H = height
	third command line argument:  N = number of vertices
	fourth command line argument: T = number of trials
	fifth command line argument:  E = whether or not to run Extra Credit (1 = run, 0 = do not run)

**NOTE**: For larger N values (10,000), expect very long times for EC to run. Even for 1 trial, it can take 25 minutes each. 

The program will display unsupressed output from the provided PerfectMatching code while it runs the trials. At the end, it will display a chart that will look something like this: 

W = 10000 H = 100000 N = 100000 Trials = 100
___________________________
MST  --> Mean: 104039  Std: 818.684
TSP-2  --> Mean: 158925  Std: 1885.5
TSP-1.5  --> Mean: 128460 Std: 1327.02
ECE 20PT-E  --> Mean: 119453 Std: 1287.09




