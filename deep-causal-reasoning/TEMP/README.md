This is a branch of the files Lekshmi is working on. 

Contains: 
1. myTest: Testing Ioannis' code on a toy network. Script to be run in python is 'mystudyfile.py'. The relevant result files are:
ds1.txt: distances
Ums1.txt and Ups1.txt : Um and Ups
xs1.txt: contains X 
us1.txt: contains the relevant reactions with Us
The constraints are printed on the main.lp file. 
Output log is outlog. 

2. myMatrix: Testing the sympy based code on the toy network and data. Relevant results to be looked with the results in the myTest folder are: 
ds_matrix.txt with myTest/ds1.txt
us_matrix.txt with myTest/us1.txt
xs_matrix.txt with myTest/xs.txt
The constraints are printed on the main.lp file. (compare with myTest/main.lp)
output log is out_matrix.txt (compare with myTest/outlog)

3. MKN1: runs the data from MKN1 cell line using both the codes. The input is the same for both files. 
a) Ioannis' code (signettrainer_subgraph.py)
INPUT: The options file. ilp_options. (The inputs file has a redundant input information removed.. KRAS has more than one entry where in some it is 1 and in one of them is is NAN. The NAN is removed. Ioannis' code handles this. But Lekshmi passes everything as matrices)
      The network (stitch_interactions_signed.sif), measurements and input files.
OUTPUT: The relevant out puts are stored in us2.txt, xs2.txt. 
logfile: outlog2

b)The matrix form of writing the constraints is run using sympybased_test2.py (generates constarints in the same order as Ioannis' code) or sympybased_test.py (Generates constraints in a different order)
INPUT: The inputs to the code are matrices :
Inputsf.mtx - inputs
Rf.mtx: reactant matrix
Pf.mtx: Product matrix
Mf.mtx: measurements
Mindf.mtx: M_ind matrix 

Detailed info about what these mean is in the Readme folder ilp_form.pdf. 
Output logs are out_matrix1 by running sympybased_test.py and out_matrix2.txt by running sympybased_test2.py

