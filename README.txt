********************************
Instructions for using toruspack
********************************


************************************
compile:

gcc -O2 toruspack.c -lm -o toruspack
************************************


*****************************************************************************************
run (without arguments):

./toruspack

expected 11 arguments: dim num size beta gam monostep itermax errstride errstop trials id


Description of arguments

dim:        (int) number of dimensions
num:        (int) number of spheres
size:       (double) width of torus
beta:       (double) RRR time-step
gam:        (double) RRR metric-update parameter
monostep:   (int) error monotonicity parameter
itermax:    (int) maximum number of RRR iterations per trial
errstride:  (int) number of iterations between printing output
errstop:    (double) stop iterations when error falls below this
trials:     (int) number of trials
id:         (char) string identifier for output files


Command line examples:

fcc.cmd     finds the fcc packing in 3D
best.cmd    finds Best's packing in 10D
14_100.cmd  finds packings of 115 spheres in 14D that have 100-monotone solutions
*****************************************************************************************


*************************************************************
Output files

id.cmd: command line of the run
id.run: summary of the run
id.err: error time series
id.sol: coordinates of sphere centers for the solutions found
id.wgt: time series of maximum metric weight
*************************************************************


For additional details, see arXiv:2305.13492