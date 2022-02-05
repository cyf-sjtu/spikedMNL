universe = vanilla
getenv   = true
executable  = /opt/matlab/bin/matlab
arguments   = -singleCompThread -nojvm -nodisplay -r testing_compare_condor($(Process))
log         = simulation.log
output      = simulation.out
queue 1001

