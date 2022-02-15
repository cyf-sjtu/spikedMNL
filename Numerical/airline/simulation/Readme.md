Run the following steps:
Simulation.m
testing_condor.m

Use condor (see execute.cmd) to run testing_condor.m to generate results for Section 8.2.3

************************************************************************************************

airline_data.mat stores the calibrated demand model (arrival process + discrete choice model).

Flight.m is the class defined to store flight info and specify flight relatived activities.

SBLP_with_spike.m SBLP_without_spike.m solved the SBLP problem tailored for the airline setting.

choose_product.m simulates customer purchase behavior.

evaluate_XXX.m improves control policies using SGD.

find_flight_product.m maps the internal product index for simulation use.

generate_xxx.m generates customer arrival processes and demands.

get_EMSRb.m computes the EMSRb booking limits.
