This code is used to perform Monte-Carlo simulations on the current Phase II drug trial for ICT-107.

There are many assumptions in this code which could be wrong, and it is possible that the code has mistakes or bugs.
USE AT YOUR OWN RISK.  There are no warranties provided either explicitly or implicitly.  The code is provided as-is.
Please attribute the use of this code or any derivatives to the user bdoglo on Github.

There are two matlab files here.

MCImucDay.m     main matlab function to call that performs multiple clinical trial simulations - 
                This function sets up the survival curves and other inputs to the simulator MCSimDay
                
MCSimDay.m      function called by MCImucDay to simulate the multiple trials

Example usage:
MCImucDay(18.8+[0 9],5000,64);

This example runs 5000 trials drawn from two arms with survival times of 18.8 months and 27.8 months. 
The event count to detect is 64.  Thus the program provides samples from the 64th event time distribution.
