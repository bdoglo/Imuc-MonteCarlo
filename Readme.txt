This code is used to perform Monte-Carlo simulations on the current Phase II drug trial for ICT-107.

There are many assumptions in this code which could be wrong, and it is possible that the code has mistakes or bugs.
USE AT YOUR OWN RISK.  There are no warranties provided either explicitly or implicitly.  The code is provided as-is.
Please attribute the use of this code or any derivatives to the user bdoglo on Github.
---------------------------------------------------------------------
There are two matlab files here.

MCImucDay.m     main matlab function to call that performs multiple clinical trial simulations - 
                This function sets up the survival curves and other inputs to the simulator MCSimDay
                
MCSimDay.m      function called by MCImucDay to simulate the multiple trials
---------------------------------------------------------------------
This runs under Matlab and Octave.

If you use Windows. Here is the link to download Octave: 
http://sourceforge.net/projects/octave/files/Octave%20Windows%20binaries/Octave%203.6.1%20for%20Windows%20Microsoft%20Visual%20Studio/

Download the file "octave-3.6.1-vs2010-setup-1.exe" then follow the installer directions.

I haven't tested on that version of Octave as I have an older version, but it should work.
---------------------------------------------------------------------
Example usage:
MCImucDay(18.8+[0 9],1000,32);

This example runs 1000 trials drawn from two arms with survival times of 18.8 months and 27.8 months. 
The event count to detect is 32.  Thus the program provides samples from the 32nd event time distribution.

This shows the results of 32 events from two arms of 18.8 months OS and 27.8 months of overall survival.  The monte carlo is run for 1000 trial simulations.  In this example you should get something close to:
    mean = 834.608  (27 mo 13 days)
    median = 831.5  (27 mo 10 days)
    std-dev = 52.5107  (1 mo 22 days)

The times provided should be relative to Feb 1, 2011.  This means the average 32nd event would occur on May 13, 2013.  Note, that the results you get will vary from run to run.  To increase the accuracy you need to run more simulations which takes longer.  But 1000 simulations should get you within a few days or the actual mean.

Also note that to keep it simple each month, is assumed to be 365.25/12 days long.  
