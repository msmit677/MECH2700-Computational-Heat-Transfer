# MECH2700-Computational-Heat-Transfer
Project that models the heat/temperature distribution in a rod using numerical methods

Running the Script:

To run the script, simply click the run button on the Python file that was submitted and all results from all tasks should be produced after about one minute of computations. If you would like to change any parameters, this can be done at the top of the script which is clearly defined.

Functions:

-A(n,dx) = creates an (n+1)x(n+1) matrix of the form presented in the report, which is tri-diagonal in its structure. The matrix begins as a zero matrix and entries are added over 3 steps

-b(n,dx) = creates an (n+1)x1 matrix of the forn presented in the report, where the first and last entries were added to a zero matrix and then the middle entries were added

-Solve(n,dx) = runs the chosen solver method, which in this case is a modified version of the Gauss-Jordan elimination method

-T_Analytic(x) = creates the analytic temperature function presented in the report which depends on x

-T_Convergence() = A function which finds the minimum n value required to produce a 0.1% error of the numerical solution from the solver method, which is achieved through a while loop

-q(n) = the numerical heat transfer rate function which depends on the value of n. It runs the solve function to first produce a set of T values

-q_Analytic(L,Ac,P) = creates the analytic heat transfer rate equation presented in the report

-q_Convergence() = Similar to the T_Convergence function, this function finds the minumum n value required to produce a 0.1% error of the heat transfer numerical solution, through a while loop

-LChange(),dChange(),tChange() = functions which calculate the required parameter values (L,d,t) to achieve the desired heat transfer rate (14kW)
