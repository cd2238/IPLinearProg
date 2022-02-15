# Linear programming with interior point methods


This program implements the interior point linear programming proposed by Nocedal and Wright (derived from Mehrotra). It uses a robust Choleski factorization described in Gill et al. It has been witten in order to better understand the interior point theory. It has 3 steps, predictor step, corrector step and centering step. It is not an operational solver !

 min c'x  
 sc   
    | Ax = b  
    | x >= 0  


## References :
* Numerical optimization, Nocedal, Wright, 2nd edition, Springer, 2000
* On the implementation of a primal-dual interior point method, Mehrotra, SIAM journal on optimization, 1992.
* Primal-dual interior-point methods, Wright, 1997
* Practical optimization, Gill, Murray, Wright, Academic press, 1984

## Install (ubuntu 18.04)
* sudo apt install gfortran liblapack-dev doxygen

## compilation
* cd make
* make all
* make all STATUS=debug (for some intermediate results and gdb possibilities)

## run
* cd bin
* linprog n10c04  (folder name of the data set in the data folder)

## To do
* Add a more usual formulation for linear programming
* Add sparse matrices management
* Add preprocessing (delete dependent rows,...)
