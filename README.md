# Linear programming with interior point methods


This program implements the interior point linear programming proposed by Nocedal and Wright. It uses a robust Choleski factorization described in Gill et al. It has been witten in order to better understand the interior point theory. It is not an operational solver !




## References :
* Numerical optimization, Nocedal, Wright, 2nd edition, Springer, 2000
* Primal-dual interior-point methods, Wright, 1997
* Practical optimization, Gill, Murray, Wright, Academic press, 1984

## Install
* sudo apt install gfortran liblapack-dev doxygen

## compilation
* cd make
* make all
* make all STATUS=debug (for some intermediate results and gdb possibilities)

## launch
* cd bin
* linprog n10c04  (folder name of the data set in the data folder)

## To do
* Add a more usual formulation for linear programming
* Add sparse matrices management
* Add preprocessing (delete dependent rows,...)
