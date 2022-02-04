!> Linear programming subroutine
!> from Nocedal Wright
!>
!> @param n the number of control variables
!> @param c the linear objective
!> @param m the number of linear constraints
!> @param A the matrix constraint
!> @param b the vector constraint
!> @param itermax the maximum number of iterations
!> @param x the  primal solution
!> @param lambda the dual solution
!> @param s the  dual solution
!> @param mu the  duality measure
!> @param iter the number of iterations
!> @param info information code

!> @todo propose other inputs possibilities : linear equalities, linear inequalities, etc.
!-----------------------------------------------------------------------
SUBROUTINE linprosimp ( n, c, m, A, b, itermax, x, lambda, s, iter,  info )
!-----------------------------------------------------------------------
IMPLICIT NONE

INTEGER          :: n         !> @info information
DOUBLE PRECISION :: c(n)      !> @c linear objective
INTEGER          :: m         !> @m size
DOUBLE PRECISION :: A(m, n)   !> @A constraint matrix
DOUBLE PRECISION :: b(m)      !> @b vector matrix
INTEGER          :: itermax   !> @itermax itermax

DOUBLE PRECISION :: x(n)      !> @x primal solution
DOUBLE PRECISION :: lambda(m) !> @lambda dual solution
DOUBLE PRECISION :: s(n)      !> @s slack solution
INTEGER          :: iter      !> @iter iter
INTEGER          :: info      !> @info information

! x calculous csup, x(n)


DOUBLE PRECISION epsi, mu
#ifdef DEBUG
INTEGER i, j
#endif

DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x0
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambda0
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: s0

PARAMETER(epsi=1.0d-15)
EXTERNAL initlp, linprosimpx


info = 0

#ifdef DEBUG
  print*, "c=", (c(i), i=1,n)
  print*, "A =", ((A(i,j), j=1,n), char(10), i=1,m)
  print*, "b =",  (b(i), i=1,m)
#endif


allocate(x0(n))
allocate(lambda0(m))
allocate(s0(n))

print*, "coucou"
! initialization
CALL initlp( m, n, A, b, c, x0, lambda0, s0, info )
IF (info .LT. 0) THEN
  RETURN
ENDIF


#ifdef DEBUG
  print*, "-><-"
  print*, "x0=", (x0(i), i=1,n)
  print*, "obj0 =", dot_product(c,x0)
  print*, "Ax0-b=", matmul(A,x0)-b
  print*, "-><-"
#endif


! algorithm
CALL linprosimpx ( n, x0, lambda0, s0, c, m, A, b, itermax, &
                   x, lambda, s, mu, iter, info  )

#ifdef DEBUG
  print*, "iter =", iter
  print*, "x =", (x(i), i=1,n)
  print*, "obj =", dot_product(c,x)
  print*, "info =", info
#endif

RETURN
END


