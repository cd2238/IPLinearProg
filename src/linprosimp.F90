!> Linear programming subroutine
!> from Nocedal Wright
!>
!> @param n the number of control variables
!> @param c the linear objective
!> @param m the number of linear constraints
!> @param A the matrix constraint
!> @param b the vector constraint
!> @param itermax the maximum number of iterations
!> @param mutol duality measure threshold
!> @param x the  primal solution
!> @param lambda the dual solution
!> @param s the  dual solution
!> @param iter the number of iterations
!> @param info information code

!> @todo propose other inputs possibilities : linear inequalities, etc.
!-----------------------------------------------------------------------
subroutine linprosimp ( n, c, m, A, b, itermax, mutol, x, lambda, s, iter,  info )
implicit none
!inputs
integer                         , intent(in)  :: n
double precision, dimension(n)  , intent(in)  :: c
integer                         , intent(in)  :: m
double precision, dimension(m,n), intent(in)  :: A
double precision, dimension(m)  , intent(in)  :: b
integer                         , intent(in)  :: itermax
double precision                , intent(in)  :: mutol
!outputs
double precision, dimension(n)  , intent(out) :: x
double precision, dimension(m)  , intent(out) :: lambda
double precision, dimension(m)  , intent(out) :: s
integer                         , intent(out) :: iter
integer                         , intent(out) :: info

! local variables
double precision epsi, mu
#ifdef DEBUG
integer i, j
#endif
double precision, dimension(:), allocatable :: x0
double precision, dimension(:), allocatable :: lambda0
double precision, dimension(:), allocatable :: s0

parameter(epsi=1.0d-15)
external initlp, linprosimpx

! init
info = 0

#ifdef DEBUG
  print*, "c=", (c(i), i=1,n)
  print*, "A =", ((A(i,j), j=1,n), char(10), i=1,m)
  print*, "b =",  (b(i), i=1,m)
#endif

! allocations
allocate(x0(n))
allocate(lambda0(m))
allocate(s0(n))

! initialization
call initlp( m, n, A, b, c, x0, lambda0, s0, info )
IF (info .LT. 0) THEN
  return
endIF


#ifdef DEBUG
  print*, "-><-"
  print*, "x0=", (x0(i), i=1,n)
  print*, "obj0 =", dot_product(c,x0)
  print*, "Ax0-b=", matmul(A,x0)-b
  print*, "-><-"
#endif


! linear program
call linprosimpx ( n, x0, lambda0, s0, c, m, A, b, itermax, mutol, &
                   x, lambda, s, mu, iter, info  )

#ifdef DEBUG
  print*, "-><-"
  print*, "results :"
  print*, "iter =", iter
  print*, "x =", (x(i), i=1,n)
  print*, "obj =", dot_product(c,x)
  print*, "info =", info
  print*, "-><-"
#endif

return
end


