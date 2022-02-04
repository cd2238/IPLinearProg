!> Provides "not too bad" initial values for primal and dual points
!>
!> @param m the number of linear constraints
!> @param n the number of control variables
!> @param A the matrix constraint
!> @param b the vector constraint
!> @param c the linear objective
!> @param x the  primal solution
!> @param lambda the dual solution
!> @param s the  dual solution
!> @param info information code
!-----------------------------------------------------------------------
subroutine initlp ( m, n, A, b, c, x, lambda, s, info )
implicit none
!inputs
integer                         , intent(in)  :: m
integer                         , intent(in)  :: n
double precision, dimension(m,n), intent(in)  :: A
double precision, dimension(m)  , intent(in)  :: b
double precision, dimension(n)  , intent(in)  :: c

!outputs
double precision, dimension(n)  , intent(out) :: x
double precision, dimension(m)  , intent(out) :: lambda
double precision, dimension(m)  , intent(out) :: s
integer                         , intent(out) :: info

! local variables
integer allocock
double precision ps, sx, ss, delta_s, delta_x, epsi
integer lwork

integer nrhs, mn

double precision, dimension(:,:), allocatable :: AA
double precision, dimension(:)  , allocatable :: bb
double precision, dimension(:)  , allocatable :: cc
double precision, dimension(:)  , allocatable :: dwork

external dgels
parameter(epsi=1.0d-15)

mn = max(m,n)


! allocations
allocate(AA(m,n), stat = allocock)
if (allocock /= 0) return    
allocate(bb(mn), stat = allocock)
if (allocock /= 0) return    
allocate(cc(mn), stat = allocock)
if (allocock /= 0) return    

! initializations
info = 0
nrhs = 1
AA = A
bb = 0.0d0
cc = 0.0d0
bb(1:m) = b
cc(1:n) = c



!! x (undetermined system)
lwork = 2*n
allocate(dwork(lwork))
call dgels('N', m, n, nrhs, AA, m, bb, mn, dwork, lwork, info)
IF (info .LT. 0) THEN
  return
ENDIF
x = bb

!! lambda (overdetermined system)
call DGELS('T', m, n, nrhs, AA, m, cc, mn, dwork, lwork, info)
IF (info .LT. 0) THEN
  print*, "info=", info
  return
ENDIF
lambda= cc

!! s
s = matmul(transpose(A), lambda)
s = c - s


! deltax,s
delta_x = max(-1.5d0*minval(x), 0.0d0)
delta_s = max(-1.5d0*minval(s), 0.0d0)

! x,s
x = x + delta_x
s = s + delta_s

! hatdeltax,s
ps = dot_product(x, s)

sx= sum(x)
ss =sum(s)
IF (abs(sx) .GT. epsi) THEN
  delta_x = 0.5d0 * ps / sx
ELSE
  delta_x = 0.0d0
ENDIF

IF (abs(ss) .GT. epsi) THEN
  delta_s = 0.5d0 * ps / ss
ELSE
  delta_s = 0.0d0
ENDIF

! x,s
x = x + delta_x
s = s + delta_s




return
end


