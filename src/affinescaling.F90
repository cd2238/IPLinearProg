!> resolves 1 linear system (affine-scaling step) for interior point linear programming, and returns also intermediary results
!> see Nocedal Wright (eq. 14.30)
!>
!> @param m the number of linear constraints
!> @param n the number of control variables
!> @param A the matrix constraint
!> @param s the current dual solution
!> @param x the current primal solution
!> @param rc see Nocedal-Wright
!> @param rb see Nocedal-Wright
!> @param rxs = XSe see Nocedal-Wright
!> @param AD2A A*D^2*A' (intermediary result to avoid redundancy)
!> @param A3 A*S^{-1} (intermediary result to avoid redundancy)
!> @param t2 -rb-AD^2rc (intermediary result to avoid redundancy)


!!
!!     Solve the system :
!!      [ 0  A' I ] [deltax     ]   [-rc ]
!!      [ A  0  0 ] [deltalambda] = [-rb ]
!!      [ S  0  X ] [deltas     ]   [-rxs]
!!  with :
!!         S = diag(s)
!!         I = Id_n
!!         X = diag(x)
!!         D^2 = X*S^{-1}
!!         rxs = XSe
!------------------------------------------------------------------------------------------------------------
subroutine affinescaling ( m, n, A, s, x, rc, rb, rxs, deltax, deltalambda, deltas, AD2A, A3, xds, t2, info )
implicit none

! inputs
integer                         , intent(in) :: m
integer                         , intent(in) :: n
double precision, dimension(m,n), intent(in) :: A
double precision, dimension(n),   intent(in) :: s
double precision, dimension(n),   intent(in) :: x
double precision, dimension(n),   intent(in) :: rc
double precision, dimension(m),   intent(in) :: rb
double precision, dimension(n),   intent(in) :: rxs

! outputs
double precision, dimension(n),   intent(out) :: deltax
double precision, dimension(m),   intent(out) :: deltalambda
double precision, dimension(n),   intent(out) :: deltas
double precision, dimension(m,m), intent(out) :: AD2A
double precision, dimension(m),   intent(out) :: t2
double precision, dimension(m,n), intent(out) :: A3
double precision, dimension(n),   intent(out) :: xds
integer                         , intent(out) :: info

! local variables
integer i,j, mn, cho
double precision prec, alpha

double precision, dimension(:),   allocatable :: b
double precision, dimension(:,:), allocatable :: d2
double precision, dimension(:,:), allocatable :: AD2
double precision, dimension(:),   allocatable :: b2
double precision, dimension(:),   allocatable :: b3
#ifdef DEBUG
double precision, dimension(:),   allocatable :: u
double precision, dimension(:),   allocatable :: v
#endif

external dpotrs, solvetriglp, modchol2
parameter(prec=1.0d-12, alpha = 1.0d-8, cho=1)
!cho = 0 : cholesky lapack
!cho = 1 : modified cholesky

! initialization
mn=max(m,n)


! allocations
allocate(b(m))
allocate(d2(n,n))
allocate(AD2(m,n))
allocate(b2(m))
allocate(b3(m))
#ifdef DEBUG
allocate(u(n))
allocate(v(m))
#endif

info = 0

#ifdef DEBUG
print*, "A=[[", ((A(i,j),j=1,n),char(10), i=1,m), "]]"
print*, "s=[", (s(i),i=1,n), "]"
print*, "x=[", (x(i),i=1,n), "]"
print*, "rc=[", (rc(i),i=1,n), "]'"
print*, "rb=[", (rb(i),i=1,m), "]'"
print*, "rxs=[", (rxs(i),i=1,n), "]'"
#endif

! calculates AD^2A', A3, t2
! AD^2
xds = x/s ! xs^-1
do j=1,n
  AD2(:,j) = A(:,j)*xds(j)
  A3(:,j)  = A(:,j)/s(j)
enddo
AD2A = matmul(AD2, transpose(A))
t2 = -rb-matmul(AD2,rc)



#ifdef DEBUG
print*, "AD2A=[[", ((AD2A(i,j),j=1,m),char(10), i=1,m), "]]"
print*, "t2=[",  (t2(i),i=1,m), "]'"
#endif

! cholesky factorization
if (cho == 0) then
  !  lapack cholesky
  call dpotrf('L', m, AD2A, m, info)
else ! cho==1
  ! modified cholesky
  call modchol2( m, AD2A, AD2A, info )
endif
if (info .NE. 0) then
  return
endif

#ifdef DEBUG
print*, "chol=[[", ((AD2A(i,j), j=1,m), char(10), i=1,m), "]]"
#endif

! performs centering/corrector step
call solvetriglp ( m, n, A, AD2A, A3, s, xds, rc, rxs, t2, deltax, deltalambda, deltas, info )
if (info .LT. 0) then
  return
endif

#ifdef DEBUG
print*, "deltax=[", (deltax(i),i=1,n), "]'"
print*, "deltalambda=[", (deltalambda(i),i=1,m), "]'"
print*, "deltas=[", (deltas(i),i=1,n), "]'"
#endif

#ifdef DEBUG
!check the large newton equation
u = matmul(transpose(A), deltalambda) + deltas
if (sum(abs(u+rc)).GT. prec) then
  print*, "prob 1:", sum(abs(u+rc))
  !STOP
endif


v = matmul(A, deltax)
if (sum(abs(v+rb)).GT. prec) then
  print*, "prob 2:", sum(abs(v+rb))
  !STOP
endif

u = s*deltax+x*deltas
if (sum(abs(u+rxs)).GT. prec) then
  print*, "prob 3:", sum(abs(u+rxs))
  !STOP
endif
#endif


return
end


