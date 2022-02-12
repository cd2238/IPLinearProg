!> resolves 1 linear system (affine-scaling step) for interior point linear programming, and returns also intermediary results

!> If it evaluates to .false.,
!> print the fail message and possibly stop execution. (A more detailed description)
!>
!> @param m the number of linear constraints
!> @param n the number of control variables
!> @param A the matrix constraint
!> @param s the matrix constraint
!
!     Solve the system :
!       [ 0  A' I ] [deltax     ]   [-rc ]
!       [ A  0  0 ] [deltalambda] = [-rb ]
!       [ S  0  X ] [deltas     ]   [-rxs]
!   with :
!          S = diag(s)
!          I = Id_n
!          X = diag(x)
!
!-----------------------------------------------------------------------
subroutine solvesyslp ( m, n, A, s, x, rc, rb, rxs, deltax, deltalambda, deltas, AD2A, A3, xds, t2, info )
!-----------------------------------------------------------------------
implicit none
integer info, m, n
double precision A(m,n), s(n), x(n), rc(n), rb(m), rxs(n), xds(n)
double precision deltax(n), deltalambda(m), deltas(n), AD2A(m,m), t2(m), A3(m,n)


integer i,j, mn
double precision prec, alpha, delta
integer cho

external dpotrs, solvetriglp, modchol2


double precision, dimension(:),   allocatable :: b
double precision, dimension(:,:), allocatable :: d2
double precision, dimension(:,:), allocatable :: AD2
double precision, dimension(:),   allocatable :: b2
double precision, dimension(:),   allocatable :: b3
double precision, dimension(:),   allocatable :: u
double precision, dimension(:),   allocatable :: v

parameter(prec=1.0d-12, alpha = 1.0d-8, cho=2, delta=1.0d-10)
!cho = 0 : cholesky lapack
!cho = 2 : modified cholesky


mn=max(m,n)

allocate(b(m))
allocate(d2(n,n))
allocate(AD2(m,n))
allocate(b2(m))
allocate(b3(m))
allocate(u(n))
allocate(v(m))


info = 0

print*, "A=[[", ((A(i,j),j=1,n),char(10), i=1,m), "]]"
print*, "s=[", (s(i),i=1,n), "]"
print*, "x=[", (x(i),i=1,n), "]"
print*, "rc=[", (rc(i),i=1,n), "]'"
print*, "rb=[", (rb(i),i=1,m), "]'"
print*, "rxs=[", (rxs(i),i=1,n), "]'"

! calculates AD^2A'
! AD^2
xds = x/s ! xs^-1
DO j=1,n
  AD2(:,j) = A(:,j)*xds(j)
  A3(:,j)  = A(:,j)/s(j)
ENDDO
!print*, "AD2=[[", ((dwork(pdA2+(j-1)*m+i-1),j=1,n),char(10), i=1,m), "]]"
AD2A = matmul(AD2, transpose(A))
print*, "AD2A=[[", ((AD2A(i,j),j=1,m),char(10), i=1,m), "]]"

! calculates deltalambda(val provisoire) = -rb-AXS^-1rc
t2 = -rb-matmul(AD2,rc)

print*, "termb=[",  (t2(i),i=1,m), "]'"

! solves AD^2A' deltalambda = -rb-AXS^-1rc+AS^-1rxs
IF (cho == 1) THEN !  chol
  CALL dpotrf('L', m, AD2A, m, info) 
ELSE ! cho==2
  ! modified choleski
  CALL modchol2( m, AD2A, delta, AD2A, info )
ENDIF
IF (info .NE. 0) THEN
  print*, "info-chol =", info
  RETURN
ENDIF

#ifdef DEBUG
PRINT*, "chol=[[", ((AD2A(i,j), j=1,m), char(10), i=1,m), "]]"
#endif

CALL SOLVETRIGLP ( m, n, A, AD2A, A3, s, xds, rc, rxs, t2, deltax, deltalambda, deltas, info )
IF (info .LT. 0) THEN
  RETURN
ENDIF

print*, "deltax=[", (deltax(i),i=1,n), "]'"
print*, "deltalambda=[", (deltalambda(i),i=1,m), "]'"
print*, "deltas=[", (deltas(i),i=1,n), "]'"

!check the large newton equation
u = matmul(transpose(A), deltalambda) + deltas
IF (sum(abs(u+rc)).GT. prec) THEN
  PRINT*, "prob 1:", sum(abs(u+rc))
  !STOP
ENDIF


v = matmul(A, deltax)
IF (sum(abs(v+rb)).GT. prec) THEN
  PRINT*, "prob 2:", sum(abs(v+rb))
  !STOP
ENDIF

u = s*deltax+x*deltas
IF (sum(abs(u+rxs)).GT. prec) THEN
  PRINT*, "prob 3:", sum(abs(u+rxs))
  !STOP
ENDIF


RETURN
END


