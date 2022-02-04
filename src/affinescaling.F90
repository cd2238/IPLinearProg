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


!>
!>     Solve the system :
!!      [ 0  A' I ] [deltax     ]   [-rc ]
!>      [ A  0  0 ] [deltalambda] = [-rb ]
!>      [ S  0  X ] [deltas     ]   [-rxs]
!>  with :
!>         S = diag(s)
!>         I = Id_n
!>         X = diag(x)
!>         D^2 = X*S^{-1}
!>         rxs = XSe
!-----------------------------------------------------------------------
subroutine affinescaling ( m, n, A, s, x, rc, rb, rxs, deltax, deltalambda, deltas, AD2A, A3, xds, t2, info )
!-----------------------------------------------------------------------
IMPLICIT NONE

INTEGER                         , INTENT(in) :: m
INTEGER                         , INTENT(in) :: n
DOUBLE PRECISION, DIMENSION(m,n), INTENT(in) :: A
DOUBLE PRECISION, DIMENSION(n),   INTENT(in) :: s
DOUBLE PRECISION, DIMENSION(n),   INTENT(in) :: x
DOUBLE PRECISION, DIMENSION(n),   INTENT(in) :: rc
DOUBLE PRECISION, DIMENSION(m),   INTENT(in) :: rb
DOUBLE PRECISION, DIMENSION(n),   INTENT(in) :: rxs

DOUBLE PRECISION, DIMENSION(n),   INTENT(out) :: deltax
DOUBLE PRECISION, DIMENSION(m),   INTENT(out) :: deltalambda
DOUBLE PRECISION, DIMENSION(n),   INTENT(out) :: deltas
DOUBLE PRECISION, DIMENSION(m,m), INTENT(out) :: AD2A
DOUBLE PRECISION, DIMENSION(m),   INTENT(out) :: t2
DOUBLE PRECISION, DIMENSION(m,n), INTENT(out) :: A3
DOUBLE PRECISION, DIMENSION(n),   INTENT(out) :: xds
INTEGER                         , INTENT(out) :: info

INTEGER i,j, mn, cho
DOUBLE PRECISION prec, alpha, delta

DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: b
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: d2
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AD2
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: b2
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: b3
#ifdef DEBUG
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: u
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: v
#endif

EXTERNAL dpotrs, solvetriglp, modchol2
PARAMETER(prec=1.0d-12, alpha = 1.0d-8, cho=1, delta=1.0d-10)
!cho = 0 : cholesky lapack
!cho = 1 : modified cholesky


mn=max(m,n)

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

! calculates AD^2A'
! AD^2
xds = x/s ! xs^-1
DO j=1,n
  AD2(:,j) = A(:,j)*xds(j)
  A3(:,j)  = A(:,j)/s(j)
ENDDO
!print*, "AD2=[[", ((dwork(pdA2+(j-1)*m+i-1),j=1,n),char(10), i=1,m), "]]"
AD2A = matmul(AD2, transpose(A))

#ifdef DEBUG
print*, "AD2A=[[", ((AD2A(i,j),j=1,m),char(10), i=1,m), "]]"
#endif

! calculates deltalambda(val provisoire) = -rb-AXS^-1rc
t2 = -rb-matmul(AD2,rc)

#ifdef DEBUG
print*, "termb=[",  (t2(i),i=1,m), "]'"
#endif

! solves AD^2A' deltalambda = -rb-AXS^-1rc+AS^-1rxs
IF (cho == 0) THEN
  !  lapack cholesky
  CALL dpotrf('L', m, AD2A, m, info) 
ELSE ! cho==1
  ! modified choleski
  CALL modchol2( m, AD2A, delta, AD2A, info )
ENDIF
IF (info .NE. 0) THEN
  RETURN
ENDIF

#ifdef DEBUG
PRINT*, "chol=[[", ((AD2A(i,j), j=1,m), char(10), i=1,m), "]]"
#endif

CALL SOLVETRIGLP ( m, n, A, AD2A, A3, s, xds, rc, rxs, t2, deltax, deltalambda, deltas, info )
IF (info .LT. 0) THEN
  RETURN
ENDIF

#ifdef DEBUG
print*, "deltax=[", (deltax(i),i=1,n), "]'"
print*, "deltalambda=[", (deltalambda(i),i=1,m), "]'"
print*, "deltas=[", (deltas(i),i=1,n), "]'"
#endif

#ifdef DEBUG
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
#endif


RETURN
END


