!> resolves 1 linear system
!>
!> @param m the number of linear constraints
!> @param n the number of control variables
!> @param A the matrix constraint
!> @param chol cholesky left factor of AD2A
!> @param A3 A*S^{-1} (intermediary result to avoid redundancy)
!> @param s the current dual solution
!> @param xds x/s

!> @param rc see Nocedal-Wright
!> @param rxs = XSe see Nocedal-Wright
!> @param AD2A A*D^2*A' (intermediary result to avoid redundancy)
!> @param A3 A*S^{-1} (intermediary result to avoid redundancy)
!> @param t2 -rb-AD^2rc (intermediary result to avoid redundancy)
!> @param deltax = xnew-x
!> @param deltalambda = lambdanew-lambda
!> @param deltas = snew-s
!> @param info


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
!-----------------------------------------------------------------------
subroutine solvetriglp ( m, n, A, chol, A3, s, xds, rc, rxs, t2, deltax, deltalambda, deltas, info )
!-----------------------------------------------------------------------
implicit none
! inputs
integer                         , intent(in) :: m
integer                         , intent(in) :: n
double precision, dimension(m,n), intent(in) :: A
double precision, dimension(m,m), intent(in) :: chol
double precision, dimension(m,n), intent(in) :: A3
double precision, dimension(n),   intent(in) :: s
double precision, dimension(n),   intent(in) :: xds
double precision, dimension(n),   intent(in) :: rc
double precision, dimension(n),   intent(in) :: rxs
double precision, dimension(m),   intent(in) :: t2

! outputs
double precision, dimension(n),   intent(out) :: deltax
double precision, dimension(m),   intent(out) :: deltalambda
double precision, dimension(n),   intent(out) :: deltas
integer                         , intent(out) :: info

! local variables
#ifdef DEBUG
integer i
#endif
double precision, dimension(:), allocatable :: righttermlineq

external dpotrs
integer nrhs

! allocations
allocate(righttermlineq(m))


! initialization
deltalambda = t2  + matmul(A3, rxs)
righttermlineq = deltalambda


! solve linear equation for deltalambda
nrhs = 1
call dpotrs('L', m, nrhs, chol, m, deltalambda, m, info) ! DPOTRS
if (info .NE. 0 ) then
#ifdef DEBUG
  print*, "info dpotrs=", info
#endif
  return
endif

! calculates deltas = -rc-A'deltalambda
deltas = matmul(transpose(A), deltalambda)
deltas = -deltas - rc

! calculates deltax = -S^-1rxs - XS^-1 Deltas
deltax = -rxs/s - xds*deltas

#ifdef DEBUG
print*, "righttermlineq=", righttermlineq
print*, "termgauche=", matmul(matmul(chol, transpose(chol)), deltalambda)

print*, "deltax=[", (deltax(i),i=1,n), "]'"
print*, "deltalambda=[", (deltalambda(i),i=1,m), "]'"
print*, "deltas=[", (deltas(i),i=1,n), "]'"
#endif




return
end


