!> resolves Cholesky robust factorization
!> see Gill, Wright, Murray
!>
!> @param n the dimension of the square matrix A
!> @param A the matrix to triangularize (n,n)
!> @param L the cholesky left triangle matrix (n,n)
!> @param info information code
!-----------------------------------------------------------------------
subroutine modchol ( n, A, delta, allowpermut, L, d, e, info )
!-----------------------------------------------------------------------
!
!     INPUTS
!            n      : size of matrix X                           integer
!            A      : symmetric matrix (n*n)                      double
!            delta  : 
!
!
!     OUTPUTS
!            L      : lower triangular matrix (n*n)               double
!            d      : diagonal of D (n)                           double
!            e      : diagonal modification     (n)               double
!            info   : diagnostic argument                        integer
!
!
!-----------------------------------------------------------------------
!
implicit none
! inputs
integer                         , intent(in)  :: n
double precision, dimension(n,n), intent(in)  :: A
double precision                , intent(in)  :: delta
integer                         , intent(in)  :: allowpermut

! outputs
double precision, dimension(n,n), intent(out) :: L
double precision, dimension(n)  , intent(out) :: d
double precision, dimension(n)  , intent(out) :: e
integer                         , intent(out) :: info

! local variables
integer i,j,s, q
double precision epsilonM, nu, gam, ksi, beta2
double precision dzero, done, summ

double precision, dimension(:,:), allocatable :: C
double precision, dimension(:,:), allocatable :: AA
double precision, dimension(:),   allocatable :: theta
double precision, dimension(:),   allocatable :: prov

! parameters
parameter(dzero= 0.0d0, done = 1.0d0, epsilonM = 1d-15)
!-----------------------------------------------------------------------
!

! Allocation
allocate(C(n,n))
allocate(theta(n))
allocate(AA(n,n))
allocate(prov(n))

! Initialization
info = 0
theta = 0.0d0
e = 0.0d0
d = 0.0d0
C = 0.0d0
AA = A
L = 0.0d0
do i=1,n
  L(i,i) = 1.0d0
enddo


!! step mc1
nu=max(1.0d0, sqrt(dble(n**2-1)))
gam = maxval(abs( (/ (AA(i,i), i=1,n) /) ))



!! step mc2
j=1
do i=1,n
  C(i,i) = AA(i,i)
enddo
ksi = maxval(abs( AA-C ))
beta2= max(gam, ksi/nu, epsilonM)
! main loop
do while (j .LE. n)
  if (allowpermut == 1) then
    !! step mc3
    q = maxloc(abs( (/ (C(i,i), i=j,n) /) ), dim=1) + j-1
    IF (q .GT. j) then
      !! permute row and column q & j
      prov = AA(j,:)
      AA(j,:) = AA(q,:)
      AA(q,:) = prov
    
      prov = AA(:,j)
      AA(:,j) = AA(:,q)
      AA(:,q) = prov
    endif
  endif
  !! mc4
  do s=1,j-1
    L(j,s) = C(j,s)/D(s)
  enddo
  do i=j+1,n
    summ = 0
    do s=1,j-1
      summ = summ + L(j,s)*C(i,s)
    enddo
    C(i,j) = AA(i,j)- summ
    if (j<n) then
      theta(j) = maxval(abs(C(j+1:n,j)))
    else
      theta(j) = 0.0d0
    endif
  enddo

  !! mc5
  D(j) = max(delta, abs(C(j,j)), theta(j)**2/beta2)
  E(j)=D(j) - C(j,j)
  if (j == n) then
    exit
  endif

  !! mc6
  do i=j+1,n
    C(i,i) = C(i,i)-C(i,j)**2/D(j)
  enddo
  j=j+1
enddo

#ifdef DEBUG
print*, "cholmod-L=["
print*,((L(i,j), j=1,n),char(10), i=1,n)
print*, "];"
print*, "cholmod-d=["
print*,(d(i), i=1,n)
print*, "];"
print*, "cholmod-e=["
print*,(e(i), i=1,n)
print*, "];"
#endif

return

end subroutine
