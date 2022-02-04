!> resolves Cholesky robust factorization - encapsulate function
!> see Gill, Wright, Murray
!>
!> @param n the dimension of the square matrix A
!> @param A the matrix to triangularize (n,n)
!> @param L the cholesky left triangle matrix (n,n)
!> @param info information code
!-----------------------------------------------------------------------
subroutine modchol2 ( n, A, L, info )
implicit none
! inputs
integer                         , intent(in)  :: n
double precision, dimension(n,n), intent(in)  :: A
! outputs
double precision, dimension(n,n), intent(out) :: L
integer                         , intent(out) :: info

! local variables
integer allowpermut, i,j
double precision delta
double precision, dimension(:),     allocatable :: d
double precision, dimension(:),     allocatable :: e
#ifdef DEBUG
integer info0
double precision, dimension(:,:),   allocatable :: A2
#endif

parameter(allowpermut=0, delta=1.0d-10)
external modchol
#ifdef DEBUG
external dpotrf
#endif


! initialization
info = 0

! allocations
allocate(d(n))
allocate(e(n))
#ifdef DEBUG
allocate(A2(n,n))
#endif

#ifdef DEBUG
print*, "modchol2-A=["
print*,((A(i,j), j=1,n),char(10), i=1,n)
print*, "];"
#endif

#ifdef DEBUG
! compare with classical cholesky
A2 = A
call dpotrf('L', n, A2, n, info0)
print*, "cholmod-L0=["
print*, ((A2(i,j), j=1,n),char(10), i=1,n)
print*, "];"
#endif


!!! cholesky factorization
call modchol ( n, A, delta, allowpermut, L, d, e, info )
!
#ifdef DEBUG
print*, "cholmod-L2=["
print*,((L(i,j), j=1,n),char(10), i=1,n)
print*, "];"
print*, "cholmod-info0-info=", info0, info
#endif

do i=1,n
  L(i:n,i) = L(i:n,i)*sqrt(d(i))
enddo

#ifdef DEBUG
print*, "cholmod-L3=["
print*,((L(i,j), j=1,n),char(10), i=1,n)
print*, "];"
print*, "cholmod-info0-info=", info0, info
#endif



return

end subroutine
