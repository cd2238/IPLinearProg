!> Linear programming needing initial point (primal and duals)
!> from Nocedal Wright
!>
!> @param n the number of control variables
!> @param x0 the initial primal solution
!> @param lambda0 the initial dual solution
!> @param s0 the initial dual solution
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

subroutine linprosimpx ( n, x0, lambda0, s0, c, m, A, b, itermax, mutol, x, lambda, s, mu, iter, info )
!-----------------------------------------------------------------------
implicit none
! inputs
integer                         , intent(in)  :: n
double precision, dimension(n),   intent(in)  :: x0
double precision, dimension(m),   intent(in)  :: lambda0
double precision, dimension(n),   intent(in)  :: s0
double precision, dimension(n),   intent(in)  :: c
integer                         , intent(in)  :: m
double precision, dimension(m,n), intent(in)  :: A
double precision, dimension(m),   intent(in)  :: b
integer                         , intent(in)  :: itermax

!outputs
double precision, dimension(n),   intent(out) :: x
double precision, dimension(m),   intent(out) :: lambda
double precision, dimension(n),   intent(out) :: s
double precision                , intent(out) :: mu
integer                         , intent(out) :: iter
integer                         , intent(out) :: info

! local variables
integer i
#ifdef DEBUG
integer j
#endif
double precision alpha_dual, alpha_dual_aff, alpha_pri, alpha_pri_aff
double precision mutol, eta, mini, mu_aff, ps, sigma
double precision, dimension(:),   allocatable :: delta_x
double precision, dimension(:),   allocatable :: delta_lambda
double precision, dimension(:),   allocatable :: delta_s
double precision, dimension(:),   allocatable :: rc
double precision, dimension(:),   allocatable :: rb
double precision, dimension(:),   allocatable :: rxs
double precision, dimension(:),   allocatable :: temp1
double precision, dimension(:),   allocatable :: temp2
double precision, dimension(:),   allocatable :: delta_x_aff
double precision, dimension(:),   allocatable :: delta_s_aff
double precision, dimension(:,:), allocatable :: AD2A
double precision, dimension(:,:), allocatable :: A3
double precision, dimension(:),   allocatable :: xds
double precision, dimension(:),   allocatable :: t2

parameter(eta=0.99d0)
external affinescaling, solvetriglp

! allocations
allocate(delta_x(n))
allocate(delta_lambda(m))
allocate(delta_s(n))
allocate(rc(n))
allocate(rb(m))
allocate(rxs(n))
allocate(temp1(n))
allocate(temp2(n))
allocate(delta_x_aff(n))
allocate(delta_s_aff(n))
allocate(AD2A(m,m))
allocate(A3(m,n))
allocate(xds(n))
allocate(t2(m))


#ifdef DEBUG
  PRINT*, "linprosimpx"
  PRINT*, "x0=", (x0(i), i=1,n)
  PRINT*, "lambda0=", (lambda0(i), i=1,m)
  PRINT*, "s0=", (s0(i), i=1,n)
  PRINT*, "c=", (c(i), i=1,n)
  PRINT*, "A=", ((A(i,j), j=1,n),char(10), i=1,m)
  PRINT*, "b=", (b(i), i=1,m)
  PRINT*, "itermax=", itermax
#endif

! init
info = 0
iter = 0
x = x0
lambda = lambda0
s = s0
ps = dot_product(x, s)
mu = ps/dble(n)

#ifdef DEBUG
  print*, "mu=", mu
  print*, "iter=", iter
  print*, "mutol=", mutol
  print*, "itermax=", itermax
#endif


! main loop
DO WHILE ((mu .GT. mutol).AND.(iter .LT.itermax))
    iter = iter + 1
#ifdef DEBUG
    print*, "iter=", iter
    print*, "mu=", mu
#endif

    ! rc=A'lambda+s-c, rb=Ax-b, rxs=XSe
    rc = matmul(transpose(A), lambda)+ s - c
    rb = matmul(A,x) - b
    rxs= x*s
#ifdef DEBUG
    print*, "rb=", rb
#endif

    ! solve 14.30
    CALL affinescaling( m, n, A, s, x, rc, rb, rxs,  &
                  delta_x, delta_lambda, delta_s, AD2A, A3, &
                  xds, t2, info)
    IF (info .NE. 0) RETURN



#ifdef DEBUG
    print*, "sum(x)=", sum(x)
    print*, "sum(x+delta_x)=", sum(x+delta_x)
    print*, "1:Ax-b=", matmul(A, x+delta_x)-b
#endif


    ! calculate alpha_pri_aff, alpha_dual_aff, mu_aff
    mini=huge(0.0d0)
    DO i=1,n
      IF (delta_x(i) .LT. 0) THEN
        mini = min(-x(i)/delta_x(i), mini)
      ENDIF
    ENDDO
    alpha_pri_aff = min(1.0, mini)

    mini=huge(0.0d0)
    DO i=1,n
      IF (delta_s(i) .LT. 0.0d0) THEN
        mini = min(-s(i)/delta_s(i), mini)
      ENDIF
    ENDDO
    alpha_dual_aff = min(1.0, mini)

    temp1 = x+delta_x*alpha_pri_aff
    temp2 = s+delta_s*alpha_dual_aff
    ps = dot_product(temp1, temp2)
    mu_aff = ps/dble(n)
#ifdef DEBUG
    print*, "mu_aff=", mu_aff
#endif

    ! centering parameter
    sigma = (mu_aff/mu)**3

    ! solve 14.35
    ! rxs = XSe + DeltaXaff DeltaSaff e - sigma mu e
    rxs = rxs + delta_x*delta_s - sigma * mu
#ifdef DEBUG
    print*, (rxs(i), i=1,n)
    print*, "-------"
#endif

    CALL SOLVETRIGLP ( m, n, A, AD2A, A3, s, xds, rc, rxs, t2, delta_x, delta_lambda, delta_s, info )
    IF (info .NE. 0) RETURN
#ifdef DEBUG
    print*, "sum(x)=", sum(x)
    print*, "2:Ax-b=", matmul(A, x+delta_x)-b
    print*, "-------"
#endif

    ! alpha_k_pri, alpha_k_dual
    mini=huge(0.0d0)
    DO i=1,n
      IF (delta_x(i) .LT. 0) THEN
        mini = min(-x(i)/delta_x(i), mini)
      ENDIF
    ENDDO
    alpha_pri = min(1.0, eta*mini)


    mini=huge(0.0d0)
    DO i=1,n
      IF (delta_s(i) .LT. 0) THEN
        mini = min(-s(i)/delta_s(i), mini)
      ENDIF
    ENDDO
    alpha_dual = min(1.0, eta*mini)
#ifdef DEBUG
    print*, "s=", s
    print*, "alpha_dual=", alpha_dual
    do i=1,n
      if (s(i) < 0) then
        print*, "s0 négatif!", s
      endif
    enddo
#endif

    ! update

    x = x + alpha_pri*delta_x
    lambda = lambda + alpha_dual*delta_lambda
    s  = s + alpha_dual*delta_s


    ! calculate mu
    ps = dot_product(x, s)
    mu = ps/dble(n)
#ifdef DEBUG
    print*, "mu=", mu
    print*, "iter=", iter
#endif

ENDDO





RETURN
END


