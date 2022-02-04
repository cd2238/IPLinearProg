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

subroutine linprosimpx ( n, x0, lambda0, s0, c, m, A, b, itermax, x, lambda, s, mu, iter, info )
!-----------------------------------------------------------------------
IMPLICIT NONE

INTEGER                         , INTENT(in)  :: n
DOUBLE PRECISION, DIMENSION(n),   INTENT(in)  :: x0
DOUBLE PRECISION, DIMENSION(m),   INTENT(in)  :: lambda0
DOUBLE PRECISION, DIMENSION(n),   INTENT(in)  :: s0
DOUBLE PRECISION, DIMENSION(n),   INTENT(in)  :: c
INTEGER                         , INTENT(in)  :: m
DOUBLE PRECISION, DIMENSION(m,n), INTENT(in)  :: A
DOUBLE PRECISION, DIMENSION(m),   INTENT(in)  :: b
INTEGER                         , INTENT(in)  :: itermax
DOUBLE PRECISION, DIMENSION(n),   INTENT(out) :: x
DOUBLE PRECISION, DIMENSION(m),   INTENT(out) :: lambda
DOUBLE PRECISION, DIMENSION(n),   INTENT(out) :: s
DOUBLE PRECISION                , INTENT(out) :: mu
INTEGER                         , INTENT(out) :: iter
INTEGER                         , INTENT(out) :: info


INTEGER i, j
DOUBLE PRECISION alpha_dual, alpha_dual_aff, alpha_pri, alpha_pri_aff
DOUBLE PRECISION epsimu, eta, mini, mu_aff, ps, sigma

DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_x
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_lambda
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_s
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: rc
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: rb
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: rxs
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: temp1
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: temp2
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_x_aff
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: delta_s_aff
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: AD2A
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: A3
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: xds
DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: t2

PARAMETER(epsimu=1.0d-20, eta=0.99d0)
EXTERNAL SOLVESPARSESYSLP, SOLVETRIGLP

! allocation
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
  print*, "epsimu=", epsimu
  print*, "itermax=", itermax
#endif


! main loop
DO WHILE ((mu .GT. epsimu).AND.(iter .LT.itermax))
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
#endif

    ! update
    do i=1,n
      if (s(i) < 0) then
        print*, "s0 négatif!"
      endif
    enddo
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


