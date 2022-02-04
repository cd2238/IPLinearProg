!> Main program in order to execute the linear program (linprosimp)
!> Example of the command line : linpro n10c04
!> the maximal number of iterations itermax, and the duality measure tolerance mutol has to be modified in the code below.
!-----------------------------------------------------------------------
program main

     implicit none
     integer :: num_args, allocock, i
     character(len=12), dimension(:), allocatable   :: args
     integer                                        :: nb_control_variables
     double precision, dimension(:), allocatable    :: linear_objective
     integer                                        :: nb_inequality_constraint
     double precision, dimension(:,:), allocatable  :: constraint_matrix
     double precision, dimension(:), allocatable    :: constraint_vector
     integer                                        :: itermax
     double precision                               :: mutol


     double precision, dimension(:), allocatable    :: x    
     double precision, dimension(:), allocatable    :: lambda   
     double precision, dimension(:), allocatable    :: s
     integer                                        :: iter
     integer                                        :: info

     ! parameters
     parameter(itermax = 1000, mutol=1.0d-10)
     ! read data
     num_args = command_argument_count()
     allocate(args(num_args))  ! 

     if (num_args >= 1) then
         call get_command_argument(1,args(1))
         print*, "dataset=", args(1)
     else
         print*, "1 argument is needed (name of the data folder)."
         return
     end if
     open(12, file="../data/"//trim(args(1))//"/nb_control_variables.txt")
     read(12,*) nb_control_variables
     print*, nb_control_variables
     close(12)
     open(12, file="../data/"//trim(args(1))//"/nb_inequality_constraint.txt")
     read(12,*) nb_inequality_constraint
     print*, nb_inequality_constraint
     close(12) 
     
     ! allocate arrays
     allocate (linear_objective(nb_control_variables), stat = allocock)
     if (allocock /= 0) return
     allocate (constraint_matrix(nb_inequality_constraint, nb_control_variables), stat = allocock)
     if (allocock /= 0) return    
     allocate (constraint_vector(nb_inequality_constraint), stat = allocock)
     if (allocock /= 0) return   
     
     open(12, file="../data/"//trim(args(1))//"/linear_objective.txt")
     read(12,*) linear_objective
     print*, linear_objective
     close(12)
     open(12, file="../data/"//trim(args(1))//"/constraint_matrix.txt")
     read(12,*) constraint_matrix
     close(12) 
     open(12, file="../data/"//trim(args(1))//"/constraint_vector.txt")
     read(12,*) constraint_vector
     print*, constraint_vector
     close(12) 

     ! outputs
     allocate (x(nb_control_variables), stat = allocock)
     if (allocock /= 0) return 
     allocate (lambda(nb_inequality_constraint), stat = allocock)
     if (allocock /= 0) return 
     allocate (s(nb_control_variables), stat = allocock)
     if (allocock /= 0) return    
     
 
     ! perform optimization
     call linprosimp ( nb_control_variables, linear_objective, &
                       nb_inequality_constraint, constraint_matrix, constraint_vector, &
                       itermax, mutol, &
                       x, lambda, s, iter, info )

     
#ifdef DEBUG
     ! write results
     print*, (x(i), i=1, nb_control_variables)
     print*, (lambda(i), i=1, nb_inequality_constraint)     
     print*, (s(i), i=1, nb_control_variables)     
#endif

     print*, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
     print*, "Results :"
     print*, "x =", (x(i), i=1,nb_control_variables)
     print*, "obj =", dot_product(linear_objective,x)
     print*, "constraints violated ? (violation if > 0)"
     do i = 1, nb_inequality_constraint
        print*, "const ", i, ":", constraint_vector(i) - dot_product(constraint_matrix(i,:), x)
     enddo
     print*, "iter/itermax :", iter, "/", itermax
     print*, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"


     deallocate(x)
     deallocate(s)
     deallocate(lambda)
     deallocate(args)
     deallocate(linear_objective)
     deallocate(constraint_matrix)
     deallocate(constraint_vector)

end program
