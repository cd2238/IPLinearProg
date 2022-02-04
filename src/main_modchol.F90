!> subroutine to test Cholesky factorization
program main_modchol
!-----------------------------------------------------------------------
implicit none



     integer :: num_args, allocock, i,j
     double precision :: delta
     character(len=12), dimension(:), allocatable   :: args
     
     integer                                        :: n
     double precision, dimension(:,:), allocatable  :: G
  
       
     
     double precision, dimension(:,:), allocatable  :: L                   
     integer                                        :: info
     
     ! read data
     num_args = command_argument_count()
     allocate(args(num_args))  ! 

     if (num_args >= 1) then
         call get_command_argument(1,args(1))
         print*, "dataset=", args(1)
     else
         print*, "1 argument is needed (name of the choleski file)."
         return
     end if
     
     open(12, file="../data/"//trim(args(1))//"/n.txt")
     read(12,*) n
     print*, n
     close(12)

     
     allocate (G(n, n), stat = allocock)
     if (allocock /= 0) return    
     
 
     open(12, file="../data/"//trim(args(1))//"/G.txt")
     read(12,*) G
     close(12) 

     
     allocate (L(n,n), stat = allocock)
     if (allocock /= 0) return   
     
 
          
     
     ! perform optimization
     delta = 1.0d-8
     call modchol2 ( n, G, delta, L, info )

     


     ! write results
     print*, ((G(i,j), j=1,n),char(10), i=1,n)
     print*, ((L(i,j), j=1,n),char(10), i=1,n)



     deallocate(L)
     deallocate(G)
     deallocate(args)

end program main_modchol

