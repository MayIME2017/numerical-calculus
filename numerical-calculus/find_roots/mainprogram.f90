!----------------------------------------------------------------------------------!
!Program that finds a single root of a function using numerical methods            !
!Written by Mayara Magalhães Carvalho                                              !
!Brazil - Military Institute of Technology                                         !
!March,2021                                                                        !
!Last edited - 03/29/2021                                                          !
!Compilation command:                                                              !
!gfortran -o rootsNumerical function.f90 dfunc.f90 RootsMethods.f90 mainprogram.f90!
!----------------------------------------------------------------------------------!

program mainProgram
    use RootsMethods
    implicit none
    integer, parameter :: DP=8, SP=4
    real(DP) :: initial_interval_value, final_interval_value, eps, root
    integer :: flag, max_iterations, l_info, key
    real(DP), external :: func, dfunc

    !Open file with the input - interval, tolerance, maximum number of iteractions
    open (1, file='input.txt', status = 'old', action = 'read', iostat=l_info)

    !Testing errors opening the input file
    if (l_info .ne. 0) then
        print*, 'Problem opening input file'
        STOP
    end if

    !Reading the data in the file
    read (1,*) !Used to skip line with data titles
    read (1,*) eps, initial_interval_value, final_interval_value, key, max_iterations

    close(1)

    !Defining the method used based on the key choosen in the input file
        if (key == 1) then
            write(*,*) ' Method - Bisection'
            call bisection(func,initial_interval_value,final_interval_value,eps,max_iterations,key,root,flag)
        else if (key == 2) then
            write(*,*) ' Method - False Position'
            call bisection(func,initial_interval_value,final_interval_value,eps,max_iterations,key,root,flag)
        else if(key == 3) then
            write (*,*) ' Method - Newton Raphson'
            call newtonraphson(func,dfunc,initial_interval_value,final_interval_value,eps,max_iterations,root,flag)
        else
            write (*,*) ' Method - Secant'
            call secant(func,initial_interval_value,final_interval_value,eps,max_iterations,root,flag)
        end if

        !Writes out the output of the subroutine if the method
        !cannot give a root in the given interval
        if(flag == 0) then
            write(*,*) ' no root found for choosen method'
            STOP
        end if

        !Writes out the output of the subroutine if the method is
        !successful
        write(*,100) eps, flag
        if(flag>0) then
            write(*,101) root
            write(*,102) func(root)
            STOP
        end if

        100 format(' tolerance   = ',1pe12.5,/, ' iterations  = ',i3)
        101 format(' root        = ', 1pe12.5)
        102 format(' func(root)  = ',1pe12.5)

end program mainProgram
