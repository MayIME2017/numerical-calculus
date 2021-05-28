!Calculates the roots in each given method

module RootsMethods
    implicit none

contains
    subroutine bisection(func,initial_interval_value,final_interval_value,eps, &
        & max_iterations,key,root,flag)
        integer, parameter :: DP=8, SP=4
        real(DP) :: func, initial_interval_value, final_interval_value, eps, root, dx
        integer:: flag, max_iterations, key
        integer :: counter
    dx=final_interval_value-initial_interval_value

    !Checking Bissection condition method
    if(func(initial_interval_value)*func(final_interval_value) .gt. 0.0) then
      flag = 0
      return
    end if

    !Method iteration
    do while(abs(dx) .gt. eps) !Checks the tolerance for the interval
        !Limiting the number of iterations
        if (counter .gt. max_iterations) then
            flag = 0
            STOP
        end if
    if (key==1) then
        root=(initial_interval_value+final_interval_value)/2.0
    else
        root=final_interval_value-func(final_interval_value)*dx/(func(final_interval_value)-func(initial_interval_value))
    end if

    if (func(initial_interval_value)*func(root) .lt. 0) then
        final_interval_value=root
    else
        initial_interval_value=root
    end if

    dx=final_interval_value-initial_interval_value
    counter=counter+1

    end do

    flag = counter
    end subroutine bisection

    subroutine newtonraphson(func,dfunc,initial_interval_value,final_interval_value,eps, &
        & max_iterations,root,flag)
        integer, parameter :: DP=8, SP=4
        real(DP) :: func, dfunc, initial_interval_value, final_interval_value, eps, root, x1, x2
        integer:: flag, max_iterations
        integer :: counter

        !Defining the firts guess
        x1 = (initial_interval_value+final_interval_value)/2.0

        ! Iterative refining the solution
        do counter=1,max_iterations
            x2 = x1 - func(x1)/dfunc(x1)
            ! condition(s) to stop iterations)
            if(abs(func(x2,0))<= eps) exit
            x1 = x2;
        end do
        Root=x2

        ! check the convergence
        if (counter .le. max_iterations) then
            flag=counter
        else
            flag = 0
        end if
    end subroutine newtonraphson

    subroutine secant(func,initial_interval_value,final_interval_value,eps, &
        & max_iterations,root,flag)
        integer, parameter :: DP=8, SP=4
        real(DP) :: func, initial_interval_value, final_interval_value, eps, root
        real (DP) :: x1, x2, x3, h, df
        integer:: flag, max_iterations
        integer :: counter

        !Initializing the calculations
        x1=(initial_interval_value+final_interval_value)/2.0
        h = 1
        x2 = x1 + h

        ! Iterative refining the solution
        do counter=1,max_iterations
        df = (func(x2)-func(x1))/(x2-x1)
        x3 = x2 - func(x2)/df
        ! condition(s) to stop iterations)
        if(abs(func(x3))<= eps) exit
        x1 = x2;
        x2 = x3;
        end do
        Root=x3

        ! check the convergence
        if (counter .le. max_iterations) then
            flag=counter
        else
            flag = 0
        end if

    end subroutine secant

end module RootsMethods
