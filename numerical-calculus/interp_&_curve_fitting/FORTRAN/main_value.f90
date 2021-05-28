program main_value
    use inter_approx
    use system_parameters
    implicit none
    real(wp), allocatable :: nodes(:,:), coefficients(:,:), points(:)
    real(wp) :: extrapolation_point
    integer  :: number_of_nodes, polynomial_order_approx
    integer  :: i,j

    ! Initialize the data
    call read_data

    allocate(coefficients(polynomial_order_approx+1,1), points(1))
    points(:) = extrapolation_point
    call OLS (nodes, coefficients,number_of_nodes,polynomial_order_approx)
    print *, isOLS(points, nodes, coefficients, number_of_nodes, polynomial_order_approx)

    deallocate(coefficients)


contains
    subroutine read_data
        implicit none
        integer :: j
        !Open file with the data for the interpolation and approximation
        open(1,file = 'Inter_Approx_Input.txt', status = 'old')
        ! Reads each piece of data that will be used
        read(1,*)
        read(1,*) number_of_nodes

        allocate(nodes(2,number_of_nodes))

        read(1,*)
        do i = 1,2
                read(1,*) (nodes(i,j), j = 1,number_of_nodes)
        end do

        read(1,*)
        read(1,*)

        read(1,*)
        read(1,*)

        read(1,*)
        read(1,*) polynomial_order_approx

        read(1,*)
        read(1,*) extrapolation_point
        ! Closes the input data file
        close(1)

    end subroutine read_data

end program main_value
