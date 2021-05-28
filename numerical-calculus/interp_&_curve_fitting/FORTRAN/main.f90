program main
    use fgsl
    use inter_approx
    use system_parameters
    implicit none
    real(wp), allocatable :: nodes(:,:), points(:), coefficients(:,:), results(:,:)
    real(wp), allocatable :: points_spline(:,:), points_approx(:,:), points_interapprox(:,:), aux(:,:)
    real(wp), allocatable :: nodes_approx(:,:)
    real(wp), allocatable :: x(:), y(:), xi(:), yi(:)
    integer  :: limit_min_points, limit_max_points
    integer  :: number_of_nodes, number_of_points, polynomial_order_approx
    integer  :: i,j

    ! Initialize the data
    call read_data

    ! Routines to valuate the points in each line of action and with FGSL
    call interp_approx_routine

    ! Saves the data in a file
    call save_results

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
        allocate(coefficients(4,number_of_nodes))

        read(1,*)
        do i = 1,2
                read(1,*) (nodes(i,j), j = 1,number_of_nodes)
        end do

        read(1,*)
        read(1,*) limit_min_points

        read(1,*)
        read(1,*) limit_max_points

        read(1,*)
        read(1,*) polynomial_order_approx
        ! Closes the input data file
        close(1)
        ! Initialzes the consequential data
        number_of_points = limit_max_points-limit_min_points+1
        allocate(points(number_of_points))
        do i = 1, number_of_points
            points(i) = limit_min_points
            limit_min_points = limit_min_points+1
        end do

    end subroutine read_data

    subroutine interp_approx_routine
        implicit none

        call spline (nodes, coefficients, number_of_nodes)
        allocate(points_spline(2,number_of_points))
        points_spline = ispline(points, nodes, coefficients, number_of_nodes)
        deallocate(coefficients)

        allocate(coefficients(polynomial_order_approx+1,1), points_approx(2,number_of_points))
        call OLS (nodes, coefficients,number_of_nodes,polynomial_order_approx)
        allocate(nodes_approx(2,number_of_nodes))
        points_approx = isOLS(points, nodes, coefficients, number_of_nodes, polynomial_order_approx)
        nodes_approx = isOLS(nodes(1,:), nodes, coefficients, number_of_nodes, polynomial_order_approx)
        deallocate(coefficients)

        allocate(points_interapprox(3,number_of_points),aux(2,number_of_points))
        aux = isLagrange(nodes_approx, points,number_of_nodes)
        points_interapprox(1,:) = aux(1,:)
        points_interapprox(2,:) = aux(2,:)

        aux =  isDDint (nodes_approx, points, number_of_nodes)
        points_interapprox(3,:) = aux(2,:)

        allocate(results(7,number_of_points))
        do i = 1,number_of_points
            results(1,i) = points(i)
            results(2,i) = points_spline(2,i)
            results(3,i) = points_approx(2,i)
            results(4,i) = points_interapprox(2,i)
            results(5,i) = points_interapprox(3,i)
        end do

        deallocate(points_interapprox, points_spline,points_approx)
        allocate(points_interapprox(2,number_of_points), points_spline(2,number_of_points))
        points_spline(1,:) = points(:)
        points_interapprox(1,:) = points(:)

        allocate(x(number_of_nodes), y(number_of_nodes), xi(number_of_points), yi(number_of_points))
        x(:) = nodes(1,:)
        y(:) = nodes(2,:)
        xi(:) = points(:)
        call interpolate_fgsl_cspline(x, y, xi, yi)
        points_spline(2,:) = yi

        y(:) = nodes_approx(2,:)
        xi(:) = points(:)
        call interpolate_fgsl_polynomial(x, y, xi, yi)
        points_interapprox(2,:) = yi

        deallocate(x,y,xi,yi)

        do i = 1,number_of_points
            results(6,i) = points_spline(2,i)
            results(7,i) = points_interapprox(2,i)
        end do

        deallocate(nodes, points)
        deallocate(points_spline, points_interapprox, nodes_approx, aux)

    end subroutine interp_approx_routine

    subroutine save_results
        implicit none
        integer :: i, j

        open(2,file = 'trabalho4_FORTRAN_results.txt', action='write')
        do i=1,number_of_points
            write(2,*) results(:,i)
        end do
        close(2)

    end subroutine save_results

end program main
