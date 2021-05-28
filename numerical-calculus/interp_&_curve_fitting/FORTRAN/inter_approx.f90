!----------------------------------------------------------------!
!Module with the interpolation and approximation                 !
!Written by Mayara Magalhães Carvalho                            !
!Brazil - Military Institute of Technology                       !
!May,2021                                                      !
!Last edited - 05/21/2021                                        !
!----------------------------------------------------------------!

module inter_approx
    use system_parameters
implicit none

contains
    subroutine spline (arg_nodes, coefficients, number_of_nodes)
 !======================================================================
 !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
 !  for cubic spline interpolation
 !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
 !  for  x(i) <= x <= x(i+1)
 !----------------------------------------------------------------------
 !  input..
 !  arg_nodes       = matrix of data abscissas and ordinates
 !  number_of_nodes = number of nodes for the interpolation
 !  output..
 !  coefficients    = matrix of spline coefficients
 !======================================================================
 implicit none
 integer number_of_nodes
 real(wp):: arg_nodes(:,:), b(number_of_nodes), c(number_of_nodes), d(number_of_nodes), coefficients(:,:)
 integer :: i, j, gap
 real(wp):: h

 gap = number_of_nodes-1
 ! check input
 if ( number_of_nodes < 2 ) return
 if ( number_of_nodes < 3 ) then
   b(1) = (arg_nodes(2,2)-arg_nodes(2,1))/(arg_nodes(1,2)-arg_nodes(1,1))   ! linear interpolation
   c(1) = 0.
   d(1) = 0.
   b(2) = b(1)
   c(2) = 0.
   d(2) = 0.
   return
 end if
 !
 ! step 1: preparation
 !
 d(1) = arg_nodes(1,2) - arg_nodes(1,1)
 c(2) = (arg_nodes(2,2) - arg_nodes(2,1))/d(1)
 do i = 2, gap
   d(i) = arg_nodes(1,i+1) - arg_nodes(1,i)
   b(i) = 2.0*(d(i-1) + d(i))
   c(i+1) = (arg_nodes(2,i+1) - arg_nodes(2,i))/d(i)
   c(i) = c(i+1) - c(i)
 end do
 !
 ! step 2: end conditions
 !
 b(1) = -d(1)
 b(number_of_nodes) = -d(number_of_nodes-1)
 c(1) = 0.0
 c(number_of_nodes) = 0.0
 if(number_of_nodes /= 3) then
   c(1) = c(3)/(arg_nodes(1,4)-arg_nodes(1,2)) - c(2)/(arg_nodes(1,3)-arg_nodes(1,1))
   c(number_of_nodes) = c(number_of_nodes-1)/(arg_nodes(1,number_of_nodes) &
   -arg_nodes(1,number_of_nodes-2)) &
   - c(number_of_nodes-2)/(arg_nodes(1,number_of_nodes-1)-arg_nodes(1,number_of_nodes-3))
   c(1) = c(1)*d(1)**2/(arg_nodes(1,4)-arg_nodes(1,1))
   c(number_of_nodes) = -c(number_of_nodes)*d(number_of_nodes-1)**2/(arg_nodes(1,number_of_nodes)-arg_nodes(1,number_of_nodes-3))
 end if
 !
 ! step 3: forward elimination
 !
 do i = 2, number_of_nodes
   h = d(i-1)/b(i-1)
   b(i) = b(i) - h*d(i-1)
   c(i) = c(i) - h*c(i-1)
 end do
 !
 ! step 4: back substitution
 !
 c(number_of_nodes) = c(number_of_nodes)/b(number_of_nodes)
 do j = 1, gap
   i = number_of_nodes-j
   c(i) = (c(i) - d(i)*c(i+1))/b(i)
 end do
 !
 ! step 5: compute spline coefficients
 !
 b(number_of_nodes) = (arg_nodes(2,number_of_nodes) - arg_nodes(2,gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(number_of_nodes))
 do i = 1, gap
   b(i) = (arg_nodes(2,i+1) - arg_nodes(2,i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
   d(i) = (c(i+1) - c(i))/d(i)
   c(i) = 3.*c(i)
 end do
 c(number_of_nodes) = 3.0*c(number_of_nodes)
 d(number_of_nodes) = d(number_of_nodes-1)

 coefficients(1,:) = arg_nodes(2,:)
 coefficients(2,:) = b
 coefficients(3,:) = c
 coefficients(4,:) = d
 end subroutine spline

 function ispline(arg_points, arg_nodes, coefficients, number_of_nodes)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point u
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! points          = the abscissa at which the spline is to be evaluated
! arg_nodes       = matrix of given data points
! coeficients     = matrix of spline coefficients computed by spline
! number_of_nodes = the number of data points
! output:
! ispline = interpolated value at required points
!=======================================================================
implicit none
integer  :: number_of_nodes
real(wp) ::  arg_points(:), arg_nodes(:,:), coefficients(:,:)
integer  :: i, j, k, q
real(wp) :: dx
integer  :: number_of_points
real(wp), allocatable :: ispline(:,:)

number_of_points = size(arg_points)
allocate(ispline(2,number_of_points))

ispline(1,:) = arg_points(:)

! if points is ouside the abscissas' nodes interval take a boundary value (left or right)
do i = 1,number_of_points
    if(arg_points(i) <= arg_nodes(1,1)) then
        ispline(2,i) = arg_nodes(2,1)
        return
    end if
    if(arg_points(i) >= arg_nodes(1,number_of_nodes)) then
        ispline(2,i) = arg_nodes(2,number_of_nodes)
        return
    end if
end do
!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
do q = 1,number_of_points
    i = 1
    j = number_of_nodes+1
    do while (j > i+1)
        k = (i+j)/2
        if(arg_points(q) < arg_nodes(1,k)) then
            j=k
        else
            i=k
        end if
    end do
    !*
    !  evaluate spline interpolation
    !*
    dx = arg_points(q) - arg_nodes(1,i)
    ispline(2,q) = coefficients(1,i) + dx*(coefficients(2,i) + dx*(coefficients(3,i) + dx*coefficients(4,i)))
end do

end function ispline

subroutine OLS (arg_nodes, coefficients,number_of_nodes,polynomial_order_approx)
    implicit none
    integer number_of_nodes
    real(wp)  :: arg_nodes(:,:), coefficients(:,:), matrix_xy(number_of_nodes,1)
    real(wp)  :: matrix_powers_x(number_of_nodes,polynomial_order_approx+1)
    real (wp) :: aa(polynomial_order_approx+1,polynomial_order_approx+1), bb(polynomial_order_approx+1,1)
    integer   :: i, j, k, info, polynomial_order_approx
    integer, allocatable :: pivot(:)

    allocate(pivot(polynomial_order_approx+1))

    ! Building the matrix of powers of x
    matrix_powers_x = 0.d0
    do j = 1,polynomial_order_approx+1
        do i = 1,number_of_nodes
            matrix_powers_x(i,j) = arg_nodes(1,i)**(polynomial_order_approx+1-j)
        end do
    end do
    aa = matmul(transpose(matrix_powers_x),matrix_powers_x)

    ! Building the matrix of products xy
    do i = 1,number_of_nodes
        matrix_xy(i,1) = arg_nodes(2,i)
    end do
    bb = matmul(transpose(matrix_powers_x),matrix_xy)

    call PLU(aa,polynomial_order_approx+1,pivot,info)

    if (info.EQ.1) then
        write(*,*) ' The system matrix is singular, no solution !'
        STOP
    else
        call LUsolver(aa,polynomial_order_approx+1,pivot,bb)
    end if

    coefficients = bb

    deallocate(pivot)

end subroutine OLS

function isOLS (arg_points, arg_nodes, coefficients, number_of_nodes, polynomial_order_approx)
    implicit none
    integer  :: number_of_nodes, polynomial_order_approx
    real(wp) ::  arg_points(:), arg_nodes(:,:), coefficients(:,:)
    integer  :: i,j
    integer  :: number_of_points
    real(wp), allocatable :: isOLS(:,:)

    number_of_points = size(arg_points)
    allocate(isOLS(2,number_of_points))

    isOLS(1,:) = arg_points(:)
    isOLS(2,:) = 0.d0

    do i = 1,number_of_points
        do j = 1,polynomial_order_approx+1
            isOLS(2,i) = isOLS(2,i) + coefficients(j,1)*arg_points(i)**(polynomial_order_approx+1-j)
        end do
    end do

end function isOLS


subroutine PLU(a,n,pivot,info)
    !In order not to conflict with matrix definition in the main program
    !a - Matrix that defines the linear system
    !b - Matrix that defines the linear system result vector
    !n - number of rows and columns of the system's matrix
    !Defining the auxiliar matrixes and the LU matrix of the method
    real(wp)              :: a(:,:)
    real(wp)              :: a_element_max, aux_sum, aux_sum2, coeff
    integer, allocatable  :: pivot(:), aux_pivot(:)
    real(wp), allocatable :: aCopy(:,:),L(:,:),U(:,:)
    integer               :: i,j,k,n,imax,info,kmax
    integer               :: err = epsilon(1.d0)

    allocate(aux_pivot(n),aCopy(n,n),L(n,n),U(n,n))

    info = 0

    aCopy = a

    !Discovering if linear system matrix is singular
    if(DET(aCopy).LT.err) then
        info = 1 !matrix a is singular
        RETURN
    end if

    aCopy = a

    !Finding the pivot indices for the open blas function
    aux_pivot = [ ( i, i=1,n ) ]
    do k = 1,n-1
        kmax = maxloc(abs(aCopy(aux_pivot(k:),k)),1) + k-1
        if (kmax /= k ) then
            aux_pivot([k, kmax]) = aux_pivot([kmax, k])
            aCopy([k, kmax],:) = aCopy([kmax, k],:)
        end if
        aCopy(k+1:,k) = aCopy(k+1:,k) / aCopy(k,k)
        forall (j=k+1:n) aCopy(k+1:,j) = aCopy(k+1:,j) - aCopy(k,j)*aCopy(k+1:,k)
    end do

    do i=1,n
        do j=1,n
            if (aux_pivot(j) == i) then
                pivot(i) = j
            end if
        end do
    end do

    aCopy = a
    !Pivoting the original matrix
    do i = 1,n
        a(i,:)=aCopy(pivot(i),:)
    end do

    !Finding the LU decomposition for the pivoted matrix
    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0

    ! step 1: forward elimination
    do k=1, n-1
        do i=k+1,n
            coeff=a(i,k)/a(k,k)
            L(i,k) = coeff
            do j=k+1,n
                a(i,j) = a(i,j)-coeff*a(k,j)
            end do
        end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,n
        L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,n
        do i=1,j
            U(i,j) = a(i,j)
        end do
    end do

    L = transpose(L)
    U = transpose(U)

    do i=1,n
        L(i,i) = 0.0
    end do

    a = L+U

    deallocate(aux_pivot,L,U,aCopy)

end subroutine PLU

subroutine LUsolver(a,n,pivot,b)
   real(wp) :: a(:,:), b(:,:)
   integer  :: pivot(:)
   integer  :: i,j,n, aux
   real(wp) :: aux_sum
   real(wp), allocatable :: L(:,:),U(:,:)

   allocate(L(n,n), U(n,n))

   !Partially pivoting b
   do i=1,n
       aux = pivot(i)
       aux_sum = b(aux,1)
       b(aux,1) = b(i,1)
   end do

   do j = 1,n
       do i = 1,n
           if (i .gt. j) then
               L(i,j) = 0.d0
               U(i,j) = a(i,j)
           else if (i .eq. j) then
               L(i,j) = 1.d0
               U(i,j) = 0.d0
           else
               L(i,j) = a(i,j)
               U(i,j) = 0.d0
           end if
       end do
   end do

   !Solving Ly=b, where y=Ux
   do i = 1,n
       b(i,1) = 2*b(i,1)-dot_product(L(:,i),b(:,1))
   end do

   !Solving Ux=y
   do i = n,1,-1
       b(i,1) = (b(i,1)-dot_product(U(:,i),b(:,1)))/a(i,i)
   end do

    deallocate(L,U)
end subroutine LUsolver

real(wp) function DET(aa)
    real(wp) :: aa(:,:)
    real(wp) :: tmp,c(size(aa,dim=1),size(aa,dim=2))
    real(wp) :: max
    integer  :: i,j,k,l,m,num(size(aa,dim=1))
    integer :: n

    n=size(aa,dim=1)
	   det=1.
	      do k=1,n
              max=aa(k,k);num(k)=k;
              do i=k+1,n
                  if(abs(max)<abs(aa(i,k))) then
                      max=aa(i,k)
                      num(k)=i
                  end if
              end do
              if (num(k)/=k) then
                  do l=k,n
                      tmp=aa(k,l)
                      aa(k,l)=aa(num(k),l)
				      aa(num(k),l)=tmp
                  end do
                  det=-1.*det
              end if
              do m=k+1,n
                  c(m,k)=aa(m,k)/aa(k,k)
                  do l=k,n
                      aa(m,l)=aa(m,l)-c(m,k)*aa(k,l)
                  end do
              end do
          end do
          do i=1,n
              det=det*aa(i,i)
          end do
          return
      end function

function isLagrange (arg_nodes, arg_points, number_of_nodes)
!=========================================================================
! Evaluates the Lagrange interpolant
!-------------------------------------------------------------------------
! inputs:
!  number_of_nodes = number of data points
!  arg_nodes       = data points (x,y)
!  points          = evaluations points
! outputs:
!  isLagrange      = interpolated values
!=========================================================================
    implicit none
    integer                :: number_of_nodes, number_of_points
    real (wp)              :: arg_nodes(:,:), arg_points(:)
    real (wp), allocatable :: isLagrange(:,:)
    integer                :: i, j, k
    real (wp), allocatable :: lb(:,:), aux_matrix(:,:)
    real (wp)              :: aux_product

    number_of_points = size(arg_points)

    allocate(isLagrange(2,number_of_points),lb(number_of_points,number_of_nodes), aux_matrix(number_of_nodes,2))

    isLagrange(1,:) = arg_points(:)

    do i = 1,number_of_points
        do j = 1,number_of_nodes
            aux_product = 1.0
            do k = 1,number_of_nodes
                if (k .ne. j) then
                    aux_product = aux_product*(arg_points(i) - arg_nodes(1,k))/(arg_nodes(1,j) - arg_nodes(1,k))
                else
                    aux_product = aux_product*1.0
                end if
            end do
            lb(i,j) = aux_product
        end do
    end do

    aux_matrix = transpose(arg_nodes)

    isLagrange(2,:)=matmul(lb,aux_matrix(:,2))

    deallocate(lb,aux_matrix)

  return

end function isLagrange

function isDDint (arg_nodes, arg_points, number_of_nodes)
!=========================================================================
! Evaluates the Newton interpolant
!-------------------------------------------------------------------------
! inputs:
!  number_of_nodes = number of data points
!  arg_nodes       = data points (x,y)
!  arg_points      = evaluations points
! outputs:
!  isLagrange      = interpolated values
!=========================================================================
    implicit none
    integer                :: number_of_nodes, number_of_points
    real (wp)              :: arg_nodes(:,:), arg_points(:)
    real (wp), allocatable :: isDDint(:,:)
    integer                :: i, j, k
    real (wp)              :: aux_sum, aux_product
    real (wp), allocatable :: cd(:,:)

    number_of_points = size(arg_points)

    allocate(isDDint(2,number_of_points),cd(number_of_nodes, number_of_nodes))

    isDDint(1,:) = arg_points(:)

    !
    !  Copy the data values
    !
    cd(:,1) = arg_nodes(2,:)

    !
    !  Compute the divided differences.
    !
      do j = 2, number_of_nodes
        do i = j, number_of_nodes
            cd(i,j) = (cd(i,j-1) - cd(i-1,j-1))/(arg_nodes(1,i)-arg_nodes(1,i+1-j))
        end do
      end do

      do i=1,number_of_points
          aux_sum = arg_nodes(2,1)
          do j = 2,number_of_nodes
              aux_product = 1.0
              do k = 1,j-1
                  aux_product = aux_product*(arg_points(i)-arg_nodes(1,k))
              end do
              aux_sum = aux_sum + aux_product*cd(j,j)
          end do
          isDDint(2,i)=aux_sum
      end do

    deallocate(cd)

  return

end function isDDint

subroutine interpolate_fgsl_cspline(x, y, xi, yi)
    ! Interpolates using FGSL
    ! inputs:
    !  x, y       = data points (x,y)
    !  xi         = evaluations points
    ! outputs:
    !  yi         = interpolated values
    use FGSL
    implicit none
    real(fgsl_double), dimension(:), intent(in) :: x, y
    real(fgsl_double), dimension(:), intent(in) :: xi
    real(fgsl_double), dimension(:), intent(inout) :: yi
    integer(fgsl_size_t) :: n
    integer(fgsl_long) :: i
    integer(fgsl_int) :: status
    type(fgsl_interp_accel) :: acc
    type(fgsl_spline) :: spline
    type(fgsl_interp_type) :: interp_type

    interp_type = fgsl_interp_cspline

    ! Prepares the interpolation/spline routine
    n = size(x)
    ! ALLOCATES AN INTERPOLATION ACCELERATOR OBJECT (OPTIONAL)
    acc = fgsl_interp_accel_alloc()
    ! ALLOCATES MEMORY AND OBJECTS TO THE SPLINE
    spline =  fgsl_spline_alloc(interp_type, n)
    ! CALCULATES THE SPLINE FUNCTIONS
    status = fgsl_spline_init(spline, x, y)
    ! Stores the interpolated values
    do i=1, size(xi)
      yi(i) = fgsl_spline_eval(spline, xi(i), acc)
    end do
    ! Frees FGSL stuff
    !call fgsl_spline_free(spline)
    call fgsl_interp_accel_free(acc)

end subroutine interpolate_fgsl_cspline

subroutine interpolate_fgsl_polynomial(x, y, xi, yi)
    ! Interpolates using FGSL
    ! inputs:
    !  x, y       = data points (x,y)
    ! outputs:
    !  yi         = interpolated values
    use FGSL
    implicit none
    real(fgsl_double), dimension(:), intent(in) :: x, y
    real(fgsl_double), dimension(:), intent(in) :: xi
    real(fgsl_double), dimension(:), intent(inout) :: yi
    integer(fgsl_size_t) :: n
    integer(fgsl_long) :: i
    integer(fgsl_int) :: status
    type(fgsl_interp_accel) :: acc_polynomial
    type(fgsl_spline) :: polynomial
    type(fgsl_interp_type) :: interp_type

    interp_type = fgsl_interp_polynomial


    ! Prepares the interpolation/spline routine
    n = size(x)
    ! ALLOCATES AN INTERPOLATION ACCELERATOR OBJECT (OPTIONAL)
    acc_polynomial = fgsl_interp_accel_alloc()
    ! ALLOCATES MEMORY AND OBJECTS TO THE SPLINE
    polynomial =  fgsl_spline_alloc(interp_type, n)
    ! CALCULATES THE SPLINE FUNCTIONS
    status = fgsl_spline_init(polynomial, x, y)
    ! Stores the interpolated values
    do i=1, size(xi)
      yi(i) = fgsl_spline_eval(polynomial, xi(i), acc_polynomial)
    end do
    ! Frees FGSL stuff
    call fgsl_interp_accel_free(acc_polynomial)

end subroutine interpolate_fgsl_polynomial

end module
