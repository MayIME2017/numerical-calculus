!Written by Mayara Magalhães Carvalho
!Introduction to scientific computing
!March, 2021
!compilation: gfortran -o machine_precision machine_precision.f90

!This program calculates the machine precision
!It returns the smallest number eps such that 1+eps>1
program epsilon
implicit none

integer, parameter :: DP=8, SP=4
!Setting the possibilities to the number of reference - real
!number single or double precision
!Epsilon must be the same type
real(DP) :: eps_d, x_d
real(SP) :: eps_s, x_s

print *, 'Define the numbers of reference:'
print *, 'Double precision'
read *,x_d
print *, 'Single precision'
read *,x_s

!Initializing epsilon values
eps_d=1.0_DP
eps_s=1.0_SP
do while (1.0_SP+eps_s .gt. 1.0_SP)
  eps_s=eps_s/2
end do
do while (1.0_DP+eps_d .gt. 1.0_DP)
  eps_d=eps_d/2
end do

!Because the loop does one more division after finding
!the naumber that satisfies the condition
!its necessary to multiply the result by 2
print *,"single precision epsilon=",2*eps_s
print *,"double precision epsilon=",2*eps_d

end program epsilon
