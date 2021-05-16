!!+leapfrog.f90
!!
!! module: provides subroutines to calculate positional/velocity updates
!!         based on the Leapfrog algortihm
!!       
!!-
MODULE leapfrog
    use parameters
    IMPLICIT NONE
    PRIVATE
  
    PUBLIC :: advance_position, update_velocities, calculate_force_acceleration
  
  CONTAINS
  
  subroutine advance_position(n_particles,x,v,a,time_step)
    implicit NONE

    integer, intent(in)     :: n_particles
    real(DP), intent(out)   :: x(n_particles, 3)
    real(DP), intent(in)    :: v(n_particles, 3), a(n_particles, 3)
    real(DP), intent(in)    :: time_step

    x = x + time_step*v + 0.5_DP * time_step*time_step * a 

   end subroutine advance_position

subroutine update_velocities(n_particles, v,a, a_old, time_step)
    implicit none
    integer, intent(in)     :: n_particles
    real(DP), intent(in)    :: a(n_particles, 3), a_old(n_particles, 3), time_step
    real(DP), intent(out)   :: v(n_particles, 3)

    v = v + 0.5_DP * time_step * (a_old+a)

end subroutine update_velocities

subroutine calculate_force_acceleration(n_particles, m, x, a, U)
    implicit none
    integer, intent(in)   :: n_particles
    real(DP), intent(in)  :: x(n_particles, 3)
    real(DP), intent(in)  :: m(n_particles)
    real(DP), intent(out) :: a(n_particles, 3),  U

    integer :: i,j
    real(DP) :: distance, distance2, fac
    real(DP), dimension(3)  :: distance_vector

    a(:,:) = 0._DP
    U = 0._DP 
    do i=1, n_particles-1
        do j=i+1, n_particles
            distance_vector = x(i,:) - x(j,:)
            distance2 = sum(distance_vector*distance_vector)
            distance = sqrt(distance2)

            fac = distance*distance*distance
            U = U - G * m(i)*m(j)/distance
            
            a(j,:) = a(j,:) + (m(i)/fac)*(distance_vector)
            a(i,:) = a(i,:) - (m(j)/fac)*(distance_vector)
        end do
    end do

end subroutine calculate_force_acceleration

  END MODULE leapfrog