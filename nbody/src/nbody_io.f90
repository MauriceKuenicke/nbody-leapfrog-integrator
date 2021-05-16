!!+nbody_io.f90
!!
!! module: provides subroutines to read nbody data, supported format
!!         is one line of header with number of particles N follow by N
!!         lines with particles data (mass, pos, vel --> 7 REALs)
!!       
!!-
MODULE nbody_io
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: load_bodies

CONTAINS

  !
  ! read particles from STDIN and load them into this also allocating
  ! the required memory
  !
  SUBROUTINE load_bodies(n_particles, x, v, a, m, a_old)
    use parameters
    IMPLICIT NONE

    integer                                               :: status, i
    integer, intent(out)                                  :: n_particles
    real(8), dimension (:,:), allocatable,  intent(out)   :: x, v, a, a_old
    real(8), dimension (:), allocatable, intent(out)      :: m
    

    ! Allocate initial particle conditions
    READ*, n_particles
    ALLOCATE(m(n_particles), &
             x(n_particles, 3), &
             a_old(n_particles, 3), & 
             v(n_particles, 3), &
             a(n_particles, 3), STAT=status)     
    IF(status/=0) STOP

    DO i=1, n_particles
        READ*, m(i), x(i, 1), x(i, 2), x(i, 3), &
               v(i, 1), v(i, 2), v(i, 3)
    END DO
  END SUBROUTINE load_bodies
END MODULE nbody_io
