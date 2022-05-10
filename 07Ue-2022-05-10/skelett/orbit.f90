PROGRAM ORBIT
   USE constants
   IMPLICIT NONE

!  Define position and velocity Vector/Arrays
   REAL(KIND=p), DIMENSION(2) :: r,v

!  Time and time step
   REAL(KIND=p) :: t, dt
!   INTEGER :: N, MAXN

!  Intermediate quantities for Runge-Kutta (RK4)
!   REAL(KIND=p), DIMENSION(2) :: 

!  Define simulation length and time step size first

   WRITE(*,*)'#Welcome to Orbit calculator!'
!   WRITE(*,*)'#Time step is set to:       ',dt,' days.'
!   WRITE(*,*)'#Number of steps is set to: ',MAXN,'.'

!  Define initial paramters
   t=0
   CALL INIT_EARTH(r,v)
   
   OPEN(unit=11,file="EarthOrbit_Euler.dat")
   OPEN(unit=12,file="EarthOrbit_Energy_Euler.dat")
   WRITE(*,*)' '
   WRITE(*,*)'#Calculating EarthOrbit using first order Euler method.'
   WRITE(*,*)'#Results will be written to EarthOrbit_Euler.dat and EarthOrbit_Energy_Euler.dat files.'
!   DO N=1,MAXN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  IMPLEMENT EXPLICIT EULER INTEGRATION HERE  !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      WRITE(11,*)t, r
      WRITE(12,*)t, ENERGY(Me,r,v)

!   ENDDO

   CLOSE(11)
   CLOSE(12)


!  Define initial paramters
   t=0
   CALL INIT_EARTH(r,v)
   
   OPEN(unit=11,file="EarthOrbit_SemiEuler.dat")
   OPEN(unit=12,file="EarthOrbit_Energy_SemiEuler.dat")
   WRITE(*,*)' '
   WRITE(*,*)'#Calculating EarthOrbit using first order SemiEuler method.'
   WRITE(*,*)'#Results will be written to EarthOrbit_SemiEuler.dat and EarthOrbit_Energy_SemiEuler.dat files.'
!   DO N=1,MAXN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  IMPLEMENT SEMI-EXPLICIT EULER INTEGRATION HERE  !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      WRITE(11,*)t, r
      WRITE(12,*)t, ENERGY(Me,r,v)

!   ENDDO

   CLOSE(11)
   CLOSE(12)



!  Define initial paramters for earth again
   t=0
   CALL INIT_EARTH(r,v)

   OPEN(unit=11,file="EarthOrbit_RK4.dat")
   OPEN(unit=12,file="EarthOrbit_Energy_RK4.dat")
   WRITE(*,*)' '
   WRITE(*,*)'#Calculating EarthOrbit using RK4 method.'
   WRITE(*,*)'#Results will be written to EarthOrbit_RK4.dat and EarthOrbit_Energy_RK4.dat files.'

!   DO N=1,MAXN
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  IMPLEMENT RK4   INTEGRATION HERE  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!     Integrate r and v using RK4

      WRITE(11,*)t, r
      WRITE(12,*)t, ENERGY(Me,r,v)
!   ENDDO
   
   CLOSE(11)
   CLOSE(12)




CONTAINS

!                                                                                           ->
! This function returns the gravitational force of the sun acting on a mass m at coordinate  r  
!

FUNCTION FORCE_SUN(m,r)
   USE constants
   IMPLICIT NONE
   REAL(KIND=p), DIMENSION(2) :: FORCE_SUN
   REAL(KIND=p), DIMENSION(2) :: r      ! position of mass
   REAL(KIND=p) :: m                    ! and mass m

   RETURN
END


!                                                                                  ->
! This function returns the kinetic and potential energy of a mass m at coordinate  r  in the grav. field of the sun
!

FUNCTION ENERGY(m,r,v)
   USE constants
   IMPLICIT NONE
   REAL(KIND=p) :: ENERGY
   REAL(KIND=p) :: m
   REAL(KIND=p), DIMENSION(2) :: r,v

!
!   GESAMTENERGIE MUSS HIER RICHTIG BERECHNET WERDEN.
!

   RETURN
END

!                                 ->    ->
! This routine defines the initial r and v values for the earth.
!
SUBROUTINE INIT_EARTH(r,v)
   IMPLICIT NONE
   REAL(KIND=p), DIMENSION(2) :: r,v

!
!   ANFANGSWERTE MUESSEN HIER EINGESETZT WERDEN
!
!   r = (/ 0.0 , 0.0 /)
!   v = (/ 0.0 , 0.0 /)

END


END PROGRAM



