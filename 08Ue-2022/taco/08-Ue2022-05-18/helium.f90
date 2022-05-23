PROGRAM HELIUM
   USE constants
   USE schroedinger
   USE poisson
   IMPLICIT NONE
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   REAL(p), ALLOCATABLE :: R(:)
   REAL(p) :: RMIN, RMAX,DR
   INTEGER :: NR, I, J

   REAL(p), ALLOCATABLE :: U(:)    ! r*Psi(r)
   REAL(p), ALLOCATABLE :: PSI(:)  ! Psi(r)
   REAL(p), ALLOCATABLE :: UPOT(:) ! r*phi(r)
   REAL(p), ALLOCATABLE :: POT(:)  ! phi(r)

   REAL(p) :: E, Z, EMAX, EMIN, EREP


   WRITE(*,*)''
   WRITE(*,*)'Welcome to the SE-Solver for Helium'
   WRITE(*,*)''

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Specify parameters of radial grid
!
   NR=2000
   RMAX=7.0_p
   RMIN=RMAX/(NR+1)
   DR=ABS(RMAX-RMIN)/(NR-1)
!
!  Allocate required memory
!
   ALLOCATE(U(NR),PSI(NR),UPOT(NR),POT(NR),R(NR))
!
!  Setup radial grid
!
   DO I=1,NR
      R(I)=RMIN+(I-1)*(RMAX-RMIN)/(NR-1)
   ENDDO
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Solve Schroedinger equation for  He+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!  parameter to solve for He+ first
   POT=0.0_p
   Z=2

!  search for lowest eigenvalue in the following interval
   EMAX = 0.0_p
   EMIN = -20._p

   WRITE(*,*) 'Solving Schroedinger equation for He+ '
   CALL SOLVE_SE(EMIN,EMAX,Z,POT,R,NR,DR,U)
   WRITE(*,*) 'Lowest EV is between ',EMIN,' and ',EMAX,' Hartree'
   !write(*,*)u
   CALL SOLVE_POISSON(U,R,NR,DR,UPOT,EREP)
   WRITE(*,*) 'Repulsion energy of electron is ',EREP,' Hartree'

   WRITE(*,*) 'Writing results to Hydrogen_orbital.dat and Hydrogen_orbital_potential.dat'
   OPEN(unit=11,file="Hydrogen_orbital.dat")
   DO I=1,NR
      WRITE(11,*) R(I),U(I)/R(I),UPOT(I)/R(I)
   ENDDO
   CLOSE(11)
   OPEN(unit=11,file="Hydrogen_orbital_potential.dat")
   DO I=1,NR
      WRITE(11,*) R(I),UPOT(I)/R(I)
   ENDDO
   CLOSE(11)


   WRITE(*,*) ''
   WRITE(*,*)''
   WRITE(*,*) 'Solving Schroedinger equation for He '
   WRITE(*,*)''
   DO J=1,100

      EMAX=0.0_p
      EMIN=-20.0_p
      CALL SOLVE_SE(EMIN,EMAX,Z,POT,R,NR,DR,U)
      WRITE(*,*) 'Lowest EV is between ',EMIN,' and ',EMAX,' Hartree'

      CALL SOLVE_POISSON(U,R,NR,DR,UPOT,EREP)
      POT=UPOT/R

   ENDDO

   WRITE(*,*) 'Groundstate energy of the electron is ',(EMAX+EMIN)-EREP,' Hartree and the repulsion energy is ', EREP


END PROGRAM
