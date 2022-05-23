MODULE SCHROEDINGER
   USE constants
   IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Main module that solves the radial Schroedinger equation
!
!
!    d^2 U(x)
!    ___________  =  ................. U(x)
!      d x^2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SOLVE_SE(EMIN,EMAX,Z,POT,R,NR,DR,U)
   USE numerow
   IMPLICIT NONE
   INTEGER, intent(in) :: NR
   REAL(p), intent(in) :: DR
   REAL(p), intent(inout) :: EMIN,EMAX,Z
   REAL(p),  intent(in) :: POT(NR)
   REAL(p),  intent(inout) :: U(NR)
   REAL(p),  intent(in) :: R(NR)
! local
   REAL(p)  :: G(NR)
   REAL(p)  :: S(NR) ! dummy S=0 for numerow
   REAL(p) :: E ! Trial energy
   INTEGER :: NODES,I
   Real(p) :: ZS


   S=0.0_p  
   E=(EMAX+EMIN)/2.0_p
   ZS = 0

!   TODO: CONSTRUCT A WHILE LOOP THAT FINDS THE SOLUTION TO THE SCHROEDINGER EQUATION:

!         CALL THE REQUIRED ROUTINES IN THE CORRECT ORDER
!         TO OBTAIN THE SOLUTION FOR THE SCHROEDINGER EQUATION FOR EACH TRIAL ENERGY
!         (about 6 lines)
!
!         ONCE THE WAVEFUNCTIOM HAS BEEN OBTAINED COMPUTE ITS NUMBER OF NODES AND
!         MODIFY EMIN and EMAX ACCORDINGLY. COMPUTE NEW TRIAL ENERGY
!         (about 4 lines lines) Wenn 0 Knoten --> Emin = E neue untere Grenze. Wenn 1 Knoten --> Emax = E neue obere Grenze

	do while (abs(E-ZS) > tol)
    
    ZS = E

    CALL INIT_U_FOR_SE(U,R,Z,NR)
    CALL CALC_G(E,Z,POT,R,NR,G)
    CALL INTEGRATE_Y_USING_NUMEROW_FROM_INSIDE(U,G,S,NR,DR)
    CALL NORMALIZE_U(U,R,NR,DR)

    NODES = CALC_NODES(U,NR)
    if(NODES .lt. 1) then
        EMIN = E
    else
        EMAX = E
    endif
        E=(EMAX+EMIN)/2.0_p
	enddo

END SUBROUTINE SOLVE_SE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Calculate function g(x) for radial Schroedinger equation in 
!   the way it is needed for the numerow method (depends on Z, E and \phi(r) ) )
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALC_G(E,Z,POT,R,NR,G)
   IMPLICIT NONE
   REAL(p), intent(in) :: E,Z
   REAL(p), intent(in) :: R(NR),POT(NR)
   INTEGER, intent(in) :: NR
   REAL(p), intent(out) :: G(NR)
! local
   INTEGER :: I
 

!   TODO: CALCULATE THE FUNCTION G(X) FOR THE RADIAL SCHROEDINGER EQUATION AS NEEDED BY THE NUMEROW METHOD
	do I=1,NR
		G(I)=2*(E+Z/R(I)-POT(I))
	enddo


END SUBROUTINE CALC_G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Set proper initial values of U at the boundary
!  as needed by the Numerow method to solve the Schroedinger equation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INIT_U_FOR_SE(U,R,Z,NR)
   IMPLICIT NONE
   REAL(p), DIMENSION(1), intent(inout) :: U(NR)
   REAL(p), DIMENSION(1), intent(in) :: R(NR)
   INTEGER  :: NR
   REAL(p)  :: Z
 

!   TODO: SET PROPEER INITIAL VALUES FOR THE WAVEFUNCTION (U)
!         NEEDED BY NUMEROW METHOD (2 lines)
U(1)=R(1)*(1-R(1)*Z) !Dürfte eh 0 sein
U(2)=R(2)*(1-R(2)*Z)

END SUBROUTINE INIT_U_FOR_SE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Normalize U properly
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE NORMALIZE_U(U,R,NR,DR)
   IMPLICIT NONE
   REAL(p), intent(inout) :: U(NR)
   REAL(p), intent(in) :: R(NR)
   INTEGER, intent(in) :: NR
   REAL(p), intent(in) :: DR
! local
   REAL(p) :: NORM
   INTEGER :: I
 
   NORM=0.0_p

!   TODO: Complete this routine to normalize U() (approximately 3 lines)
	do i = 1,NR
		NORM = NORM + DR*U(i)*U(i)
        !NORM = NORM + 4*acos(-1._p)*DR*U(i)*U(i)
	enddo
	U = U/sqrt(NORM)

END SUBROUTINE NORMALIZE_U

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This function computes the number of nodes of the radial wavefunction U
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION CALC_NODES(U,NR)
   IMPLICIT NONE
   INTEGER :: CALC_NODES
   INTEGER, intent(in) :: NR
   REAL(p), intent(in) :: U(NR)
! local
   INTEGER :: I
 
   CALC_NODES=0

!   TODO: Complete this function to compute the number of nodes of the function U() (approximately 3 lines)
do I=1,NR
if(U(I)*U(I+1)<0) CALC_NODES=CALC_NODES+1
enddo

   RETURN
END FUNCTION CALC_NODES


END MODULE
