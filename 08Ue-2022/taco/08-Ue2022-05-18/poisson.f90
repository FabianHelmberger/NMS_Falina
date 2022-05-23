MODULE POISSON
   USE constants
   IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Main module that solves the  Poisson equation
!
!    d^2 UPOT(x)
!    ___________  =  - URHO(X)
!      d x^2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SOLVE_POISSON(U,R,NR,DR,UPOT,EREP)
   USE numerow
   IMPLICIT NONE
   INTEGER, intent(in) :: NR
   REAL(p), intent(in) :: DR
   REAL(p), intent(in) :: U(NR),R(NR) 
   REAL(p), intent(out) :: UPOT(NR)
   REAL(p), intent(out) :: EREP
! local
   INTEGER :: I
   REAL(p), DIMENSION(1) :: URHO(NR),G(NR) ! dummy function G=0 for Numerow routine

   G=0.0_p


!   TODO: CALL THE THREE REQUIRED ROUTINES IN THE CORRECT ORDER
!         TO OBTAIN THE SOLUTION FOR THE POISSON EQUATION HERE (3 lines)
	call CALC_URHO(U,R,NR,URHO)
	call INIT_UPOT(UPOT,R,NR)
	call INTEGRATE_Y_USING_NUMEROW_FROM_OUTSIDE(UPOT,G,-URHO,NR,DR)
	
   EREP=0.0_p
   DO I=1,NR

!   TODO: ENTER APPROPRIATE EXPRESSION FOR
!         REPULSION ENERGY HERE (1 line)
		EREP=EREP+URHO(I)*UPOT(I)*DR

   ENDDO

END SUBROUTINE SOLVE_POISSON

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Calculate charge density from U and write it to URHO
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALC_URHO(U,R,NR,URHO)
   IMPLICIT NONE
   INTEGER, intent(in) :: NR
   REAL(p), intent(in) :: U(NR)
   REAL(p), intent(in) :: R(NR)
   REAL(p), intent(out) :: URHO(NR)
! local
   INTEGER :: I

    !write(*,*)u
!   TODO: CALCULATE CHARGE DENSITY HERE (3 lines)
	do I=1,NR
		URHO(I)=U(I)*U(I)/R(i)
        !URHO(I)=16*acos(-1._p)*acos(-1._p)*U(I)*U(I)/R(i)
    !write(*,*)urho(i)
	enddo

END SUBROUTINE CALC_URHO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Set proper initial values of UPOT at tfhe boundary
!  as needed by the Numerow method to solve the Poisson equatio
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INIT_UPOT(UPOT,R,NR)
   IMPLICIT NONE
   INTEGER, intent(in) :: NR
   REAL(p),intent(inout) :: UPOT(NR)
   REAL(p),intent(in) :: R(NR)
 

!   TODO: SET PROPEER INITIAL VALUES FOR THE POTENTIAL
!         NEEDED BY NUMEROW METHOD (2 lines)
	UPOT(NR)=1._p
	UPOT(NR-1)=1._p
	!UPOT(1)=0.0_p

END SUBROUTINE INIT_UPOT


END MODULE
