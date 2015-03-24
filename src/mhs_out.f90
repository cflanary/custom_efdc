SUBROUTINE MHS_OUT

! Output MHS station data in calibration files
!***************************************************************
!***************************************************************

USE GLOBAL
IMPLICIT NONE

INTEGER::I,J,L,K
INTEGER::ITEMP1,ITEMP2,JTEMP1,JTEMP2,SIGN
REAL::NAN
REAL,DIMENSION(LCM)::zeta,vel_maxc,vel_max,tau_max
REAL::time_efdc
LOGICAL,SAVE::FIRSTTIME=.FALSE.	
!REAL Waterflow
!REAL,DIMENSION(8)::Waterflowtot


! Create files if first call
IF(.NOT.FIRSTTIME)THEN

    ! Velocity calibration file with surface elevation and velocity components at all MHS stations
    OPEN (UNIT=112,FILE='vel_cal.dat', STATUS='REPLACE')
    !WRITE(112,*)'Time,H1,U1,V1,H2,U2,V2,H5,U5,V5,H6,U6,V6,H7,U7,V7'

    ! Flow calibration file with flow rates at all MHS stations
    !OPEN (UNIT=113,FILE='flow_cal.dat', STATUS='REPLACE')
    !WRITE(113,*)'Time,Flow_1,Flow_2,Flow_5,Flow_6a,Flow_6'

    ! Tracer calibration file with salinity and dye conc. at all MHS stations
    OPEN (UNIT=115,FILE='tracer_cal.dat', STATUS='REPLACE')
    !WRITE(115,*)'Time,Salt1,Dye1,Salt2,Dye2,Salt5,Dye5,Salt6,Dye6,Salt7,Dye7'

    ! TSS calibration file with TSS concentrations at all MHS stations
    OPEN (UNIT=105, FILE='tss_cal.dat', STATUS='REPLACE')

    ! Boundary calibration file
    OPEN (UNIT=106, FILE='bry_cal.dat', STATUS='REPLACE')

    ! Shear stress calibration file
    OPEN (UNIT=107, FILE='shear_cal.dat', STATUS='REPLACE')

    ! Initialize variable for first time step
    vel_maxc=0.0
    vel_max=0.0
    tau_max=0.0

	FIRSTTIME=.TRUE.
ENDIF

! EFDC time parameter
time_efdc=DT*FLOAT(N)+TCON*TBEGIN 
time_efdc=time_efdc/86400.

!  Water flux at cross sections
!DO LOC=1,5
!    Waterflowtot(LOC)=0.0
!
!    SELECT CASE (LOC)
!    CASE(1)
!    ITEMP1=118;ITEMP2=118;JTEMP1=120;JTEMP2=95
!    SIGN=-1 !! Positive is out
!    CASE(2)
!    ITEMP1=49;ITEMP2=49;JTEMP1=13;JTEMP2=15
!    SIGN=1
!    CASE(3)
!    ITEMP1=37;ITEMP2=42;JTEMP1=202;JTEMP2=202
!    SIGN=-1
!    CASE(4)
!    ITEMP1=108;ITEMP2=110;JTEMP1=311;JTEMP2=311
!    SIGN=-1
!    CASE(5)
!    ITEMP1=118;ITEMP2=118;JTEMP1=311;JTEMP2=316
!    SIGN=-1
!    END SELECT
                               
!     DO I=ITEMP1,ITEMP2
!       DO J=JTEMP1,JTEMP2
!         IF(LMASKDRY(LIJ(I,J))) THEN
!           Waterflow=SIGN*U(LIJ(I,J),1)*HP(LIJ(I,J))*DXU(LIJ(I,J))+SIGN*V(LIJ(I,J),1)*HP(LIJ(I,J))*DYV(LIJ(I,J))
!           Waterflowtot(LOC)=Waterflowtot(LOC)+Waterflow*ISHPRT*DELTAT
!         ENDIF
!       ENDDO
!     ENDDO
!ENDDO

DO  L=2,LA
    ! Sediment thickness
    TSET0T(L)=SUM(TSED0(1:KB,L)/BULKDENS(1:KB,L))
	TSEDT(L)=SUM(TSED(1:KB,L)/BULKDENS(1:KB,L))
    THCK(L)=TSEDT(L)-TSET0T(L)

    ! Calculate surface elevations and velocities for all cells and levels
    zeta(L)=(HP(L)+BELV(L))
ENDDO 

! Write velocity calibration data each call
WRITE(112,'(16F7.3)') time_efdc,zeta(LIJ(122,114)), &
    U(LIJ(122,114),1),V(LIJ(122,114),1),zeta(LIJ(45,29)),U(LIJ(45,29),1),V(LIJ(45,29),1), &
    zeta(LIJ(39,202)),U(LIJ(39,202),1),V(LIJ(39,202),1),zeta(LIJ(119,312)),U(LIJ(119,312),1), &
    V(LIJ(119,312),1),zeta(LIJ(130,349)),U(LIJ(130,349),1),V(LIJ(130,349),1)

! Write flow calibration data each call
!WRITE(113,'(6F7.3)')  time_efdc,Waterflowtot(1), &
!    Waterflowtot(2),Waterflowtot(3),Waterflowtot(4),Waterflowtot(5)

! Write tracer calibration data each call
WRITE(115,'(11F7.3)') time_efdc,SAL(LIJ(122,114),1), &
    DYE(LIJ(122,114),1),SAL(LIJ(45,29),1),DYE(LIJ(45,29),1),SAL(LIJ(39,202),1),DYE(LIJ(39,202),1), &
    SAL(LIJ(119,312),1),DYE(LIJ(119,312),1),SAL(LIJ(130,349),1),DYE(LIJ(130,349),1)

! If sediment is activated
IF(ISTRAN(6).EQ.1) THEN
    ! TSS calibration file with SEDZLJ (hard coded for 1 water layer (KC))
    DO K=1,NSCM
       WRITE(105,299)  time_efdc, K, SED(LIJ(122,114),1,K), &
        SED(LIJ(45,29),1,K), SED(LIJ(39,202),1,K), SED(LIJ(119,312),1,K), SED(LIJ(130,349),1,K)
    ENDDO

    ! Boundary calibration file with SEDZLJ
    WRITE(106, '(7F10.3)') time_efdc, zeta(LIJ(136,92)), &
     SAL(LIJ(136,92),1), SED(LIJ(136,92),1,1)+SED(LIJ(136,92),1,2), zeta(LIJ(77,8)), SAL(LIJ(77,8),1), &
     SED(LIJ(77,8),1,1)+SED(LIJ(77,8),1,2)

! If sediment is NOT activated, create files but fill with -7999 value
ELSE
    ! TSS calibration file
    WRITE(105,299)  time_efdc, 1, -7999.0, -7999.0, -7999.0, -7999.0, -7999.0

    ! Boundary calibration file
    WRITE(106, '(7F10.3)') time_efdc, zeta(LIJ(136,92)), &
        SAL(LIJ(136,92),1),-7999.0, zeta(LIJ(77,8)), SAL(LIJ(77,8),1),-7999.0
ENDIF

! Define not a number.  Only compatible with ifort
NAN=1.0/0.0

DO J=3,JC-2
	DO I=3,IC-2
	    IF(LIJ(I,J)>0) THEN
            IF(LMASKDRY(L).AND.HP(L).GT.0.3) THEN

                IF(vel_maxc(L).GT.vel_max(L).AND.vel_maxc(L).NE.NAN) THEN
                    vel_max(L)=vel_maxc(L)
                ENDIF

                IF(TAU(L).GT.tau_max(L).AND.TAU(L).NE.NAN) THEN
                    tau_max(L)=TAU(L)
                ENDIF

            ENDIF
	    ENDIF
    ENDDO
ENDDO

! Shear stress calibration file
WRITE(107,'(11F10.3)') time_efdc,TAU(LIJ(122,114)),tau_max(LIJ(122,114)),TAU(LIJ(45,29)),tau_max(LIJ(45,29)), &
    TAU(LIJ(39,202)),tau_max(LIJ(39,202)),TAU(LIJ(119,312)),tau_max(LIJ(119,312)), &
    TAU(LIJ(130,349)),tau_max(LIJ(130,349))

! Format for tss_cal.dat file
299 FORMAT(F7.3,2X,I1,F10.3,F10.3,F10.3,F10.3,F10.3)

RETURN
END SUBROUTINE MHS_OUT