SUBROUTINE MHS_OUT

! Output MHS station data in calibration files
!***************************************************************
!***************************************************************

USE GLOBAL
IMPLICIT NONE

INTEGER::I,J,L,K,LOC
INTEGER::ITEMP1,ITEMP2,JTEMP1,JTEMP2
REAL,DIMENSION(LCM)::zeta
REAL*8::time_efdc
REAL::tss_flux_u_tmp, tss_flux_v_tmp, flow_u_tmp, flow_v_tmp
REAL,DIMENSION(8)::tss_flux_u, tss_flux_v, flow_u, flow_v
LOGICAL,SAVE::FIRSTTIME=.FALSE.	

! Create files if first call
IF(.NOT.FIRSTTIME)THEN

! Initialize max. shear stress value to be computed in s_sedzlj.f90
TAUMAX = 0.0

    ! Velocity calibration file with surface elevation and velocity components at all MHS stations
    OPEN (UNIT=112,FILE='vel_cal.dat', STATUS='REPLACE')
	
	! Flow rate
	OPEN (UNIT=109, FILE='flow.dat', STATUS='REPLACE')

    ! Tracer calibration file with salinity and dye
    IF (ISTRAN(1).EQ.1.OR.ISTRAN(3).EQ.1) OPEN (UNIT=115,FILE='tracer_cal.dat', STATUS='REPLACE')

    IF (ISTRAN(6).EQ.1) THEN
        ! TSS calibration file
        OPEN (UNIT=105, FILE='tss_cal.dat', STATUS='REPLACE')
    
        ! Bed thickness calibration file
        OPEN (UNIT=106, FILE='thick_cal.dat', STATUS='REPLACE')
    
        ! TSS flux file
        OPEN (UNIT=108, FILE='tss_flux.dat', STATUS='REPLACE')
    ENDIF

    ! Shear stress calibration file
    OPEN (UNIT=107, FILE='shear_cal.dat', STATUS='REPLACE')

    FIRSTTIME=.TRUE.
ENDIF

! EFDC time parameter
time_efdc=DT*FLOAT(N)+TCON*TBEGIN 
time_efdc=time_efdc/86400.

! Calculate flow rate and TSS flux at MHS01, MHS02, MHS05, MHS06, MHS07, and the LBC-BCC connector
DO LOC=1,6
    ! Initialize temporary variables
    flow_u(LOC)=0.0
    flow_v(LOC)=0.0
    tss_flux_u(LOC)=0.0
    tss_flux_v(LOC)=0.0

    SELECT CASE (LOC)
        CASE(1)
		    ! MHS01 flux line
            ITEMP1=116;ITEMP2=129;JTEMP1=114;JTEMP2=114

        CASE(2)
		    ! MHS02 flux line
            ITEMP1=43;ITEMP2=46;JTEMP1=23;JTEMP2=23

        CASE(3)
		    ! MHS05 flux line
            ITEMP1=35;ITEMP2=41;JTEMP1=203;JTEMP2=203

        CASE(4)
		    ! MHS06 flux line
            ITEMP1=119;ITEMP2=119;JTEMP1=311;JTEMP2=314
			
        CASE(5)
		    ! MHS07 flux line
		    ITEMP1=128;ITEMP2=131;JTEMP1=349;JTEMP2=349
			
        CASE(6)
		    ! LBC-BCC connector flux line
		    ITEMP1=40;ITEMP2=40;JTEMP1=179;JTEMP2=179
    END SELECT

    DO I=ITEMP1,ITEMP2
        DO J=JTEMP1,JTEMP2
            IF(LMASKDRY(LIJ(I,J))) THEN
                ! Compute flow rate
                flow_u_tmp=U(LIJ(I,J),1)*HP(LIJ(I,J))*DXU(LIJ(I,J))
                flow_u(LOC)=flow_u(LOC) + flow_u_tmp
                flow_v_tmp=V(LIJ(I,J),1)*HP(LIJ(I,J))*DYV(LIJ(I,J))
                flow_v(LOC)=flow_v(LOC) + flow_v_tmp

                ! Only compute TSS flux if SED ISTRAN flag is activated in efdc.inp
                IF (ISTRAN(6).EQ.1) THEN
                    DO K=1,NSCM
                        ! Compute TSS flux for all size classes
                        tss_flux_u_tmp=U(LIJ(I,J),1)*HP(LIJ(I,J))*DXU(LIJ(I,J))*SED(LIJ(I,J),1,K)
                        tss_flux_u(LOC)=tss_flux_u(LOC)+tss_flux_u_tmp
                        tss_flux_v_tmp=V(LIJ(I,J),1)*HP(LIJ(I,J))*DYV(LIJ(I,J))*SED(LIJ(I,J),1,K)
                        tss_flux_v(LOC)=tss_flux_v(LOC)+tss_flux_v_tmp                    
                    ENDDO
                ENDIF
            ENDIF
        ENDDO
    ENDDO
ENDDO


DO  L=2,LA
    ! Sediment thickness
    TSET0T(L)=SUM(TSED0(1:KB,L)/BULKDENS(1:KB,L))
    TSEDT(L)=SUM(TSED(1:KB,L)/BULKDENS(1:KB,L))
    THCK(L)=TSEDT(L)-TSET0T(L)

    ! Calculate surface elevations and velocities for all active cells
    zeta(L)=(HP(L)+BELV(L))
ENDDO 

! Write velocity calibration file
WRITE(112,298) time_efdc,zeta(LIJ(122,114)), &
    U(LIJ(122,114),1),V(LIJ(122,114),1),zeta(LIJ(45,29)),U(LIJ(45,29),1),V(LIJ(45,29),1), &
    zeta(LIJ(39,202)),U(LIJ(39,202),1),V(LIJ(39,202),1),zeta(LIJ(119,312)),U(LIJ(119,312),1), &
    V(LIJ(119,312),1),zeta(LIJ(130,349)),U(LIJ(130,349),1),V(LIJ(130,349),1)
    
    ! Write flow rate file
    WRITE(109,299)  time_efdc, flow_u(1),flow_v(1), flow_u(2),flow_v(2), &
    flow_u(3),flow_v(3), flow_u(4),flow_v(4), flow_u(5),flow_v(5), &
	flow_u(6),flow_v(6)

IF (ISTRAN(1).EQ.1.OR.ISTRAN(3).EQ.1) THEN
    ! Write tracer calibration file
    WRITE(115,300) time_efdc,SAL(LIJ(122,114),1), &
    DYE(LIJ(122,114),1),SAL(LIJ(45,29),1),DYE(LIJ(45,29),1),SAL(LIJ(39,202),1),DYE(LIJ(39,202),1), &
    SAL(LIJ(119,312),1),DYE(LIJ(119,312),1),SAL(LIJ(130,349),1),DYE(LIJ(130,349),1)
ENDIF

! If sediment is activated
IF(ISTRAN(6).EQ.1) THEN
    ! TSS calibration file with SEDZLJ (hard coded for 1 water layer (KC))
    DO K=1,NSCM
        WRITE(105,301)  time_efdc, K, SED(LIJ(122,114),1,K), &
        SED(LIJ(45,29),1,K), SED(LIJ(39,202),1,K), SED(LIJ(119,312),1,K), SED(LIJ(130,349),1,K)
    ENDDO

    ! Write thickness file for MHS01,MHS02,MHS05,MHS06, and MHS07
    WRITE(106,302)  time_efdc, THCK(LIJ(122,114)), THCK(LIJ(45,29)), &
    THCK(LIJ(39,202)), THCK(LIJ(119,312)), THCK(LIJ(130,349))

    ! Write TSS flux file if SED is activated
    WRITE(108,303)  time_efdc, tss_flux_u(1),tss_flux_v(1), tss_flux_u(2),tss_flux_v(2), &
    tss_flux_u(3),tss_flux_v(3), tss_flux_u(4),tss_flux_v(4), tss_flux_u(5),tss_flux_v(5), &
	tss_flux_u(6),tss_flux_v(6)
	
ENDIF

! Shear stress calibration file
WRITE(107,304) time_efdc,TAU(LIJ(122,114)),TAU(LIJ(45,29)), &
    TAU(LIJ(39,202)),TAU(LIJ(119,312)),TAU(LIJ(130,349))
	
! Write formats 
298 FORMAT(F9.5,15F9.3)
299 FORMAT(F9.5,12F10.3)
300 FORMAT(F9.5,10F9.3)
301 FORMAT(F9.5,2X,I1,5F12.3)
302 FORMAT(F9.5,5F9.3)
303 FORMAT(F9.5,12F15.3)
304 FORMAT(F9.5,5F9.3)

FLUSH(112)
FLUSH(109)
IF (ISTRAN(1).EQ.1.OR.ISTRAN(3).EQ.1) FLUSH(115)

IF (ISTRAN(6).EQ.1) THEN
    FLUSH(105)
    FLUSH(106)
    FLUSH(108)
ENDIF

FLUSH(107)

RETURN
END SUBROUTINE MHS_OUT
