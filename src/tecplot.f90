SUBROUTINE TECPLOT

!
! Tecplot customized output
! Calculate shear stress here even for non-sediment cases
!
! CALL SEDZLJ_SHEAR
!
! REVISION DATE:  May, 2014
! Craig Jones and Scott James
!***************************************************************
!*********************************************
USE GLOBAL
IMPLICIT NONE

INTEGER::I,J,LN,L,K,KK,ITEMPMSK,LOC
INTEGER::ITEMP1,ITEMP2,JTEMP1,JTEMP2,SIGN
INTEGER,DIMENSION(LCM)::NDAVG
REAL::UTMP,VTMP,TEMPMSK,FLUXLOAD,NAN
REAL::UTMPA,VTMPA,TEMPSTINC,LONTEMP,LATTEMP,WVTMP
REAL,DIMENSION(KC)::CTEMP1
REAL,DIMENSION(LCM)::UTMPS,VTMPS,VMAG,SURFEL
REAL,DIMENSION(LCM)::DMAX,DAVG,VMAGC,CSMAX
INTEGER,SAVE::nstep
REAL::deltat,CPGV
LOGICAL,SAVE::FIRSTTIME=.FALSE.	

!REAL Waterflow
!REAL,DIMENSION(8)::Waterflowtot

REAL,DIMENSION(LCM)::VMAX,TAUMAX,TAUAVG

! Create files if first call
IF(.NOT.FIRSTTIME)THEN
	
!  This opens the Tecplot output file
!        IF(MAXVAL(MVEGL(2:LA))>90)THEN !MHK devices exist
!			OPEN(UNIT=222,FILE='powerout.dat')
!			FORALL(I=1:TCOUNT)ICOUNT(I)=I
!			WRITE(222,'("TURBINE",100(I6,6X))')(ICOUNT(I),I=1,TCOUNT)
!			WRITE(222,'(7X,100(3X,I3,3X,I3))')((IJLTURB(I,1),IJLTURB(I,2)),I=1,TCOUNT)
!		ENDIF
    OPEN (UNIT=111,FILE='tecplot2d.dat',STATUS='REPLACE')
	WRITE(111,'(A30)')'TITLE = "EFDC 2D Tecplot Data"'

    ! Header for Sandia Coastal variables
    WRITE(111,*)'VARIABLES= "I","J","X","Y","U","V","HP","TAU","D50","THCK","WvHt","DYE"'
    
    ! Header for BCSA variables
    !WRITE(111,*)'VARIABLES= "I","J","X","Y","TAU","TAUAVG","VMAX"'

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

    ! Initialize variable for first time step
    DMAX=0.1
    CMAX=0.1
    VMAGC=0.0
    VMAX=0.0
    TAUMAX=0.0
    NDAVG=1
    DAVG=0.0
    TAUAVG=0.0

	FIRSTTIME=.TRUE.
ENDIF

! Timing parameters
deltat=tidalp/float(ntsptc)
nstep=nstep+1

! Tecplot customized output
TEMPMSK=-1.0
ITEMPMSK=-1

!  Calculate shear stress here even for non-sediment cases
CALL SEDZLJ_SHEAR

! DETERMINE INITIAL SED. THICKNESS
! CALC. FINAL SEDIMENT THICKNESS, INITIAL SEDIMENT-WATER INTERFACE
! IS AT ZERO
FORALL(L=2:LA)
	TSET0T(L)=SUM(TSED0(1:KB,L)/BULKDENS(1:KB,L))
	TSEDT(L)=SUM(TSED(1:KB,L)/BULKDENS(1:KB,L))
ENDFORALL
FORALL(L=2:LA)THCK(L)=TSEDT(L)-TSET0T(L)

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

!*********************************************
! Write tecplot2d.dat header at each call
WRITE(111,*)'ZONE T="',tbegin+float(nstep-1)*deltat*float(ishprt)/86400.0,'" I= ' ,IC-4,' J= ' ,JC-4,' F=POINT'

IF(STINC.LT.1)THEN
	TEMPSTINC=1
ELSE
	TEMPSTINC=STINC
ENDIF

DO  L=2,LA
	I=IL(L)
	J=JL(L)

	IF(IWRSP(1)<1) NCORENO(I,J)=0
	    CBLTOT(L)=1.0E6*SUM(CBL(1,L,1:NSCM))
	    FORALL(KK=1:KC)CAVG(L,KK)=SUM(SED(L,KK,1:NSCM))
	    FORALL(K=1:KC-1)HEIGHT(I,J,K)=-(HP(L)-HP(L)*Z(K))
	    HEIGHT(I,J,KC)=0.0
	    FORALL(K=1:KC)
		    CTEMP1(K)=SUM(SED(L,K,1:NSCM))
	    ENDFORALL
	    CAVGT(L)=SUM(CTEMP1(1:KC)*DZC(1:KC))
	    FORALL(KS=1:NSCM)
    	    CAVGS(L,KS)=SUM(SED(L,1:KC,KS)*DZC(1:KC))
	    ENDFORALL
ENDDO 

! Calculate surface elevations and velocities for all cells and levels
DO L=2,LA
	LN=LNC(L)
    SURFEL(L)=(HP(L)+BELV(L))
    VMAG(L)=SQRT(U(L,1)**2+V(L,1)**2) 

	DO K=1,KC
!		UTMPS=U(LIJ(I,J),K) ! m/s
!		VTMPS=V(LIJ(I,J),K)  
		UTMPS(L)=0.5*STCUV(L)*(RSSBCE(L)*U(L+1,K)+RSSBCW(L)*U(L,K))  ! m/s
		VTMPS(L)=0.5*STCUV(L)*(RSSBCN(L)*V(LN ,K)+RSSBCS(L)*V(L,K)) 
!		UTECPLOT(L,K)=CUE(L)*UTMPS+CVE(L)*VTMPS  
!		VTECPLOT(L,K)=CUN(L)*UTMPS+CVN(L)*VTMPS
	ENDDO
ENDDO

! Write velocity calibration data each call
WRITE(112,'(16F7.3)') tbegin+float(nstep-1)*deltat*float(ishprt)/86400.0,SURFEL(LIJ(122,114)), &
    U(LIJ(122,114),1),V(LIJ(122,114),1),SURFEL(LIJ(45,29)),U(LIJ(45,29),1),V(LIJ(45,29),1), &
    SURFEL(LIJ(39,202)),U(LIJ(39,202),1),V(LIJ(39,202),1),SURFEL(LIJ(119,312)),U(LIJ(119,312),1), &
    V(LIJ(119,312),1),SURFEL(LIJ(130,349)),U(LIJ(130,349),1),V(LIJ(130,349),1)

! Write flow calibration data each call
!WRITE(113,'(6F7.3)')  tbegin+float(nstep-1)*deltat*float(ishprt)/86400.0,Waterflowtot(1), &
!    Waterflowtot(2),Waterflowtot(3),Waterflowtot(4),Waterflowtot(5)

! Write tracer calibration data each call
WRITE(115,'(11F7.3)') tbegin+float(nstep-1)*deltat*float(ishprt)/86400.0,SAL(LIJ(122,114),1), &
    DYE(LIJ(122,114),1),SAL(LIJ(45,29),1),DYE(LIJ(45,29),1),SAL(LIJ(39,202),1),DYE(LIJ(39,202),1), &
    SAL(LIJ(119,312),1),DYE(LIJ(119,312),1),SAL(LIJ(130,349),1),DYE(LIJ(130,349),1)

! Write TSS calibration data each call (hard coded for single water layer (KC))
DO K=1,NSCM
    WRITE(105, 299) tbegin+float(nstep-1)*deltat*float(ishprt)/86400.0, K, SED(LIJ(122,114),1,K), &
        SED(LIJ(45,29),1,K), SED(LIJ(39,202),1,K), SED(LIJ(119,312),1,K), SED(LIJ(130,349),1,K)
ENDDO

! 2D ouput of tecplot2d.out for all cells
NAN=1.0/0.0
TAUAVG=0.0

DO J=3,JC-2
	DO I=3,IC-2
		
	    IF(LIJ(I,J)>0) THEN
            L=LIJ(I,J)	
            UTMPA=0.0
            VTMPA=0.0
            WVTMP=SQRT(8/G*WVENEP(L))
            VMAGC(L)=SQRT(UTMPA**2+VTMPA**2)

            ! Calculate velocities for each water layer
	        DO K=1,KC
		        UTMPS=U(LIJ(I,J),K) ! m/s
		        VTMPS=V(LIJ(I,J),K) 
		        UTMPA=UTMPS(LIJ(I,J))*DZC(K)+UTMPA
		        VTMPA=VTMPS(LIJ(I,J))*DZC(K)+VTMPA
	        ENDDO	

            ! 
            IF(LMASKDRY(L).AND.HP(L).GT.0.3) THEN
                CSMAX(L)=MAX(CSMAX(L),CAVGT(L))
                DMAX(L)=MAX(DMAX(L),DYE(L,1))

                DTOT(L)=DTOT(L)+DYE(L,1)
                DAVG(L)=DTOT(L)/NDAVG(L)

                TAUTTOT(L)=TAUTTOT(L)+TAU(L)
                TAUAVG(L)=TAUTTOT(L)/NDAVG(L)

                IF(I.EQ.130.AND.J.EQ.366)THEN
                    TAUAVG(L)=0.0
                ENDIF

                NDAVG(L)=NDAVG(L)+1
                UTMPA=0.0
                VTMPA=0.0

	            DO K=1,KC
		            UTMPS=U(L,K) ! m/s
		            VTMPS=V(L,K)
                    UTMPA=UTMPS(L)*DZC(K)+UTMPA
		            VTMPA=VTMPS(L)*DZC(K)+VTMPA
	            ENDDO

                VMAGC(L)=SQRT(UTMPA**2+VTMPA**2)

                IF(VMAGC(L).GT.VMAX(L).AND.VMAGC(L).NE.NAN) THEN
                    VMAX(L)=VMAGC(L)
                ENDIF

                IF(TAU(L).GT.TAUMAX(L).AND.TAU(L).NE.NAN) THEN
                    TAUMAX(L)=TAU(L)
                ENDIF
            ENDIF
            
            ! Write tecplot data line for active cell    
            WRITE(111,'(I4,1X,I4,1X,10E17.7)')I,J,DLON(LIJ(I,J)),DLAT(LIJ(I,J)), &
			    UTMPA,VTMPA,HP(LIJ(I,J)),TAU(LIJ(I,J)),D50AVG(LIJ(I,J)),THCK(LIJ(I,J)), &
                WVTMP,TAUMAX(L) !Sandia Coastal

!           WRITE(111,'(I4,1X,I4,1X,10E17.7)')I,J,DLON(LIJ(I,J)),DLAT(LIJ(I,J)), &
!				TAUMAX(L),TAUAVG(L),CSMAX(L) !BCSA Average variables
                
        ELSE

            ! Write tecplot data line for inactive cell
			WRITE(111,'(I4,1X,I4,1X,12E13.4)')I,J,TEMPMSK,TEMPMSK,TEMPMSK,TEMPMSK,TEMPMSK, &
				TEMPMSK,0,0,TEMPMSK,TEMPMSK ! Sandia Coastal

!			WRITE(111,'(I4,1X,I4,1X,12E13.4)')I,J,TEMPMSK,TEMPMSK, &
!			    TEMPMSK,TEMPMSK,TEMPMSK ! BCSA
	    ENDIF
    ENDDO
ENDDO

! Format for tss_cal.dat file
299 FORMAT(F7.3,2X,I1,F7.3,F7.3,F7.3,F7.3,F7.3)

RETURN
END SUBROUTINE TECPLOT
