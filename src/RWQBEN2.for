      SUBROUTINE RWQBEN2 (TIMTMP)  
C  
C M. MORTON 01/30/98: CHANGED CODE TO ALLOW FOR TEMPORALLY  
C   VARYING BENTHIC FLUXES IN THE BENFN FILE.  PREVIOUS VERSION ONLY  
C   PROVIDED SPATIALLY VARYING FLUX (NO PROVISION FOR TIME VARYING).  
C CHANGE RECORD  
C READ IN SPATIALLY AND/OR TEMPORALLY VARYING PARAMETERS FOR BENTHIC  
C FLUXES OF PO4D, NH4, NO3, SAD, COD, O2  
C FORMAT OF BENFN FILE IS:  
C    TITLE 1  
C    TITLE 2  
C    TITLE 3  
C  270.00000  <-- DAY AT WHICH FOLLOWING FLUXES BECOME ACTIVE  
C  350.00000  <-- DAY AT WHICH FOLLOWING FLUXES BECOME ACTIVE  
C 9999.99999 <-- ENTER LARGE DAY AT END OF FILE  
C  
      USE GLOBAL  

      CHARACTER TITLE(3)*79, CCMRM*1  
      INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)::IZONE  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::XBFCOD  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::XBFNH4  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::XBFNO3  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::XBFO2  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::XBFPO4D  
      REAL,SAVE,ALLOCATABLE,DIMENSION(:)::XBFSAD  

      IF(.NOT.ALLOCATED(IZONE  ))THEN
		ALLOCATE(IZONE(NSMZM))
		ALLOCATE(XBFCOD(NSMZM))
		ALLOCATE(XBFNH4(NSMZM))
		ALLOCATE(XBFNO3(NSMZM))
		ALLOCATE(XBFO2(NSMZM))
		ALLOCATE(XBFPO4D(NSMZM))
		ALLOCATE(XBFSAD(NSMZM))
	    IZONE=0
	    XBFCOD=0.0 
	    XBFNH4=0.0 
	    XBFNO3=0.0 
	    XBFO2=0.0 
	    XBFPO4D=0.0 
	    XBSFAD=0.0 
	ENDIF
C  
      OPEN(1,FILE=BENFN,STATUS='UNKNOWN')  
      OPEN(2,FILE='WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')  
C  
C SKIP OVER THREE HEADER RECORDS:  
C  
      READ(1,50) (TITLE(M),M=1,3)  
      WRITE(2,999)  
      WRITE(2,50) (TITLE(M),M=1,3)  
C  
C SKIP OVER ALL COMMENT CARDS AT BEGINNING OF FILE:  
C  
      REWIND(1)  
      CCMRM = '#'  
      CALL SKIPCOMM(1, CCMRM)  
      READ(1, *) IBENZ  
      WRITE(2, 65) TIMTMP, IBENZ  
   65 FORMAT(' * BENTHIC FLUXES AT     ', F10.5,' DAYS OF MODEL RUN',/,  
     &    ' NUMBER OF BENTHIC FLUX ZONES = ', I4)  
C  
C SEQUENTIALLY READ THROUGH BENTHIC FLUX FILE UNTIL THE APPROPRIATE  
C TIME IS FOUND:  
C   BDAY   = CURRENT DAY AT WHICH BENTHIC FLUX IS IN EFFECT  
C   BENDAY = NEXT DAY AT WHICH BENTHIC FLUX CHANGES (PASSED TO MAIN PROG  
C  
   10 READ(1, *, END=15) BENDAY  
      IF(BENDAY .GT. TIMTMP) GOTO 20  
      BDAY = BENDAY  
      DO I=1,IBENZ  
        READ(1,*,END=15) MM, XBFPO4D(MM), XBFNH4(MM), XBFNO3(MM),  
     &      XBFSAD(MM), XBFCOD(MM), XBFO2(MM)  
        IZONE(I) = MM  
      ENDDO  
      GOTO 10  
C  
C UNEXPECTED END-OF-FILE ENCOUNTERED:  
C  
   15 WRITE(2,16) BENFN  
   16 FORMAT(//,' ************* WARNING *************',/,  
     &    ' END-OF-FILE ENCOUNTERED IN FILE: ', A20,/,/  
     &    ' BENTHIC FLUXES SET TO VALUES CORRESPONDING ',  
     &    ' TO LAST DAY IN FILE.',/)  
      BENDAY=(TCON*TBEGIN + NTC*TIDALP)/86400.0  ! *** PMC SINGLE LINE
   20 CONTINUE  
      WRITE(2, 48) BDAY  
   48 FORMAT(/,' DAY IN BENTHIC FLUX FILE: ',F10.5,/,  
     &    '    ZONE    FPO4    FNH4    FNO3    FSAD    FCOD    FSOD')  
      DO I=1,IBENZ  
        MM = IZONE(I)  
        WRITE(2,51) MM,XBFPO4D(MM),XBFNH4(MM),XBFNO3(MM),XBFSAD(MM),  
     &      XBFCOD(MM),XBFO2(MM)  
      ENDDO  
C  
C DETERMINE BENTHIC FLUX FOR EACH CELL (L) BY INTERPOLATING BETWEEN  
C THE MUD AND SAND FLUXES.  XBENMUD(L) IS THE PERCENT MUD FOR EACH  
C CELL.  
C  
      DO L=2,LA  
        IZM = IBENMAP(L,1)  
        IZS = IBENMAP(L,2)  
        XM = XBENMUD(L)  
        WQBFPO4D(L) = XM*XBFPO4D(IZM) + (1.0-XM)*XBFPO4D(IZS)  
        WQBFNH4(L)  = XM*XBFNH4(IZM)  + (1.0-XM)*XBFNH4(IZS)  
        WQBFNO3(L)  = XM*XBFNO3(IZM)  + (1.0-XM)*XBFNO3(IZS)  
        WQBFSAD(L)  = XM*XBFSAD(IZM)  + (1.0-XM)*XBFSAD(IZS)  
        WQBFCOD(L)  = XM*XBFCOD(IZM)  + (1.0-XM)*XBFCOD(IZS)  
        WQBFO2(L)   = XM*XBFO2(IZM)   + (1.0-XM)*XBFO2(IZS)  
      ENDDO  
      CLOSE(1)  
      CLOSE(2)  
  999 FORMAT(1X)  
   50 FORMAT(A79)  
   51 FORMAT(I8, 10F8.3)  
   52 FORMAT(I7, 1X, A3)  
   60 FORMAT(/, A24, I5, A24)  
      RETURN  
      END  

