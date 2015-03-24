C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE LTMT
C
C **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
C
C----------------------------------------------------------------------C
C
C CHANGE RECORD
C DATE MODIFIED   July 2014 Craig Jones
C
C----------------------------------------------------------------------C
C
C **  SUBROUTINE LTMT EXECUTES A LONG-TERM MASS TRANSPORT
C **  TIME INTEGRATION 
C
C**********************************************************************C
C
      USE GLOBAL 
C
C**********************************************************************C
C
      OPEN(99,FILE='RESTRAN.INP',STATUS='UNKNOWN')
C
C**********************************************************************C
C
C **  READ FIRST MEAN MASS TRANSPORT FIELD AND INITIALIZED NEXT
C
C----------------------------------------------------------------------C
C
      CALL RESTRAN
C
C
C
C **  ASSIGN VOLUMETRIC SOURCES AND SINKS
C
      do_110: DO NS=1,NQSER
        do_120: DO K=1,KC
           QSERT(K,NS)=QSRTLPP(K,NS)+QSRTLPN(K,NS)
        END DO do_120
      END DO do_110

C **  ASSIGN TOTAL WATER DEPTH
C
      do_130: DO L=2,LA
C
        H1P(L)=HP_OLD(L)
        HP(L)=HP_NEW(L)
C
      END DO do_130

        DO K=1,KS
            DO L=2,LA
                AB(L,K)=ABLPF(L,K)
            ENDDO
        ENDDO

C
C **  ASSIGN HORIZONTAL VOLUME FLUX
C
      do_160: DO K=1,KC
        do_170: DO L=2,LA
          UHDY2(L,K)=UHDY2LPF(L,K)
          VHDX2(L,K)=VHDX2LPF(L,K)
        END DO do_170
	END DO do_160

C     CALCULATE U & V
C
      do_165: DO K=1,KC
        do_175: DO L=2,LA
          LS=LSC(L)
          U(L,K)=UHDY2(L,K)/(DYU(L)*0.5*(HP(L)+HP(L-1)))
          V(L,K)=VHDX2(L,K)/(DXV(L)*0.5*(HP(L)+HP(LS)))
        END DO do_175
      END DO do_165

C--------------------------- -------------------------------------------C

C**********************************************************************C
C
C **  INITIALIZE COUNTERS
C
      NCTS=0
      N=0
C
C**********************************************************************C
C**********************************************************************C
C
C **  BEGIN THE TIME INTEGRATION LOOP
C
C----------------------------------------------------------------------C
C Start at first averaging period after beginning (NTSMMT)
C
      DO 1000 N=1,NTS
        NCTS=NCTS+1
C
C**********************************************************************C
C
C **  ADVANCE CONCENTRATION FIELDS
C
C----------------------------------------------------------------------C
C

      do_190: DO L=2,LA
            HPI(L)=1./HP(L)
            HU(L)=0.5*(HP(L)+HP(L-1))
            HUI(L)=1./HU(L)
            LS=LSC(L)
            HV(L)=0.5*(HP(L)+HP(LS))
            HVI(L)=1./HV(L)
      END DO do_190

      CALL CALCSER (2) !Updates time series at boundaries
      CALL CALCONC (2,0) ! Calculates advection of scalars.

      IF(ISHOW.GT.0) CALL SHOWVAL
C
C----------------------------------------------------------------------C
C
C **  CHECK RANGE OF SALINITY AND DYE CONCENTRATION
C
!      IF(ISMMC.EQ.1)THEN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(SAL(L,K).GT.SALMAX)THEN
       SALMAX=SAL(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(SAL(L,K).LT.SALMIN)THEN
       SALMIN=SAL(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
C
!      WRITE(6,6001)N
!      WRITE(6,6002)SALMAX,IMAX,JMAX,KMAX
!      WRITE(6,6003)SALMIN,IMIN,JMIN,KMIN
C
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
      DO L=2,LA
      IF(DYE(L,K).GT.SALMAX)THEN
       SALMAX=DYE(L,K)
       IMAX=IL(L)
       JMAX=JL(L)
       KMAX=K
      ENDIF
      IF(DYE(L,K).LT.SALMIN)THEN
       SALMIN=DYE(L,K)
       IMIN=IL(L)
       JMIN=JL(L)
       KMIN=K
      ENDIF
      ENDDO
      ENDDO
C
!      WRITE(6,6004)SALMAX,IMAX,JMAX,KMAX
!      WRITE(6,6005)SALMIN,IMIN,JMIN,KMIN
C
!      ENDIF
C
 6001 FORMAT('  N=',I10)
 6002 FORMAT('  SALMAX=',F14.4,5X,'I,J,K=',(3I10))
 6003 FORMAT('  SALMIN=',F14.4,5X,'I,J,K=',(3I10))
 6004 FORMAT('  DYEMAX=',F14.4,5X,'I,J,K=',(3I10))
 6005 FORMAT('  DYEMIN=',F14.4,5X,'I,J,K=',(3I10))
C
C**********************************************************************C
C
C **  Update Scalar
C
C----------------------------------------------------------------------C
        do_240: DO K=1,KC
          do_250: DO L=2,LA
            SAL1(L,K)=SAL(L,K)
            TEM1(L,K)=TEM(L,K)
            DYE1(L,K)=DYE(L,K)
          END DO do_250
        END DO do_240



C**********************************************************************C
C
C **  ADVANCE MMT FIELDS
C
C----------------------------------------------------------------------C
C
      IF(ISSSMMT.EQ.0.AND.N.LT.NTS-1)THEN ! Normal LTMT transport
C
      IF(NCTS+1.GT.NTSMMT)THEN
      NCTS=NCTS-NTSMMT
C
C
          CALL RESTRAN
C
C
C **        UPDATE TOTAL WATER DEPTH FOR CURRENT & NEXT TIME STEP
C
            do_300: DO L=2,LA
C
              H1P(L)=HP(L)
              HP(L)=HP_OLD(L)+(HP_NEW(L)-HP_OLD(L))*FLOAT(NCTS+NFSTP)/
     +                                                     FLOAT(NTSMMT)
C 
	      END DO do_300
C
C **        UPDATE VERTICAL DIFFUSIVITY & VOLUME FLUX FOR NEXT AVERAGE PERIOD
C
            do_310: DO K=1,KS
              do_320: DO L=2,LA
                AB(L,K)=SQRT(AB(L,K)*ABLPF(L,K))   ! SQRT FILTER
              END DO do_320
	      END DO do_310
C
C **        UPDATE HORIZONTAL VOLUME FLUX FOR NEXT AVERAGE PERIOD
C
            do_330: DO K=1,KC
              do_340: DO L=2,LA
                UHDY2(L,K)=UHDY2LPF(L,K)
                VHDX2(L,K)=VHDX2LPF(L,K)
              END DO do_340
            END DO do_330
C
          ELSE
C
C **        UPDATE TOTAL WATER DEPTH FOR CURRENT & NEXT TIME STEP
C
            do_350: DO L=2,LA
C
              H1P(L)=HP(L)
              HP(L)=HP_OLD(L)+(HP_NEW(L)-HP_OLD(L))*FLOAT(NCTS+NFSTP)/
     +                                                     FLOAT(NTSMMT)
C
            END DO do_350
C
C     
C         CALCULATE U & V
C
          do_355: DO K=1,KC
            do_365: DO L=2,LA
              LS=LSC(L)
              U(L,K)=UHDY2(L,K)/(DYU(L)*0.5*(HP(L)+HP(L-1)))
              V(L,K)=VHDX2(L,K)/(DXV(L)*0.5*(HP(L)+HP(LS)))
            END DO do_365
	    END DO do_355  
C
C
C
      ENDIF
C
      ENDIF
C
C**********************************************************************C
C
 1000 CONTINUE
C
C**********************************************************************C
C**********************************************************************C
C
C **  WRITE GRAPHICS FILES FOR RESIDUAL VARIABLES
C

      CLOSE(99) 
C
C**********************************************************************C
C
      RETURN
      END
