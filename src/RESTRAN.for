      SUBROUTINE RESTRAN  
C  
C CHANGE RECORD  
C **  SUBROUTINE RESTRAN READS A RESIDUAL TRANSPORT FILE  
C  
      USE GLOBAL
     
      IF(NTSMMT.LT.NTSPTC)THEN  
        DO L=2,LA  
          READ(99)HMP(L),HP_OLD(L),HP_NEW(L)  
          READ(99)(UHDY2LPF(L,K),K=1,KC)
          READ(99)(VHDX2LPF(L,K),K=1,KC) 
          READ(99)(AHULPF(L,K),K=1,KC)  
          READ(99)(AHVLPF(L,K),K=1,KC)  
          READ(99)(SALLPF(L,K),K=1,KC)  
          READ(99)(ABLPF(L,K),K=1,KS)  
          READ(99)(ABEFF(L,K),K=1,KS)  
        ENDDO  

        do_120: DO NS=1,NQSER
          READ(99) (QSRTLPP(K,NS),K=1,KC)
          READ(99) (QSRTLPN(K,NS),K=1,KC)
        END DO do_120

        do_125: DO L=2,LA
          READ(99) TAULPF(L)
        END DO do_125

        do_126: DO L=2,LA
          TAU(L)=TAULPF(L)
        END DO do_126

      ELSE  
        DO L=2,LA  
          READ(99)HMP(L),HLPF(L),QSUMELPF(L)  
          READ(99)(UHLPF(L,K),K=1,KC)  
          READ(99)(VHLPF(L,K),K=1,KC)  
          READ(99)(VPZ(L,K),K=1,KC)  
          READ(99)(AHULPF(L,K),K=1,KC)  
          READ(99)(AHVLPF(L,K),K=1,KC)  
          READ(99)(SALLPF(L,K),K=1,KC)  
          READ(99)(VPX(L,K),K=1,KS)  
          READ(99)(VPY(L,K),K=1,KS)  
          READ(99)(ABLPF(L,K),K=1,KS)  
C  
C      READ(99,907)(ABEFF(L,K),K=1,KS)  
C  
        ENDDO  
      ENDIF  
      DO K=1,KC  
        DO L=2,LA  
          AHULPF(L,K)=AHULPF(L,K)+AHO  
          AHVLPF(L,K)=AHVLPF(L,K)+AHO  
        ENDDO  
      ENDDO  
      DO K=1,KC  
        DO L=2,LA  
          AH(L,K)=0.25*(AHULPF(L,K)+AHULPF(L+1   ,K)  
     &        +AHVLPF(L,K)+AHVLPF(LNC(L),K))  
        ENDDO  
      ENDDO  
C      IF(NTSMMT.LT.NTSPTC.OR.ISLTMT.EQ.2)THEN  
      IF(NTSMMT.LT.NTSPTC)THEN  
        DO K=1,KC  
          DO L=2,LA  
            VPZ(L,K)=0.  
          ENDDO  
        ENDDO  
        DO K=1,KS  
          DO L=2,LA  
            VPX(L,K)=0.  
            VPY(L,K)=0.  
          ENDDO  
        ENDDO  
      ENDIF  
  907 FORMAT(12E12.4)  
      RETURN  
      END  

