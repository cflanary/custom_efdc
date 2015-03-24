      FUNCTION CSNDEQC(SNDDIA,SSG,WS,TAUR,TAUB,D50,SIGPHI,  
     &    SNDDMX,IOPT,ISNDAL)  
      INTEGER::IOPT,ISNDAL
      REAL::CSNDEQC,REY,SNDDIA,SSG,DFAC,D50,RLAM,SIGPHI,USTAR,TAUB
      REAL::VAL,WS,TMP,TAUR,TAURS,REY3,RATIO,SNDDMX
C CHANGE RECORD  
C  
C  
C **  CALCULATES NEAR BED REFERENCE CONCENTRATION FOR NONCOHESIVE  
C **  SEDIMENT  
C **  IOPT=1  BASED ON  
C **  
C **  GARCIA, M., AND G. PARKER, 1991: ENTRAINMENT OF BED SEDIMENT  
C **  INTO SUSPENSION, J. HYDRAULIC ENGINEERING, 117, 414-435.  
C  
      IF(IOPT.EQ.1)THEN  
        REY=1.E6*SNDDIA*SQRT( 9.8*(SSG-1.)*SNDDIA )  
        REY=REY**0.6  
        DFAC=1.  
        IF(ISNDAL.GE.1) DFAC=(SNDDIA/D50)**0.2  
        RLAM=1.-0.29*SIGPHI  
        USTAR=SQRT(TAUB)  
        VAL=DFAC*RLAM*REY*USTAR/WS  
        VAL=1.3E-7*(VAL**5)  
        TMP=VAL/(1+3.33*VAL)  
        CSNDEQC=1.E6*SSG*TMP  
        IF(USTAR.LT.WS) CSNDEQC=0.  
      ENDIF  
C  
C **  IOPT=2  BASED ON  
C **  
C **  SMITH, J. D., AND S. R. MCLEAN, 1977: SPATIALLY AVERAGED FLOW  
C **  OVER A WAVY SURFACE, J. GEOPHYSICAL RESEARCH, 82, 1735-1746.  
C  
      IF(IOPT.EQ.2)THEN  
        VAL=2.4E-3*( (TAUB/TAUR)-1. )  
        VAL=MAX(VAL,0.)  
        TMP=0.65*VAL/(1.+VAL)  
        CSNDEQC=1.E6*SSG*TMP  
        USTAR=SQRT(TAUB)  
        IF(USTAR.LT.WS) CSNDEQC=0.  
      ENDIF  
C  
C **  IOPT=3  BASED ON  
C **  
C **  VAN RIJN, L. C., 1984: SEDIMENT TRANSPORT, PART II: SUSPENDED  
C **  LOAD TRANSPORT, J. HYDRAULIC ENGINEERING, 110, 1623-1641.  
C  
      IF(IOPT.EQ.3)THEN  
        REY=1.E4*SNDDIA*( (9.8*(SSG-1.))**0.333 )  
        IF(REY.LE.10.) TAURS=(4.*WS/REY)**2  
        IF(REY.GT.10.) TAURS=0.016*WS*WS  
        REY3=REY**0.3  
        VAL=(TAUB/TAURS)-1.  
C  
C        VAL=(TAUB/TAUR)-1.  
C  
        VAL=MAX(VAL,0.)  
        VAL=VAL**1.5  
        RATIO=SNDDIA/(3.*SNDDMX)  
        TMP=0.015*RATIO*VAL/REY3  
        CSNDEQC=1.E6*SSG*TMP  
        USTAR=SQRT(TAUB)  
        IF(USTAR.LT.WS) CSNDEQC=0.  
      ENDIF  
  600 FORMAT(10E12.4)  
      RETURN  
      END  

