      FUNCTION FSTRSE(VOID,BMECH1,BMECH2,BMECH3)  
      REAL::FSTRSE,BMECH1,TMP,VOID,BMECH2,BMECH3,FSTRSELOG
C  ADDED STANDARD EXPONENTIAL FORM CONSTITUTIVE RELATIONSHIP  
C CHANGE RECORD  
C  
C  
C **  FSTRSE IS WATER SPECIFIC WEIGHT NORMALIZED EFFECTIVE STRESS  
C  
      IF(BMECH1.GT.0.0)THEN  
        TMP=-(VOID-BMECH2)/BMECH3  
        FSTRSE=BMECH1*EXP(TMP)  
      ELSE  
C  
C        FSTRSELOG=-0.0147351*(VOID**3)+0.333957*(VOID**2)-3.28661*VOID  
C  
        FSTRSELOG=-0.0147351*(VOID**3)+0.311854*(VOID**2)-2.96371*VOID  
     &      +7.34698  
        FSTRSE=EXP(FSTRSELOG)  
      END IF  
      RETURN  
      END  

