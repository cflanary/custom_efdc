      FUNCTION FDSTRSE(VOID,BMECH1,BMECH2,BMECH3)  
      REAL::FDSTRSE,BMECH1,TMP,VOID,BMECH2,BMECH3,FSTRSEL,DFSTRSEL
C CHANGE RECORD  
C  ADDED STANDARD EXPONENTIAL FORM CONSTITUTIVE RELATIONSHIP  
C  
C  
C **  FDSTRSE IS COMPRESSION LENGTH SCALE  
C        STRESS WITH RESPECT TO VOID RATIO  
C  
      IF(BMECH1.GT.0.0)THEN  
        TMP=-(VOID-BMECH2)/BMECH3  
        TMP=-(VOID-BMECH2)/BMECH3  
        FDSTRSE=(BMECH1/BMECH3)*EXP(TMP)  
      ELSE  
C  
C        FDSTRSEL=-0.00835084*(VOID**3)+0.207061*(VOID**2)  
C  
        FSTRSEL=-0.0147351*(VOID**3)+0.311854*(VOID**2)  
     &      -2.96371*VOID+7.34698  
        DFSTRSEL=-0.0442053*(VOID**2)+0.623708*VOID  
     &      -2.96371  
        FDSTRSE=DFSTRSEL*EXP(FSTRSEL)  
      END IF  
      RETURN  
      END  
