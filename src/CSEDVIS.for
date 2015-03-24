      REAL FUNCTION CSEDVIS(SED)  
      REAL::SED,VISR,WTL,WTH,VISL,VISH
C CHANGE RECORD  
C  
C  
C **  CALCULATES KINEMATIC VISCOSITY OF HIGH CONCENTRATION COHESIVE  
C **  SEDIMENT-WATER MIXTURE BASED ON  
C **  
C **  MEHTA, A. J., AND F.JIANG, 1990: SOME OBSERVATIONS ON BOTTOM  
C **  MUD MOTION DUE TO WAVES. COASTAL AND OCEANOGRAPHIC ENGINEERING  
C **  DEPARTMENT, UNIVERSITY OF FLORIDA, GAINESVILLE, FL32661  
C  
      IF(SED.LE.25667.) VISR=0.116883E-3*SED  
      IF(SED.GE.36667.) VISR=1.52646E-6*SED+3.125  
      IF(SED.GT.25667.0.AND.SED.LT.36667.0)THEN  
        WTL=(36667.-SED)/11000.  
        WTH=(SED-25667.)/11000.  
        VISL=0.116883E-3*25667  
        VISH=1.52646E-6*36667.+3.125  
        VISR=WTL*VISL+WTH*VISH  
      ENDIF  
      CSEDVIS=1.E-6*(10.**VISR)  
      RETURN  
      END  

