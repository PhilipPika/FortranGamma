C+
C Project     :	Course FORTRAN Codes
C
C Name        :	GAMMLN
C
C Purpose     :	Numerical Recipes natural log of the gamma function.
C
C Explanation :	This function calculates the natural log of the gamma 
C               function for XX > 0, where the gamma function, GAMMA(a), 
C               is given by:
C
C                   GAMMA(XX) = INT_0^INFTY t**(a-1) EXP(-t) dt.
C
C Use         :	Result = GAMMLN(XX)
C
C Inputs      : XX:     Scalar containing the value passed to the gamma function.
C
C Outputs     :	Result: ln(GAMMA(XX)).
C
C Calls       :	None.
C
C Common      :	None.
C
C Restrictions:	FORTRAN 77 coding.
C
C Side effects:	None.
C
C Category    :	Data fitting.
C
C Prev. Hist. :	Based on the GAMMLN function of Numerical Recipes for FORTRAN 77.
C
C Written     :	Donald G. Luttermoser, ETSU/Physics, 2 October 2013.
C
C Modified    :	Version 1, Donald G. Luttermoser, ETSU/Physics, 2 Oct 2013
C			Initial program.
C
C Version     :	Version 1,  2 October 2013.
C
C-
      FUNCTION GAMMLN(XX)
      REAL GAMMLN, XX
      INTEGER J
      DOUBLE PRECISION SER, STP, TMP, X, Y, COF(6)
      SAVE COF, STP
C
      DATA COF, STP/76.18009172947146D0, -86.50532032941677D0,
     +    24.01409824083091D0, -1.231739572450155D0, 
     +    0.1208650973866179D-2, -0.5395239384953D-5, 
     +    2.5066282746310005D0/
C
      X = XX
      Y = X
      TMP = X + 5.5D0
      TMP = (X + 0.5D0) * LOG(TMP) - TMP
      SER = 1.000000000190015D0
C
      DO 11 J = 1, 6
          Y = Y + 1.D0
          SER = SER + COF(J)/Y
 11   CONTINUE
C
      GAMMLN = TMP + LOG(STP * SER / X)
C
      RETURN
      END
