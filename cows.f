      PROGRAM cows
C*********************************************************************C
C       Program to calcualte the Gamma function for the MultiG        C
C       organic matter degradation in the bioturbated layer in        C
C       the BRNS model.                                               C
C       May 2017, PP                                                 C
C*********************************************************************C
         IMPLICIT NONE
         REAL :: gammln,lll,ll,gammp,pp,GAMMQ
         PRINT *,"gamma ln"
         PRINT *,gammln(0.125)
         ll = gammln(0.125)
         lll = exp(ll)
         PRINT *,lll
         PRINT *,"gammaq"
         PRINT *,GAMMQ(0.125,0.0003)
         PRINT *,"gammap"
         PRINT *,gammp(0.125,0.0003)

         PRINT *,"Here as 1-GAMMQ"
         pp = 1 - GAMMQ(0.125,0.0003)
         PRINT *,pp
         call MultiG()

      END PROGRAM
C
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
C
C
        FUNCTION gammp(a,x)
        REAL a,gammp,x
C       USES gcf,gser
C       Returns the incomplete gamma function P(a,x).
        REAL gammcf,gamser,gln
C       if(x.lt.0..or.a.le.0.)!pause ’bad arguments in gammp’ 
        if(x.lt.a+1.)then !Use the series representation.
                call gser(gamser,a,x,gln)
                gammp=gamser
        else !Use the continued fraction representation
                call gcf(gammcf,a,x,gln)
                gammp=1.-gammcf ! and take its complement
        endif
        return 
        END
C
C
C+
C Project     :	Course FORTRAN Codes
C
C Name        :	GAMMQ
C
C Purpose     :	Numerical Recipes incomplete gamma function Q.
C
C Explanation :	This function actually returns the complement of the 
C               the incomplete gamma function, 1 - P(a,x), where P(a,x) is 
C               actually the incomplete gamma function defined by
C
C                   P(a,x) = (1./GAMMA(a)) INT_0^x EXP(-t) t**(a-1) dt,
C
C               where GAMMA(a) is the gamma function defined by
C
C                   GAMMA(a) = INT_0^INFTY t**(a-1) EXP(-t) dt.
C
C Use         :	Result = GAMMQ(A, X)
C
C Inputs      :	A:      Scalar containing the value of the exponential of 
C                       the gamma function.
C
C		X:	Scalar containing the upper limit of the integral of 
C                       the incomplete gamma function.
C
C Outputs     :	Result: Returns the value of the compliment of the incomplete
C                       gamma function.
C
C Calls       :	Subroutines GCF, GSER.
C
C Common      :	None.
C
C Restrictions:	FORTRAN 77 coding.
C
C Side effects:	None.
C
C Category    :	Data fitting.
C
C Prev. Hist. :	Based on the GAMMQ function of Numerical Recipes for FORTRAN 77.
C
C Written     :	Donald G. Luttermoser, ETSU/Physics, 2 October 2013.
C
C Modified    :	Version 1, Donald G. Luttermoser, ETSU/Physics, 2 Oct 2013
C			Initial program.
C
C Version     :	Version 1,  2 October 2013.
C
C-
      FUNCTION GAMMQ(A, X)
C
      REAL A, GAMMQ, X
C
C Define internal parameters.
C
      REAL GAMMCF, GAMSER, GLN
C
C Check passed parameters.
C
      IF ((X .LT. 0.) .OR. (A .LE. 0.)) THEN
          PRINT *, 'Bad arguments in GAMMQ.'
          GAMMQ = 0.
          RETURN
      ENDIF
C
      IF (X .LT. A+1.) THEN
          CALL GSER(GAMSER, A, X, GLN)
          GAMMQ = 1. - GAMSER
      ELSE
          CALL GCF(GAMMCF, A, X, GLN)
          GAMMQ = GAMMCF
      ENDIF
C
      RETURN
      END
C
C
C+
C Project     :	Course FORTRAN Codes
C
C Name        :	GCF
C
C Purpose     :	Numerical Recipes incomplete gamma function continued fraction.
C
C Explanation :	This subroutine calculates the continued fraction solution for 
C               the incomplete gamma function given by Eq. (6.2.6) in Numerical
C               Recipes-- The Art of Scientific Computing (1986) by Press, et al.
C
C Use         :	CALL GCF(GAMMCF, A, X, GLN)
C
C Inputs      :	A:      Scalar containing the value of the exponential of 
C                       the gamma function.
C
C		X:	Scalar containing the upper limit of the integral of 
C                       the incomplete gamma function.
C
C Outputs     :	GAMMCF: Returns the value of the continued fraction solution 
C                       for the incomplete gamma function.
C
C               GLN:    Natural log of the gamma function.
C
C Calls       :	Function GAMMLN.
C
C Common      :	None.
C
C Restrictions:	FORTRAN 77 coding.
C
C Side effects:	None.
C
C Category    :	Data fitting.
C
C Prev. Hist. :	Based on the GSER subroutine of Numerical Recipes for FORTRAN 77.
C
C Written     :	Donald G. Luttermoser, ETSU/Physics, 2 October 2013.
C
C Modified    :	Version 1, Donald G. Luttermoser, ETSU/Physics, 2 Oct 2013
C			Initial program.
C
C Version     :	Version 1,  2 October 2013.
C
C-
      SUBROUTINE GCF(GAMMCF, A, X, GLN)
C
      INTEGER ITMAX
      REAL A, GAMMCF, GLN, X, EPS, FPMIN
      PARAMETER (ITMAX=100, EPS=3.E-7, FPMIN=1.E-30)
C
      INTEGER I
      REAL AN, B, C, D, DEL, H, GAMMLN
C
      GLN = GAMMLN(A)
      B = X + 1. - A
      C = 1. / FPMIN
      D = 1. / B
      H = D
      DO 11 I = 1, ITMAX
          AN = -I * (I - A)
          B = B + 2.
          D = AN*D + B
          IF (ABS(D) .LT. FPMIN) D = FPMIN
          C = B + AN/C
          IF (ABS(C) .LT. FPMIN) C = FPMIN
          D = 1. / D
          DEL = D * C
          H = H * DEL
          IF (ABS(DEL-1.) .LT. EPS) GOTO 1
 11   CONTINUE
C
      PRINT *, 'A too large, ITMAX too small in GCF.'
  1   GAMMCF = EXP(-X + A*LOG(X) - GLN) * H
C
      RETURN
      END
C
C
C+
C Project     :	Course FORTRAN Codes
C
C Name        :	GSER
C
C Purpose     :	Numerical Recipes incomplete gamma function series solution.
C
C Explanation :	This subroutine calculates the series solution for the 
C               incomplete gamma function:
C
C                   GAMSER = EXP(-x) x**a * 
C                       SUM_(n=0)^INFTY (GAMMA(A)/GAMMA(a+1+n)) x**n
C
C               where GAMMA(a) is the gamma function defined by
C
C                   GAMMA(a) = INT_0^INFTY t**(a-1) EXP(-t) dt.
C
C Use         :	CALL GSER(GAMSER, A, X, GLN)
C
C Inputs      :	A:      Scalar containing the value of the exponential of 
C                       the gamma function.
C
C		X:	Scalar containing the upper limit of the integral of 
C                       the incomplete gamma function.
C
C Outputs     :	GAMSER: Returns the value of the series defined above.
C
C               GLN:    Natural log of the gamma function.
C
C Calls       :	Function GAMMLN.
C
C Common      :	None.
C
C Restrictions:	FORTRAN 77 coding.
C
C Side effects:	None.
C
C Category    :	Data fitting.
C
C Prev. Hist. :	Based on the GSER subroutine of Numerical Recipes for FORTRAN 77.
C
C Written     :	Donald G. Luttermoser, ETSU/Physics, 2 October 2013.
C
C Modified    :	Version 1, Donald G. Luttermoser, ETSU/Physics, 2 Oct 2013
C			Initial program.
C
C Version     :	Version 1,  2 October 2013.
C
C-
      SUBROUTINE GSER(GAMSER, A, X, GLN)
C
      INTEGER ITMAX
      REAL A, GAMSER, GLN, X, EPS
C
      PARAMETER (ITMAX=100, EPS=3.E-7)
C
      INTEGER N
      REAL AP, DEL, SUM, GAMMLN
C
C Call the GAMMLN function.
C
      GLN = GAMMLN(A)
      IF (X .LE. 0.) THEN
          IF (X .LT. 0.) PRINT *, 'X < 0 in GSER.'
          GAMSER = 0.
          RETURN
      ENDIF
C
      AP = A
      SUM = 1. / A
      DEL = SUM
C
      DO 11 N = 1, ITMAX
          AP = AP + 1.
          DEL = DEL * X / AP
          SUM = SUM + DEL
          IF (ABS(DEL) .LT. ABS(SUM)*EPS) GOTO 1
 11   CONTINUE
C
      PRINT *, 'A too large, ITMAX too small in GSER.'
C
  1   GAMSER = SUM * EXP(-X + A*LOG(X) - GLN)
C
      RETURN
      END
