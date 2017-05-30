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
