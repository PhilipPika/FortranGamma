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
