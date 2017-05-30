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
