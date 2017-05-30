C+
C Project     :	Course FORTRAN Codes
C
C Name        :	FIT
C
C Purpose     :	Numerical Recipes linear least squares fitting subroutine.
C
C Explanation :	Given a set of NDATA points X(I), Y(I), with standard 
C               deviations SIG(I), fit them to a straight line y = a + bx by 
C               minimizing chi**2.  Returned are A, B and their respective 
C               probable uncertainties SIGA, SIGB, the chi-squared (CHI2), 
C               and the goodness-of-fit probability Q (that the fit would 
C               have chi**2 this large or larger).  If MWT=0 on input, then 
C               the standard deviations are assumed to be unavailable: 
C               Q is returned as 1.0 and the normalization of CHI2 is to unit 
C               standard deviation on all points.
C
C Use         :	CALL FIT(X, Y, NDATA, SIG, MWT, A, B, SIGA, SIGB, CHI2, Q)
C
C Inputs      :	X:      One-dimensional array of single precision reals of 
C                       size NDATA containing the independent variable data.
C
C		Y:	One-dimensional array of single precision reals of 
C                       size NDATA containing the dependent variable data.  
C                       There must be a one-to-one correspondance between the 
C			array elements in X and Y.
C
C               NDATA:  Integer containing the number of elements in the X 
C                       and Y arrays.
C
C               SIG:    The uncertainty in the Y data.
C
C               MWT:    A switch to let the routine know if standard deviations
C                       (i.e., uncertainties) are known on input for the Y data.
C
C                       MWT = 0:  Standard deviations are unavailable.
C                       MWT > 0:  Standard deviations are available and stored 
C                                 in SIG.
C
C Outputs     :	A:      The y-intercept of the fitted straight line.
C
C               B:      The slope of the fitted straight line.
C
C               SIGA:   The uncertainty in the value of A:  A+/-SIGA.
C
C               SIGB:   The uncertainty in the value of B:  B+/-SIGB.
C
C               CHI2:   The chi-square of the goodness-of-fit.
C
C               Q:      The probability that the goodness-of-fit of the data to  
C                       a straight line is believable.  Note that if 
C
C                       Q > 0.1:    The goodness-of-fit is believable.
C                       0.001 < Q < 0.1:  The fit may be acceptable if the errors 
C                                   are nonnormal or have been moderately 
C                                   underestimated.
C                       Q < 0.001:  The goodness-of-fit is questionable.
C
C Calls       :	Function GAMMQ.
C
C Common      :	None.
C
C Restrictions:	FORTRAN 77 coding.
C
C Side effects:	None.
C
C Category    :	Data fitting.
C
C Prev. Hist. :	Based on the FIT subroutine of Numerical Recipes for FORTRAN 77.
C
C Written     :	Donald G. Luttermoser, ETSU/Physics, 2 October 2013.
C
C Modified    :	Version 1, Donald G. Luttermoser, ETSU/Physics, 2 Oct 2013
C			Initial program.
C
C Version     :	Version 1,  2 October 2013.
C
C-
      SUBROUTINE FIT(X, Y, NDATA, SIG, MWT, A, B, SIGA, SIGB, CHI2, Q)
C
      INTEGER MWT, NDATA
      REAL A, B, CHI2, Q, SIGA, SIGB, SIG(NDATA), X(NDATA), Y(NDATA)
C
      INTEGER I
      REAL SIGDAT, SS, ST2, SX, SXOSS, SY, T, WT, GAMMQ
C
C Initialize sums to zero.
C
      SX = 0.
      SY = 0.
      ST2 = 0.
      B = 0.
C
C Accumulate sums
C
      IF (MWT .NE. 0) THEN
          SS = 0.
C
C  with weights:
C
          DO 11 I = 1, NDATA
              WT = 1. / (SIG(I)**2)
              SS = SS + WT
              SX = SX + X(I)*WT
              SY = SY + Y(I)*WT
 11       CONTINUE
      ELSE
C
C  or without weights:
C
          DO 12 I = 1, NDATA
              SX = SX + X(I)
              SY = SY + Y(I)
 12       CONTINUE
          SS = FLOAT(NDATA)
      ENDIF
C
      SXOSS = SX / SS
C
      IF (MWT .NE. 0) THEN
          DO 13 I = 1, NDATA
              T = (X(I) - SXOSS) / SIG(I)
              ST2 = ST2 + T*T
              B = B + T*Y(I) / SIG(i)
13        CONTINUE
      ELSE
          DO 14 I = 1, NDATA
              T = X(I) - SXOSS
              ST2 = ST2 + T*T
              B = B + T*Y(I)
 14       CONTINUE
      ENDIF
C
C Solve for A, B, SIGA, and SIGB.
C
      B = B / ST2
      A = (SY - SX*B) / SS
      SIGA = SQRT((1. + SX*SX / (SS*ST2)) / SS)
      SIGB = SQRT(1. /ST2)
C
C Calculate chi-square.
C
      CHI2 = 0.
      Q = 1.
      IF (MWT .EQ. 0) THEN
          DO 15 I = 1, NDATA
              CHI2 = CHI2 + (Y(I) - A - B*X(I))**2
 15       CONTINUE
          SIGDAT = SQRT(CHI2 / (FLOAT(NDATA-2)))
          SIGA = SIGA * SIGDAT
          SIGB = SIGB * SIGDAT
      ELSE
          DO 16 I = 1, NDATA
              CHI2 = CHI2 + ((Y(I) - A - B*X(I)) / SIG(I))**2
 16       CONTINUE
          IF (NDATA .GT. 2) Q = GAMMQ(0.5*(FLOAT(NDATA-2)), 0.5*CHI2)
      ENDIF
C
      RETURN
      END
