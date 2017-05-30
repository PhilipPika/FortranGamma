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
