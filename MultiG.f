c
c      SUBROUTINE MULTIG
c
      subroutine MultiG()
c        include 'common_geo.inc'
c        include 'common.inc'
      real*8 kvalue(14,1),fvalue(14,1), POC(1000,14), SUMPOC(1000,1),
     +POC0,a,b,Aa,Bb,w0,db0,z,zz,SPOC(1000,1)
        integer nG,xx,i
c_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
c_/                        Loading in data                           _/        
c_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
      write(*,*) 'This is MultiG being called!'
      open (unit=1, file='f_POC_Values.dat')
              read(1,*) (fvalue(i,1), i=1,14)
              write(*,*) "Loading f: ",i ," done"
      close(1)
      open (unit=1, file='k_POC_Values.dat')
              read(1,*) (kvalue(i,1), i=1,14)
              write(*,*) "Loading k: ",i ," done"
      close(1)
c_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
c_/                        Loading in parameter                      _/        
c_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

      POC0=13.3
      w0=0.11
      db0=27

c_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
c_/                        Start depth for-loop                      _/        
c_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
      
      do xx=1,1000
      zz=(xx-1)
      z=zz/100
c_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
c_/                        Start fractions for-loop                  _/        
c_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/

      SUMPOC=0.0
      do nG=1,14

!        a=(w0-(w0**2+4*db0*kvalue(nG,1))**(0.5))/(2*db0)!Works OK
!        a=(w0-(w0**2+4*db0*kvalue(nG,1))**(1/2))/(2*db0)!Doesn't work
        a=(w0-sqrt(w0**2+4*db0*kvalue(nG,1)))/(2*db0)
!        b=(w0+(w0**2+4*db0*kvalue(nG,1))**(0.5))/(2*db0)!Works OK
!        b=(w0+(w0**2+4*db0*kvalue(nG,1))**(1/2))/(2*db0)!Doesn't work
        b=(w0+sqrt(w0**2+4*db0*kvalue(nG,1)))/(2*db0)
C Specific to boundary conditions        
        Aa=-(fvalue(nG,1)*POC0*exp(b*z)*b)/(exp(a*z)*a-exp(b*z)*b)
        Bb=(exp(a*z)*POC0*a*fvalue(nG,1))/(exp(a*z)*a-exp(b*z)*b)

        POC(xx,nG)=Aa*exp(a*z)+Bb*exp(b*z);

C       Debugging part
!        write(*,*) 'z',z,'nG',nG,'a',a,'b',b,'Aa',Aa,'Bb',Bb,'POC',POC(x
!     +x,nG)
!      write(*,*) 'z',z,'nG',nG, 'POC', POC(xx,nG),'SPOC',SUMPOC(xx,1)
!      write(*,*) 'k', nG,'=',kvalue(nG,1),'F',nG,'=',fvalue(nG,1)
      

       SUMPOC(xx,1)=SUMPOC(xx,1)+POC(xx,nG)! Works now as as those below
c       work too
      enddo
      write(33,*) SUMPOC(xx,1)
      
!      SUMPOC(xx,1)=POC(xx,1)+POC(xx,2)+POC(xx,3)+POC(xx,4)+POC(xx,5)+POC
!     +(xx,6)+POC(xx,7)+POC(xx,8)+POC(xx,9)+POC(xx,10)+POC(xx,11)+POC(
!     +xx,12)+POC(xx,13)+POC(xx,14)
!      write(33,*) SUMPOC(xx,1)
!      SPOC=0.0
!      DO i=1,14
!      SPOC(xx,1)=SPOC(xx,1)+POC(xx,i)
!      ENDDO
!      write(44,*) SPOC(xx,1)
     
      enddo
      
      write(22,*) POC

      END
