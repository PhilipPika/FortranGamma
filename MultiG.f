c
c      SUBROUTINE MULTIG
c
      subroutine MultiG()
!        include 'common_geo.inc'
!        include 'common.inc'
        real*8 kvalue(14,1), fvalue(14,1), POC(2,14), SPOC(14,1), SUMPOC
     +(3,1),POC0,z,a,b,Aa,Bb,w0,db0
        integer nG,xx,i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      open (unit = 1, file = 'f_POC_Values.dat')
              read(1,*) (fvalue(i,1), i = 1,14)
              write(*,*) "Loading f: ",i ," done"
      close(1)

      open (unit = 1, file = 'k_POC_Values.dat')
              read(1,*) (kvalue(i,1), i = 1,14)
              write(*,*) "Loading k: ",i ," done"
      close(1)

      write(*,*) 'This is MultiG being called!'
!     if(xx_pos.gt.xxbiot)then
!     elseif(xx_poxx.le.xxboit)
 
      POC0=13.4
      w0=0.11
      db0=27
      
      do xx = 2,2
      z = (xx-1)/100
     
      SUMPOC=0.0
      do nG = 1,14

        a = (w0-(w0**2+4*db0*kvalue(nG,1))**(1/2))/(2*db0);
c        write(*,*) 'a',a
        b = (w0+(w0**2+4*db0*kvalue(nG,1))**(1/2))/(2*db0);
c        write(*,*) 'b',b
        !b.c. specific
        Aa = -(fvalue(nG,1)*POC0*exp(b*z)*b)/(exp(a*z)*a-exp(b*z)*b);
c        write(*,*) 'Aa',Aa
        Bb = (exp(a*z)*POC0*a*fvalue(nG,1))/(exp(a*z)*a-exp(b*z)*b);
c        write(*,*) 'Bb',Bb

        POC(xx,nG)=Aa*exp(a*z) + Bb*exp(b*z);

        SUMPOC(xx,1) = SUMPOC(xx-1,1) + POC(xx,nG)
      enddo
!        SUMPOC(xx,1) = SUMPOC(xx-1,1) + POC(xx,nG)
      enddo
      write(33,*) SUMPOC
      write(22,*) POC
      END

!      a11=(w0-(w0**2+4*db0(j)*kk)**(1/2))/(2*db0(j));
!      a11_(i)=a11;
!      b11=(w0+(w0**2+4*db0(j)*kk)**(1/2))/(2*db0(j));
!      b11_(i)=b11;
!      
!      a21=(-kk/w0);
!      a21_(i)=a21;
!      
!      A11=-(FF*POCt*b11*exp(b11*zbio))/(a11*exp(a11*zbio)-b11*exp(b11*zbio)); % The same as A above!!
!      A11_(pp,i) = A11;
!      
!      valid = exp(a21*zbio);
!      if valid <= 10**-21
!          valid = 10**-20;
!      end
!      
!      A21=(A11*(exp(a11*zbio)-exp(b11*zbio))+FF*POCt*exp(b11*zbio))/valid;
!      
!      %         A21_(pp,i)=A21;
!      %
!      %         if isnan(A21) == 1 || A21 == Inf || A21 == -Inf
!      %             A21 = 1e-20;
!      %         elseif  A21 >= 0.5 || A21 <= -0.5
!      %             A21 = 1e-20;
!      %         end
!      %         A21__(pp,i)=A21;
!      
!      if(xx(pp)<zbio)
!          POC(pp,i)=A11*(exp(a11*xx(pp))-exp(b11*xx(pp)))+FF*POCt*exp(b11*xx(pp));
!      else
!          
!          %             valevale=exp(a21*xx(pp));
!          %             if valevale == 0
!          %                 valevale = 1e-20;
!          %             end
!          
!          POC(pp,i)=A21*exp(a21*xx(pp));
!      END
