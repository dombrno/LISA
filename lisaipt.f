C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: lisaipt.f
C NOTICE : This program accompanies the revue article:
C
C           The Local Impurity Self Consistent Approximation (LISA)
C                   to Strongly Correlated Fermion Systems
C                   and the Limit of Infinite Dimensions
C
C                                   by
C
C          A. Georges, G. Kotliar, W. Krauth, M. Rozenberg
C
C            to be published in: Reviews of Modern Physics (1996)
C
C          (the paper will be referred to as ``GKKR'').
C          you are kindly asked to cite the paper 
C                   if you use this program.
C
C          Unfortunately, we cannot guarantee correctness of the
C          programs,
C          but they have been thoroughly tested on SUN Sparks, HP 9000,
C          IBM RS6000 stations.
C
C TYPE   : main
C PURPOSE: Iterated perturbation theory at finite temperature.
C          An essentially identical version of the program was used
C          in: A. Georges W. Krauth, Phys. Rev. B 48 (1993). 
C          In this version of the program only the paramagnetic 
C          half-filled case is treated
C
C I/O    : 
C VERSION: 30-Nov-95
C AUTHOR : W. Krauth (krauth@physique.ens.fr)
C COMMENT: Even though FORTRAN is case-insensitive, I have
C          exclusively capitalized all
C          the global variables, i. e. the variables appearing in the
C          COMMON block /global/ (cf. file lisaipt.dat)
C
C          Convention used for imaginary-time Green's functions:
C                      G(tau)=-T<    > 'Landau-Lifshitz'
C
C          (this is the opposite of what is done in lisaqmc.f )         
C
C          How to choose the parameters L and iwmax: Do
C          your own extrapolation. In the above paper, L =2^9, iwmax=2^13.
C          this gave convergence down to temperatures of T= 0.05. 
C          Remember that omega_max should be much larger than the 
C          bandwidth. In any case, to go down to yet lower (but finite)
C          temperatures, you may have to treat the asymptotic tail differently.
C          It would also be useless to do an explicit calclulation for
C          every single Matsubara frequency
C          to run the program, compile the Fortran file, and execute
C          the 2 commands:
C 
C          rm -f lisaipt.G lisaipt.Greent lisaipt.Greenw
C          lisaipt.out
C
C========+=========+=========+=========+=========+=========+=========+=$
      include 'lisaipt.dat'
      complex*16 xi,argu,omega,wfun
      external dens
      Zero=0
      One=1
      Two=2
      Pi=acos(-One)
      xi=cmplx(Zero,One)
      ic=0
cc
cc    (ic 0:Gauss; 1:square; 2:semicirc; 3:half circle;
cc
      open(1,file='lisaipt.input',status='old')
      read(1,*)U,Beta
      open (unit=3,file='lisaipt.Greent',form='formatted',status='new')
      open (unit=4,file='lisaipt.Greenw',form='formatted',status='new')
cc
cc    Here we initialize Green0t with a 'fantasy'- function.
cc    The interested user may want to replace this function with
cc    a more educated guess.
cc
      do  i=1,L
          Green0t(i)=-.5 +i/real(L)/5.
      end do
      call fourier(Green0t,Green0w)

cc
cc    do any loop you want...here we produce a table of energy
cc    versus temperature at fixed U. Derivation of the curve
cc    gives the specific heat.
cc
      delt=-.01
      do Temp=.4,.02,delt
         Beta=1./temp
	 print*,U,1./Beta,'  =============================U temp'
         do idiff=1,80
c           Define sigmat and FT
cc
cc          The fourier transform program fourier is for functions with a jump 
cc          condition of 1 (Green's functions). Now, we have to 
cc          transform Sigmat, which has a jump which will be
cc          called  xjump (xjump=.25 at half-filling). 
cc          We scale Sigmat, Fourier transform it, 
cc          and scale back immediately afterwards. 
cc      
cc
            xjump=.25
            do i=1,L
               Sigmat(i) = Green0t(i)**3/xjump
            end do 
            call fourier(Sigmat,Sigmaw)
cc
cc          The following 4 lines code the IPT approximation
cc          and do the rescaling
cc
            do i=0,Iwmax
               Sigmaw(i)=-conjg(Sigmaw(i))*U**2*xjump
               omega=(Two*i+One)*Pi/Beta*xi
            end do
cc
cc  Calculate Greenw and update Green0w
cc
            energy=Zero
            do i=0,Iwmax
               omega=(Two*i+One)*Pi/Beta*xi
               argu=omega-Sigmaw(i)
               f=One
               if (imag(argu).le.0) f=-One
cc
cc             Gaussian
cc
               if (ic.eq.0) then
   	          argu=argu*f
                  Greenw(i)=-f*xi*sqrt(Pi)*wfun(argu)
cc
cc                Square
cc
               else if (ic.eq.1) then
	             argu=(argu+One)/(argu-One)
	             Greenw(i)=log(argu)/Two
cc
cc                semi circular
cc
               else if (ic.eq.2) then
	             do iii=1,8
	                Greenw(i)=One/(-Greenw(i)/Two+argu)
	             end do
               end if
cc
cc             use one of the formulas below...the second one is 
cc             designed to break 2-cycles
cc

 	       Green0w(i)=One/(One/Greenw(i)+Sigmaw(i))
c              Green0w(i)=(Green0w(i)+One/(One/Greenw(i)+Sigmaw(i)))/Two
cc
cc             energy calculation
cc
               energy=energy - 4.+4.*Greenw(i)*(omega-Sigmaw(i)) 
     &         + Two*Sigmaw(i)*Greenw(i)
            end do 
            energy=energy/Beta+U/4.
            call invfourier(Green0w,Green0t)
cc 
cc  Symmetrisation of Green0t: We only treat the half-filled case here
cc
	    Green0t(1)=-One/Two
	    do i=2,L/2	
	       dummy=Green0t(i)+Green0t(L+2-i)
               Green0t(i)=dummy/Two
               Green0t(L+2-i)=dummy/Two
            end do
            do i=1,L
               Green0t(i)=-abs(Green0t(i))
            end do
c
cc  inverse Fourier transform to get Greent and symmetrize
            call invfourier(Greenw,Greent)
            write(*,'(2f20.10)')One/Beta,-Greent(L/2)
            if (abs(oldgreen-Greent(L/2)).lt.1.e-6) then
               goto 913
            else
               oldgreen=Greent(L/2)
            end if
	 end do
         print*, ' ATTENTION : NO CONVERGENCE OBTAINED'
cc
cc       notice that we have made no effort to obtain maximum
cc       precision....sometimes, cyles appear... they can be
cc       treated with standard procedures. 
cc
913      continue
cc
cc       plot Greenw for the first Matsubara frequencies...Adapt as 
cc       you like
cc
	 write(4,*)'" U Beta',U,Beta,'"'
  	 do i=0,min(100,iwmax)
            om=(Two*i+One)*Pi/Beta
  	    write(4,*)om,imag(Greenw(i))
         end do
         write(4,*)
	 write(45,*)'# U Beta',U,Beta
c        do eps=-1.95,1.95,.1
c           write(45,*)eps,dens(eps)
c        end do
         write(3,*)'#'
         do i=1,L
c           write(3,*)(i-1)/real(L)*Beta,-Greent(i)
            write(3,*)i,-Greent(i)
         end do
	 rewind(2)
         call wheader(2)
	 do i=1,L
 	    write(2,*)(i-1)/real(L)*Beta, Green0t(i)
         end do
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: dens
C TYPE   : function
C PURPOSE: calculate N(epsilon) from Green's function
C I/O    :
C VERSION: 6-16-92
C COMMENT: calculate difference with free fermi gas
C========+=========+=========+=========+=========+=========+=========+=$
      function dens(eps)
      include 'lisaipt.dat'
      complex*16 sigma, xi,omega
      xi=Cmplx(Zero,One)
      dens=0.
      do i=0,Iwmax
         omega=(Two*i+One)*Pi/Beta*xi
         sigma=One/Green0w(i)-One/Greenw(i) 
         dens=dens+Two/(omega-sigma-eps)-Two/(omega-eps)
      end do
      dens=dens/beta+One/(exp(Beta*eps)+One)
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: fourier
C TYPE   : subroutine
C PURPOSE: fourier-transform the function Green(tau) in order to
C          calculate function Green(omega)
C I/O    :
C VERSION: 4-16-92
C COMMENT: It may be noticed that in this program,
C          a different interpolation than in lisaself.f is used
C          (here: piecewise linear interpolation of Greent))
C          (in lisaself: spline interpolation of Greent)
C          In lisaself.f, we are limited in the number of
C          timeslices that we can use: L is fixed, increasing L
C          is difficult, and we should use the best interpolation
C          possible. This is the spline intepolation used there.
C          Here, increasing L is very simple, so we use a
C          piece-wise linear interpolation, in combination with
C          a fast Fourier transform. 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine fourier(rindata,coutdata)
      include 'lisaipt.dat'
      dimension rindata(L)
      complex*16 data(2*L),coutdata(0:Iwmax),xi,cdummy
      xi=(0.,1.)
      do  i=1,L
         data(L+i)=rindata(i)
         data(i)=-rindata(i)
      end do
      call fft(data,2*L,1)
      do  i=0,L/2-1
         coutdata(i)=-data(2*i+2)*beta/2./real(L) 
	 coutdata(L-i-1)=conjg(coutdata(i))
      end do
cc
cc    extend
cc
      do i=L,Iwmax
           coutdata(i)=coutdata(i-L)
      end do
cc
cc    put in attenuation factors: cf Stoer and Bulirsch p. 88 ff
cc
cc    if you want to understand the following formula without
cc    going through Stoer and Bulirsch, simply write down the 
cc    analytic formula for the Fourier transform of the 
cc    piecewise linear interpolation of G(tau)  with the 
cc    correct jump condition at tau=0. One can also try 
cc    spline interpolation, as mentioned in the header of this
cc    subroutine
cc
      do i=0,Iwmax
         om=(2.*i+1.)*Pi/Beta
         delta=Beta/L
         cdummy=xi*om*delta
         coutdata(i)=1./xi/om -(exp(-cdummy)-1.)/om**2/delta
     1   + 4.*sin(om*delta/2.)**2/om**2/delta**2*coutdata(i)
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: invfourier
C TYPE   : subroutine
C PURPOSE: inverse fourier transform
C          Greent, Greenw use physical definition
C          Greent(i)=G((i-1)*deltau) for i=1,...,T
C          Greenw(n)=G(i w_n), for n=0,T/2-1
C                 w_n=(2*n+1)pi/beta
C          Symmetry property: 
C          G(iw_(-n)=G(iw_(n-1))*
C          coupled to the impurity
C I/O    :
C VERSION: 11-18-91
C COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine invfourier(cindata,routdata)
      include 'lisaipt.dat'
      dimension routdata(L)
      complex*16 cindata(0:Iwmax),cinprime(2*Iwmax)
      do j=0,Iwmax-1
         cinprime(2*j+1)=0.
	 cinprime(2*j+2)=cindata(j)
      enddo
      call fft(cinprime,2*Iwmax,-1)
      do i=1,L
         chu=cinprime((i-1)*Iwmax/L+1)*2/beta
         routdata(i)=chu
      end do
cc
cc    special treatment for tau=0
cc
      routdata(1)=-1./2.+routdata(1)
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: fft.f
C TYPE   : subroutine
C PURPOSE: discrete fast fourier transform
C I/O    :
C VERSION: 11-13-91
C COMMENT: 
C          standard storage of frequencies: first pos frequ, 
C          then neg frequ.
C          NN must be power of 2!
C          isign= 1:     Fourier transform
C          isign=-1: inverse Fourier transform
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine fft(data,nn,isign)
      implicit double precision (a-h,o-z)
      dimension data(2*nn)
      n=2*nn
      j=1
      do i=1,n,2
         if (j.gt.i) then
            tempr=data(j)
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
         end if
         m=n/2
1        if ((m.ge.2).and.(j.gt.m)) then
            j=j-m
            m=m/2
            go to 1
         end if
         j=j+m
      end do
      mmax=2
2     if (n.gt.mmax) then
         istep=2*mmax 
         theta=6.28318530717959d0/(isign*mmax)
         wpr=-2.d0*sin(0.5d0*theta)**2
         wpi=dsin(theta)
         wr=1.d0
         wi=0.d0
         do m=1,mmax,2
            do i=m,n,istep
               j=i+mmax
               tempr=(wr)*data(j)-(wi)*data(j+1)
               tempi=(wr)*data(j+1)+(wi)*data(j)
               data(j)=data(i)-tempr
               data(j+1)=data(i+1)-tempi
               data(i)=data(i)+tempr
               data(i+1)=data(i+1)+tempi
            end do
            wtemp=wr
            wr=wr*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
         end do
         mmax=istep
         go to 2
      endif
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: rheader
C TYPE   : subroutine
C PURPOSE: read the program's header from unit k 
C I/O    :
C VERSION: 11-18-91
C COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine rheader(k) 
      do 1 i=1,5
         read(k,*)
1     continue
      end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: wheader.f  
C TYPE   : subroutine
C PURPOSE: write the program's header onto standard output 
C I/O    : 
C VERSION: 11-18-91
C COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine wheader(k) 
      include 'lisaipt.dat'
      write(k,'(a55)')'========================================'
      write(k,'(a55)')' LISAIPT                                '
      write(k,'(a55)')'========================================'
      write(k,'(4(a6,I4))') 'L=',L
      write(k,'(2(a6,f15.7),a6,I6)')'Beta=',Beta,'U= ',U,
     &	'Iwmax',Iwmax
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: wfun.f
C TYPE   : subroutine
C PURPOSE: calculate exp rescaled error function w (Abr-Stegun)
C I/O    :
C VERSION: 11-28-91
C COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
      function wfun(z)
      complex*16 z,wfun
      double precision x,y,u,v
      logical flag
      x=dreal(z)
      y=imag(z)
      call wofz(x,y,u,v,flag)
      wfun=cmplx(u,v)
      end  
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: wofz
C TYPE   : subroutine
C PURPOSE: calculate exp rescaled error function w (Abr-Stegun)
C I/O    :
C VERSION: 11-28-91
C COMMENT: the following program was obtained from NETLIB
C========+=========+=========+=========+=========+=========+=========+=$
C      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 16, NO. 1, PP. 47.
      SUBROUTINE WOFZ (XI, YI, U, V, FLAG)
*
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
*
      LOGICAL A, B, FLAG
      PARAMETER (FACTOR   = 1.12837916709551257388D0,
     *           RMAXREAL = 1.340780792994259D+154,
     *           RMAXEXP  = 709.0895657128241D0,
     *           RMAXGONI = 0.6746518850690209D10)
*
      FLAG = .FALSE.
*
      XABS = DABS(XI)
      YABS = DABS(YI)
      X    = XABS/6.3
      Y    = YABS/4.4
*
C
C     THE FOLLOWING IF-STATEMENT PROTECTS
C     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
C
      IF ((XABS.GT.RMAXREAL).OR.(YABS.GT.RMAXREAL)) GOTO 100
*
      QRHO = X**2 + Y**2
*
      XABSQ = XABS**2
      XQUAD = XABSQ - YABS**2
      YQUAD = 2*XABS*YABS
*
      A     = QRHO.LT.0.085264D0
*
      IF (A) THEN
C
C  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
C  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
C  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
        QRHO  = (1-0.85*Y)*DSQRT(QRHO)
        N     = IDNINT(6 + 72*QRHO)
        J     = 2*N+1
        XSUM  = 1.0/J
        YSUM  = 0.0D0
        DO 10 I=N, 1, -1
          J    = J - 2
          XAUX = (XSUM*XQUAD - YSUM*YQUAD)/I
          YSUM = (XSUM*YQUAD + YSUM*XQUAD)/I
          XSUM = XAUX + 1.0/J
 10     CONTINUE
        U1   = -FACTOR*(XSUM*YABS + YSUM*XABS) + 1.0
        V1   =  FACTOR*(XSUM*XABS - YSUM*YABS)
        DAUX =  DEXP(-XQUAD)
        U2   =  DAUX*DCOS(YQUAD)
        V2   = -DAUX*DSIN(YQUAD)
*
        U    = U1*U2 - V1*V2
        V    = U1*V2 + V1*U2
*
      ELSE
C
C  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
C  CONTINUED FRACTION
C  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
C  ACCURACY
C
C  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
C  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
C  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
C  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
C  TO OBTAIN THE REQUIRED ACCURACY
C  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
C  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
C
*
        IF (QRHO.GT.1.0) THEN
          H    = 0.0D0
          KAPN = 0
          QRHO = DSQRT(QRHO)
          NU   = IDINT(3 + (1442/(26*QRHO+77)))
        ELSE
          QRHO = (1-Y)*DSQRT(1-QRHO)
          H    = 1.88*QRHO
          H2   = 2*H
          KAPN = IDNINT(7  + 34*QRHO)
          NU   = IDNINT(16 + 26*QRHO)
        ENDIF
*
        B = (H.GT.0.0)
*
        IF (B) QLAMBDA = H2**KAPN
*
        RX = 0.0
        RY = 0.0
        SX = 0.0
        SY = 0.0
*
        DO 11 N=NU, 0, -1
          NP1 = N + 1
          TX  = YABS + H + NP1*RX
          TY  = XABS - NP1*RY
          C   = 0.5/(TX**2 + TY**2)
          RX  = C*TX
          RY  = C*TY
          IF ((B).AND.(N.LE.KAPN)) THEN
            TX = QLAMBDA + SX
            SX = RX*TX - RY*SY
            SY = RY*TX + RX*SY
            QLAMBDA = QLAMBDA/H2
          ENDIF
 11     CONTINUE
*
        IF (H.EQ.0.0) THEN
          U = FACTOR*RX
          V = FACTOR*RY
        ELSE
          U = FACTOR*SX
          V = FACTOR*SY
        END IF
*
        IF (YABS.EQ.0.0) U = DEXP(-XABS**2)
*
      END IF
*
*
C
C  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
C
*
      IF (YI.LT.0.0) THEN
*
        IF (A) THEN
          U2    = 2*U2
          V2    = 2*V2
        ELSE
          XQUAD =  -XQUAD
*
C
C         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
C         AGAINST OVERFLOW
C
          IF ((YQUAD.GT.RMAXGONI).OR.
     *        (XQUAD.GT.RMAXEXP)) GOTO 100
*
          W1 =  2*DEXP(XQUAD)
          U2  =  W1*DCOS(YQUAD)
          V2  = -W1*DSIN(YQUAD)
        END IF
*
        U = U2 - U
        V = V2 - V
        IF (XI.GT.0.0) V = -V
      ELSE
        IF (XI.LT.0.0) V = -V
      END IF
*
      RETURN
*
  100 FLAG = .TRUE.
      RETURN
*
      END
