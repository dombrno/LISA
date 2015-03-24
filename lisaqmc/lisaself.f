C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: lisaself.f
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
C          you are kindly asked to cite the paper (and, if applicable,
C              the original works) if you use this program.
C
C          Unfortunately, we cannot guarantee correctness of the
C          programs,
C          but they have been thoroughly tested on SUN Sparks, HP 9000,
C          IBM RS6000 stations.
C
C TYPE   : main 
C PURPOSE: second main program for the Hubbard program. The 
C          function G(tau) of program kondo.f is fourier-
C          transformed
C          to yield G(omega), and a new function GO(omega) is
C          produced. New input G0(tau) for program kondo.f 
C          by inverse fourier transform.
C I/O    : parameters in lisaqmc.dat 
C VERSION: 30-Sep-95
C          29-Nov-95 removal of minimal bug in nfourier
C AUTHOR : W. Krauth (krauth@physique.ens.fr)
C COMMENT: Even though FORTRAN is case-insensitive, I have 
C          capitalized all
C          the global variables, i. e. the variables appearing in the
C          COMMON block /global/
C          there is one complication: 
C
C          For historical reasons, and in order to conform both
C          to the QMC-convention, and the one used by other people,
C          Green's functions in programs lisaqmc.f and lisaself.f
C          are defined differently.
C          in lisaqmc.f    :  G(tau)= T<    > 'anti Landau-Lifshitz'
C          everywhere else :  G(tau)=-T<    > 'Landau-Lifshitz'
C
C          
C========+=========+=========+=========+=========+=========+=========+=$
      program lisaself
      include 'lisaqmc.dat' 
      character*80 xyz 
      dimension vdummyup(L), vdummydo(L)
      complex*16 dumcup,dumcdo,gdumup,gdumdo,argu,arguup,argudo, 
     $xi,omega,wfun,ar,cdummy
      logical extended,half 
      external wfun
      Zero=0
      One=1
      Two=2
      xpi=acos(-One)
      xi=cmplx(Zero,One)
cc
cc    open parameter file
cc
      open (unit=1,file='lisaqmc.input',form='formatted',status='old')
C========+=========+=========+=========+=========+=========+=========+=$
c     initial set-up
C========+=========+=========+=========+=========+=========+=========+=$
      read (1,*)
      read (1,*)
      read (1,*)
      read (1,*)
cc
cc    xmut = mu - Uo/2
cc
      read (1,*)Beta,Uo,Vt,xmut,h
      read (1,'(a80)') xyz 
      read (1,'(a80)') xyz 
      read (1,'(a80)') xyz 
      read (1,*) ic,jc
      extended= (abs(H).gt.1.e-7.or.jc.ne.0)
      half= (abs(xmut).lt.1.e-7)
cc
cc     (ic 0:Gauss; 1:square; 2:semicirc; 3:half circle;
cc      jc 0: para; 1: antiferro) 
cc
        open (unit=14,file='lisaqmc.init',
     &            form='formatted',status='old')
        open (unit=17,file='lisaqmc.result',form=
     &           'formatted',status='old')
        call rheader(17)
        call rheader(14)
cc
cc     read in and change sign of function green(tau)
cc     the following is in fact Green0t (which is not needed)
cc
       do 98 i=1,L
          if (extended) then
             read(14,*)dummy1, dummy2
             vdummyup(i)=-dummy1
             vdummydo(i)=-dummy2
          else
             read(14,*)dummy1
             vdummyup(i)=-dummy1
          end if
98     continue
cc
cc     read Is(i) I=1,L, so that it can be put into the file lisaqmc.init
cc     for the next iteration
cc
       do 1748 i=1,L
          read(14,*)Is(i)
1748   continue
       read(14,*)idum
cc
cc     symmetrize the function Green0t(tau) and Green0t(beta-tau)
cc     for the half-filled case
cc
       if (.not.extended.and.half) then
          vdummyup(1)=-One/Two
          do 1789 i=2,L/2
             dummy=vdummyup(i)+vdummyup(L+2-i)
             vdummyup(i)=dummy/Two
             vdummyup(L+2-i)=dummy/Two
1789      continue
       end if
       call nfourier(vdummyup,Green0wup)

       if (extended) call nfourier(vdummydo,Green0wdo) 
       do i=1,L 
          if (extended) then
             read(17,*)dummy1, dummy2
             Greentup(i)=-dummy1
             Greentdo(i)=-dummy2
          else
             read(17,*)dummy1
             Greentup(i)=-dummy1
          end if
       end do
cc
cc     symmetrize the function Greent(tau) and Greent(beta-tau)
cc     for the half-filled case
cc
       if (.not.extended.and.half) then
          do 1848 i=2,L/2
             dummy=Greentup(i)+Greentup(L+2-i)
             Greentup(i)=dummy/Two
             Greentup(L+2-i)=dummy/Two
1848      continue
          Greentup(1)=-One/Two
       end if
       call nfourier(Greentup,Greenwup)
       if (extended)call nfourier(Greentdo,Greenwdo)
       print*,Greenwup(3),'up  qq'
       print*,Greenwdo(3),'do  qq'
cc
cc     update Green0wup, see notes 15/11/91, 22/3/92
cc
       do 818 i=0,Iwmax
          omega=(Two*i+One)*xpi/Beta*xi

          dumcup=One/Greenwup(i)-One/Green0wup(i) 
          arguup=omega+dumcup+xmut+h
          fup=One
          if (imag(arguup).le.0) fup=-One
cc
cc    Gaussian
cc
          if (ic.eq.0) then
             arguup=arguup*fup
             gdumup=-fup*xi*sqrt(xpi)*wfun(arguup)
cc
cc    Square
cc
          else if (ic.eq.1) then
             arguup=(arguup+One)/(arguup-One)
              gdumup=log(arguup)/Two
cc
cc     semi circular
cc
          else if (ic.eq.2) then
             gdumup=-Two/Greenwup(i)
             dumcup=-omega - xmut
          end if
cc
cc        all cases in which the down-Green's function is needed
cc
          if (extended) then

             dumcdo=One/Greenwdo(i)-One/Green0wdo(i) 
             argudo=omega+dumcdo+xmut-h
                fdo=One
                if (imag(argudo).le.0) fdo=-One
                if (ic.eq.0) then
                   argudo=fdo*argudo
                   gdumdo=-fdo*xi*sqrt(xpi)*wfun(argudo)
                   else if (ic.eq.1) then
                   argudo=(argudo+One)/(argudo-One)
                   gdumdo=log(argudo)/Two
                end if
cc
cc        treat separately the antiferromagnetic case (h==0!)
cc
             if (jc.eq.1) then
                if (abs(h).gt.1.e-6) pause 'h=0!'
                dumcup=1./Greenwup(i)-1./Green0wup(i)
                dumcdo=1./Greenwdo(i)-1./Green0wdo(i)
                arguup=omega+dumcup+xmut
                argudo=omega+dumcdo+xmut
                argu=sqrt(arguup*argudo)
                if (ic.eq.0) then
                   f=One
                   if (imag(argu).le.0) f=-One
                   argu=f*argu
                   gdumup=-f*xi*sqrt(xpi)*wfun(argu)*
     &             argudo/argu
                   gdumdo=-f*xi*sqrt(xpi)*wfun(argu)*
     &             arguup/argu
                else if (ic.eq.1) then
                   ar=(argu+One)/(argu-One)
                   gdumup=log(ar)/Two*argudo/argu
                   gdumdo=log(ar)/Two*arguup/argu 
                else if (ic.eq.2) then
                   cdummy=argu**2 - Two
                   gdumup=argudo/argu*(argu - sqrt(cdummy))
                   gdumdo=arguup/argu*(argu - sqrt(cdummy))
                end if
             end if
          end if
cc
cc        now use self-consistency condition
cc        choose one of the following two lines
cc
c         Green0wup(i)=(Green0wup(i)+One/(-dumcup+One/gdumup))/Two
          Green0wup(i)=One/(-dumcup+One/gdumup)
          if (i.lt.10)   print*,Green0wup(i),i,'up'

          if (extended) then
c            Green0wdo(i)=(Green0wdo(i)+One/(-dumcdo+One/gdumdo))/Two
             Green0wdo(i)=One/(-dumcdo+One/gdumdo)
             if (i.lt.10)print*,Green0wdo(i),i,'do'
          end if
818    continue
cc
cc     inverse fourier transform and sign change, in order to
cc     conform to Kondo-definition
cc
       call invfourier(Green0wup,Green0tup) 
       if (extended) call invfourier(Green0wdo,Green0tdo) 
       rewind(14) 
       call wheader(14) 
       do i=1,L 
          if (extended) then
             write(14,'(2f20.10)')-Green0tup(i), -Green0tdo(i)
          else
             write(14,'(f20.10)')-Green0tup(i)
          end if
       end do
       do 1849 i=1,L
          write(14,'(I5)')Is(i)
1849   continue
       write(14,*)idum
C========+=========+=========+=========+=========+=========+=========+=$
c      END OF Iteration 
C========+=========+=========+=========+=========+=========+=========+=$
       write(6,'(a20)')'         Summary - LISASELF'
       write(6,'(a60)')'========================================'
       end
C========+=========+=========+=========+=========+=========+=========+=$
C       PROGRAM: nfourier
C       TYPE   : subroutine
C       PURPOSE: fourier-transform the natural-spline interpolation
C                of function Green(tau) 
C                calculate function Green(omega)
C       I/O    :
C       VERSION: 2-16-92
C                29-Nov-95 removal of minimal bug concerning 
C                          dimension of rindata
C       COMMENT: cf J. Stoer R. Bulirsch, Introduction to numerical
C                analysis (Springer, New York, 1980)
C========+=========+=========+=========+=========+=========+=========+=$
       subroutine nfourier(rindata,coutdata)
       include 'lisaqmc.dat'
       dimension rindata(L)
       dimension rincopy(L+1),a(L),b(L),c(L),d(L), 
     & u(L+1), q(L+1),xm(L+1)
       complex*16 coutdata(0:Iwmax),xi,cdummy,explus,ex
       xi=cmplx(Zero,One)
       xpi=acos(-One)
       delta=Beta/L
       do i=1,L
          rincopy(i)=rindata(i)
       end do
       rincopy(L+1)=-1-rindata(1)
       three=Two+One
       six=2*three
     
cc
cc     spline interpolation:  the spline is given by
cc     G(tau)=a(i) + b(i) (tau-tau_i) + c(i) ( )^2 + d(i) ( )^3
cc     The following formulas are taken directly from  Stoer and
cc     Bulirsch p. 102
cc
       q(1)=Zero
       u(1)=Zero
       do k=2,L
          p=q(k-1)/Two+Two
          q(k)=-One/Two/p
          u(k)=three/delta**2*(rincopy(k+1)+rincopy(k-1)-Two*rincopy(k))
          u(k)=(u(k)-u(k-1)/Two)/p
       end do
       XM(L+1)=0
       do k=L,1,-1
          XM(k)=q(k)*XM(k+1)+u(k)
       end do
cc
cc     The following formulas are taken directly from  Stoer and
cc     Bulirsch p. 98
cc
       do j = 1, L
          a(j)=rincopy(j)
          c(j)= XM(j)/Two
          b(j)=(rincopy(j+1)-rincopy(j))/delta - 
     &       (Two*XM(j)+XM(j+1))*delta/6.
          d(j)=(XM(j+1)-XM(j))/(6.*delta)
       end do
cc
cc     The Spline multiplied by the exponential can now be exlicitely
cc     integrated. The following formulas were obtained using
cc     MATHEMATICA
cc
        do i=0,Iwmax
           om=(Two*i+One)*xpi/Beta
           coutdata(i)=Zero
           do j=1,L
              cdummy=xi*om*delta*j
              explus=exp(cdummy)
              cdummy=xi*om*delta*(j-1)
              ex=exp(cdummy)
              coutdata(i)=coutdata(i) + explus*(
     &         ( -six* d(j) )/om**4 + 
     &         ( Two*xi*c(j) + six*delta*xi*d(j)  )/om**3 +
     &         ( b(j)+ Two*delta*c(j)+ three*delta**2*d(j) )/om**2 +
     &         (- xi*a(j) - delta*xi*b(j) - delta**2*xi*c(j) -
     &         delta**3*xi*d(j))/om)
 
              coutdata(i) = coutdata(i) + ex*(
     &        six*d(j)/om**4 - Two*xi*c(j)/om**3 
     &        -b(j)/om**2 + xi*a(j)/om)
           end do
        end do
        end
C========+=========+=========+=========+=========+=========+=========+=$
C       PROGRAM: invfourier
C       TYPE   : subroutine
C       PURPOSE: inverse fourier transform
C                Greent, Greenw use physical definition
C                Greent(i)=G((i-1)*deltau) for i=1,...,L
C                Greenw(n)=G(i w_n), for n=0,L/2-1
C                       w_n=(2*n+1)pi/beta
C                Symmetry property: 
C                G(iw_(-n)=G(iw_(n-1))*
C                coupled to the impurity
C       I/O    :
C       VERSION: 6-16-92
C       COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
       subroutine invfourier(cindata,routdata)
       include 'lisaqmc.dat'
       dimension routdata(L)
       complex*16 cindata(0:Iwmax),cdummy
       xpi=acos(-One) 
       do 1 i=1,L
       routdata(i)=Zero
          tau=(i-1)*beta/real(L)
          do 2 j=0,Iwmax
               om=mod((2*j+One)*xpi/Beta*tau,2*xpi)
               cdummy=cmplx(Zero,One)*om
               dummy=cindata(j)*exp(-cdummy)
         routdata(i)=routdata(i)+Two/beta*dummy
2         continue
1      continue
cc
cc     special treatment for tau=0
cc
       routdata(1)=-One/Two+routdata(1)
       end
C========+=========+=========+=========+=========+=========+=========+=$
C       PROGRAM: rheader.f  
C       TYPE   : subroutine
C       PURPOSE: read the program's header from unit k 
C       I/O    :
C       VERSION: 6-16-92
C       COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine rheader(k) 
        do 1 i=1,10
        read(k,*)
1       continue
        end 
C========+=========+=========+=========+=========+=========+=========+=$
C       PROGRAM: wheader.f  
C       TYPE   : subroutine
C       PURPOSE: write the program's header onto standard output 
C       I/O    : 
C       VERSION: 6-16-92
C       COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine wheader(k) 
        include 'lisaqmc.dat'
        character *80 xyz
        write(k,'(a55)')'========================================'
        write(k,'(a55)')'lisaself   : version 30-Sep-95 '
        write(k,'(a55)')'========================================'
        rewind(1)
        write(k,'(4(a6,I4))') 'L=',L
        read(1,'(a60)')xyz
        read(1,'(a60)')xyz
        read(1,'(a60)')xyz
        do 3 i=1,6
        read(1,'(a60)')xyz
        write(k,'(a60)')xyz
3       continue 
        rewind(1)
        end
C========+=========+=========+=========+=========+=========+=========+=$
C       PROGRAM: wfun.f
C       TYPE   : subroutine
C       PURPOSE: calculate exponentially  rescaled error function w 
C       I/O    :
C       VERSION: 6-16-92
C       COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
       function wfun(z)
       complex*16 z,wfun
       double precision xi,yi,u,v
       logical flag
       xi=dreal(z)
       yi=imag(z)
       call wofz(xi,yi,u,v,flag)
       wfun=dcmplx(u,v)
        end  
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
