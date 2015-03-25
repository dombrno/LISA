C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: lisaqmc.f  
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
C          The programs have  been thoroughly tested on SUN Sparks, 
C          HP 900, IBM RS6000 stations.
C
C TYPE   : main 
C PURPOSE: program for qmc-simulation of Anderson impurity
C          problem
C I/O    : cf file README_lisaqmc
C VERSION: 30-Sep-95 
C AUTHOR : W. Krauth (krauth@physique.ens.fr)
C COMMENT: Even though FORTRAN is case-insensitive, I have capitalized 
C          all the global variables, i. e. the variables appearing 
C          in the common/global/ block (cf file lisaqmc.dat).
C          At the end of the program, a number of SLATEC routines
C          have been appended. You may not want to  print these.
C========+=========+=========+=========+=========+=========+=========+=$
      program lisaqmc
      include 'lisaqmc.dat' 
      real ranw
      logical stochastic
      dimension greenm(-L+1:L-1),greenm2(-L+1:L-1),greenupsum(0:L-1)
      dimension greendosum(0:L-1)
      integer itau(0:L+1),iss(L)
      nran(i)=mod(int(i*ranw(Idum)),i) + 1 
      Zero=0
      One=1
      Two=2
      stochastic=.true.
cc      if (L.le.18) stochastic=.false.
      open (unit=1,file='lisaqmc.input',form='formatted',status='old')
cc     
cc    write header onto standard output 
cc
      call wheader(6) 
c======================================================================
c     initial set-up
c======================================================================
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      read (1,*)Beta,U,xmut,H
cc
cc    cosh(Xlambda)=exp(deltau*U/2)     cf. eq. (\ref{hirschdecoup})
cc    cf. eq. 119
cc
      deltau=Beta/real(L)
      dummy=exp(deltau*U/Two)
      Xlambda=log(dummy+sqrt(dummy**2-One))
      read (1,*)
      read (1,*)maxt
      read (1,*) 
      read (1,*)ic,jc
      if (jc.eq.0.and.abs(H).lt.1.e-7) then
         Paramagnet=.true.
      else
         Paramagnet=.false.
      end if    
      open (unit=14,file='lisaqmc.init',form='formatted',status='old')
      print*,' begin initi paa', paramagnet
      call initial 
      print*,' end init'
      print*, ' green0up'
       do i=1,L
             write(*,'(8f10.5)')(green0up(i,j),j=1,8)
       end do
       print*, ' green0do'
       do i=1,L
             write(*,'(8f10.5)')(green0do(i,j),j=1,8)
       end do
      naccept=0
      isign=1
      do 810 i=0,L-1
         greenupsum(i)=Zero
         greendosum(i)=Zero
810   continue
      dsum=Zero
      upsum=Zero
      dosum=Zero
      iitime=0
      iter=0
c***********************************************************************
c		start simulation  
c***********************************************************************
cc
cc    SIMULATION BY MONTE CARLO   (stochastic.eq.true)
cc
      print*, ' stochastic', stochastic
      if (stochastic) then 
         nfastup=100
         do 1000 iii=1,1000000000
cc
cc       we make nfastup sweeps with the fast update (eq. (\ref{fastupdate}))
cc       until checking with subroutine 'update' (eq. 
cc       whether precision hasn't 
cc       deteriorated. If that is the case, we make nfastup smaller, 
cc       otherwise bigger.
cc       The MC algorithm is driven by a heat bath algorithm.
cc       Acceptance is decided according to Eqs. 134, 
cc       with the ratios of determinants given by eq.131
cc       The fast update algorithm relates to eq. 130.
cc
         do 2000 kkk=1,nfastup
         iter=iter+1
         if (iter.gt.maxt) goto 3120
         do 1900 kk=1,L
            k=nran(L)
c
c           try flipping spin k
c
            isignnew=abs(detrat(k))/detrat(k)
            if (isignnew.ne.isign) then
cc                print*,'signchange',isignnew
               isign=isignnew
            end if
            dummy=abs(detrat(k))
            if (ranw(Idum).lt.dummy/(One+dummy)) then
cc
cc             accept flip, 
cc
               call record(k)
               naccept=naccept+1
            end if 
1900        continue
cc
cc          end of sweep: The green's function,
cc          which are functions of imaginary time tau,
cc          are "averaged", in order to reduce the statistical
cc          noise, by forcing the translation invariance in imaginary
cc          time of the physical Green's functions.
cc          Cf remark, last paragraph of section VI A 1 c,
cc          left column of p38
cc
            iitime=iitime+1
            do 192 idel=-L+1,L-1
               isum=0
               greenm(idel)=Zero
               greenm2(idel)=Zero
               inumb=min(L-idel,L)-max(1,1-idel)+1
               do 193 i=max(1,1-idel),min(L-idel,L)
                  dummy=inumb	
                  greenm(idel)=greenm(idel)+greenup(i+idel,i)/dummy
                  greenm2(idel)=greenm2(idel)+greendo(i+idel,i)/dummy
	          isum=isum+1
193            continue 
192         continue
            do 191 i=1,L
               dummy=L
cc            dsum is used to track the number of doubly
cc            occupied sites. An average over tau is done, but 
cc            this point is not discussed in the review.
               dsum=dsum+greenup(i,i)*greendo(i,i)/dummy
191         continue
cc         dosum/upsum are the down and up density
cc         Needs tgo be justified...
            dosum=dosum+greenm2(0)
            upsum=upsum+greenm(0)
            do 818 i=0,L-1
               greenupsum(i)=greenupsum(i)+greenm(i)
               greendosum(i)=greendosum(i)+greenm2(i)
818         continue
2000     continue
cc
cc       update greens functions from scratch
cc
         diff=-100000
         call update
         do 887 i=1,L
            do 887 j=1,L
               diff=max(diff,abs(Greenup(i,j)-Gnewup(i,j)),
     &         abs(Greendo(i,j)-Gnewdo(i,j)))
               Greenup(i,j)=Gnewup(i,j)
               Greendo(i,j)=Gnewdo(i,j)
887      continue
         if (diff.gt.0.0005) then
            print*, nfastup,diff,'   nfastup diff'
            nfastup=max(1,nfastup/2)
         end if
         if (diff.lt.0.0005) nfastup=nfastup*2
1000     continue
3120  continue
      else
cc
cc       exact enumeration for L <= 18  (stochastic=.false.)
cc       initialization of Gray code variables (itau is a variable used
cc       for the construction of the code cf. (Reingold  et al 1977)
cc
         do 3001 j=0,L
            itau(j)=j+1
3001     continue
         do 3002 j=1,L
            is(j)=-1
3002     continue
         call update

         do 3199 n=1,L
            do 3199 m=1,L
               Greenup(n,m)=Gnewup(n,m)
               Greendo(n,m)=Gnewdo(n,m)
3199     continue
c        print*,' Greenup'
c        do i=1,L
c           write(*,'(8f10.5)')(Greenup(i,j),j=1,8)
c        end do
c        print*,' Greendo'
c        do i=1,L
c           write(*,'(8f10.5)')(Greendo(i,j),j=1,8)
c        end do
         do 910 i=0,L-1
            greenupsum(i)=greenup(i+1,1)
   	    greendosum(i)=greendo(i+1,1)
910      continue
         dsum=greenup(1,1)*greendo(1,1)
         upsum=greenup(1,1)
         dosum=greendo(1,1)
         det=One
cc
cc       Note that det=One gives the determinant of the initial configuration
cc       only up to a multiplicative factor. If you are interested in the 
cc       numerical value of the partition function, or the free energy,
cc       you will have to replace the above line by  the following:
cc
cc       det=One/(determinant(greenup,L)*determinant(greendo,L))
cc
cc       (the external function determinant is provided below, but not
cc       actually used
cc
         partition=det
         itup=0
         do 3000 iii=1,2**L-1
            print*, iii
            itup=itup+1
cc
cc          use Gray code to calculate index of spin which has to be flipped
cc
            k=itau(0)
            itau(k-1)=itau(k)
            itau(k)=k+1
            if (k.ne.1) itau(0)=1
cc
cc          flip spin k
cc
            det=detrat(k)*det
	    partition=partition+det
            isignnew=abs(detrat(k))/detrat(k)
            if (isignnew.ne.isign) then
cc               print*,'signchange',isignnew
               isign=isignnew
            end if
            call record(k)
            do 920 i=0,L-1
	       greenupsum(i)=greenupsum(i)+greenup(i+1,1)*det
	       greendosum(i)=greendosum(i)+greendo(i+1,1)*det
920         continue
            dsum=dsum + greenup(1,1)*greendo(1,1)*det
            upsum=upsum + greenup(1,1)*det
            dosum=dosum + greendo(1,1)*det
cc
cc          check that precision of determinant is not degraded
cc          (the last configuration of M.C. spins is again recreated from
cc          (-1,-1,...-1) and the determinant is recomputed).
cc
            if (itup.gt.100) then
               itup=0
               detcheck=One
               do 3007 j=1,L
                  iss(j)=is(j)
                  is(j)=-1
3007           continue
               call update
               do 3799 n=1,L
                  do 3799 m=1,L
                     Greenup(n,m)=Gnewup(n,m)
                     Greendo(n,m)=Gnewdo(n,m)
3799           continue
               do 881 k=1,L
                  if (iss(k).ne.-1) then
                     detcheck=detcheck*detrat(k)
                     call record(k)
                  end if
881            continue
               det=detcheck
            end if
3000     continue
      end if
c====================================================================
c                      END OF SIMULATION
c====================================================================
      if (stochastic) then 
         facnorm=iitime
      else
         facnorm=partition
      end if
      open (unit=17,file='lisaqmc.result',form='formatted',status='new')
      call wheader(17)
      if (Paramagnet) then
         do 717 i=0,L-1
         write(17,'(f20.10)')(greenupsum(i)+greendosum(i))/facnorm/Two
717      continue
      else 
         do 718 i=0,L-1
            write(17,'(2f20.10)')greenupsum(i)/facnorm,
     &      greendosum(i)/facnorm
718   continue
      end if 
      write(*,'(a20,f15.7)')' prob double occ=',dsum/facnorm
      write(*,'(a20,f15.7)')' density    up  =',
     & One-upsum/facnorm 
      write(*,'(a20,f15.7)')' density    do  =',
     &  One-dosum/facnorm 
      write(6,'(a60)')'========================================'
      write(6,'(a20)')'         Summary'
      write(6,'(a20,I15)')'Length of Simul = ',iter
      write(6,'(a20,f15.4)')'Acceptance Prob = ',
     &    naccept/real(maxt*L)
      write(6,'(a60)')'========================================'
      open (unit=15,file='lisaqmc.end',form='formatted',status='new')
      call wheader(15)
      if (Paramagnet) then
         do 616 i=1,L
            write(15,'(f20.10)')Green0up(i,1)
616      continue
      else
         do 618 i=1,L
            write(15,'(f20.10,f20.10)')Green0up(i,1),Green0do(i,1)
618      continue
      end if
      do 617 i=1,L
         write(15,'(i5)')Is(i)
617   continue
cc
cc    write seed for next run onto file 15
cc
      write(15,'(I10)')int(-ranw(Idum)*100000)
      end
C========+=========+=========+=========+=========+=========+=========+=$
C     PROGRAM: detrat.f  
C     TYPE   : function
C     PURPOSE: calculate the ratio of the new and
C              old determinants (cf. eq. (\ref{detrat}) i.e. eq. 131
C     I/O    : 
C     VERSION: 30-Sep-95
C     COMMENT:
C========+=========+=========+=========+=========+=========+=========+=$
      function detrat(k)
      include 'lisaqmc.dat'
      rup=One+(One-Greenup(k,k))*(exp(-Two*Xlambda*real(Is(k)))-One)
      rdo=One+(One-Greendo(k,k))*(exp( Two*Xlambda*real(Is(k)))-One)
      detrat=rup*rdo
      end
C========+=========+=========+=========+=========+=========+=========+=$
C     PROGRAM: initial.f  
C     TYPE   : subroutine 
C     PURPOSE: read in initial configuration of bath Green's function
C              and of Ising spins, expand G(i-j) into matrix G(i,j). 
C              invoke subroutine Update to calculate 
C              Green's function for the initial choice of
C              Ising spins.
C     I/O    : 
C     VERSION: 30-Sep-95
C     COMMENT:  
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine initial()
      include 'lisaqmc.dat'
      dimension gtempup(-L+1:L-1)
      dimension gtempdo(-L+1:L-1)
      print*,L,paramagnet,' L paramagnet'
cc
cc    read in G(sigma=0) from file lisaqmc.init
cc
      call rheader(14)
      if (Paramagnet) then
         do 177 i=0,L-1
            read(14,*) gtempup(i)
            gtempdo(i)=gtempup(i)
177      continue
      else
         do 171 i=0,L-1
            read(14,*) gtempup(i),gtempdo(i)
171      continue
      end if
      print*,' 171'
cc
cc    reflection of G to calculate Greens function for neg. arguments	
cc
      do 123 i=1,L-1
         gtempup(-i)=-gtempup(L-i)
         gtempdo(-i)=-gtempdo(L-i)
123   continue

      print*,' 123'
cc
cc    choice of spins
cc
      do 11 n=1, L
         read(14,*) Is(n)
11    continue
      print*,' 11'
      read(14,*)Idum
cc
cc    calculation of Green's function
cc    update puts results in array Gnewup, Gnewdo, which is then
cc    transcribed into array Greenup, Greendo (use translation invariance).
cc
      dummy=L
      deltau=Beta/dummy
      do 99 i=1,L
         do 99 j=1,L
            Green0up(i,j)=gtempup(i-j)
            Green0do(i,j)=gtempdo(i-j)
99    continue
      do 2 n=1,L
         do 2 m=1,L
            Greenup(n,m)=Green0up(n,m)
            Greendo(n,m)=Green0do(n,m)
2     continue
      call update
      do 3 n=1,L
         do 3 m=1,L
            Greenup(n,m)=Gnewup(n,m)
            Greendo(n,m)=Gnewdo(n,m)
3     continue
      end
C========+=========+=========+=========+=========+=========+=========+=$
C     PROGRAM: ranw.f
C     TYPE   : real function
C     PURPOSE: produce uniformly distributed random numbers
C              following the algorithm of Mitchell and Moore
C     I/O    :
C     VERSION: 30-Sep-95
C     COMMENT: cf. D. E. Knuth, Seminumerical Algorithms, 2nd edition
C              Vol 2 of  The Art of Computer Programming (Addison-Wesley,
C              1981) pp 26f. (Note: the procedure ran3 in 
C              W. H. Press et al,  Numerical
C              Recipes in FORTRAN, 2nd edition (Cambridge University
C              Press 1992)  is based on the same algorithm). 
C              I would suggest that you make sure for yourself that 
C              the quality of the random number generator is sufficient, 
C              or else replace it!
C========+=========+=========+=========+=========+=========+=========+=$
      real function ranw(idum)
      Parameter (Mbig=2**30-2, Xinvers=1./Mbig)
      Integer IX(55)
      data ibit/ 1/
      save
      if (ibit.ne.0) then
         ibit=0
cc
cc       fill up the vector ix with some random integers, which are
cc       not all even
cc       
         if (idum.eq.0) pause 'use nonzero value of idum'
         idum=abs(mod(idum,Mbig))
         ibit=0
         Ix(1)=871871
         Do i=2,55
            Ix(i)=mod(Ix(i-1)+idum,Ix(i-1))
            Ix(i)=max(mod(Ix(i),Mbig),idum)
         end do
         j=24
         k=55
cc
cc       warm up the generator     
cc
         do i=1,1258
            Ix(k)=mod(Ix(k)+Ix(j),Mbig)
            j=j-1
            if (j.eq.0) j=55 
            k=k-1
            if (k.eq.0) k=55 
         end do
      end if
cc
cc    this is where execution usually starts:
cc
      Ix(k)=mod(Ix(k)+Ix(j),Mbig)
      j=j-1
      if (j.eq.0) j=55 
      k=k-1
      if (k.eq.0) k=55 
      ranw=Ix(k)*Xinvers 
      end
C========+=========+=========+=========+=========+=========+=========+=$
C       PROGRAM: record  
C       TYPE   : subroutine 
C       PURPOSE: record changes of accepted move on the Green's function
C                (cf  eq. (\ref{fastupdate}))
C       I/O    : 
C       VERSION: 30-9-95
C       COMMENT: k is the index of the spin which has been flipped
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine record(k)
        include 'lisaqmc.dat'
cc
cc    update Green's function (implementation of  eq. (\ref{fastupdate}))
cd    i.e. eq. 130.
cc
        del=-Two*Xlambda*Is(k)
        do 1 i=1,L
           idel=0
           if (i.eq.k) idel=1
           do 1 j=1,L
              Gnewup(i,j)=Greenup(i,j)+(Greenup(i,k)-idel)*
     $        (exp(del)-1.)/(1.+(1.-Greenup(k,k))*(exp(del)-One))*
     $        Greenup(k,j)
              Gnewdo(i,j)=Greendo(i,j)+(Greendo(i,k)-idel)*
     $        (exp(-del)-One)/(One+(One-Greendo(k,k))*(exp(-del)-One))*
     $        Greendo(k,j)
1       continue
        do 2 j=1,L
           do 2 i=1,L
              Greenup(i,j)=Gnewup(i,j)
              Greendo(i,j)=Gnewdo(i,j)
2       continue  
cc
cc   update spin
cc
        Is(k)=-Is(k)
        end 

C========+=========+=========+=========+=========+=========+=========+=$
C     PROGRAM: update.f  
C     TYPE   : subroutine 
C     PURPOSE: calculate the Green's function 
C              for a given configuration of spins 
C              (in vector Is) from the Green's function
C              for spins set equal to zero  (eq. (\ref{inversion})) i.e. eq. 128, 122
C     I/O    : 
C     VERSION: 30-Sep-95
C     COMMENT: can be used to initialize run
C                     (subroutine initial),
C              or to check for deterioration of precision 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine update
      include 'lisaqmc.dat'
      dimension a(L,L),b(L,L),ainv(L,L),binv(L,L)
cc
cc    calculate the matrix a=1-(g-1)(exp(v')-1)
cc    Cf eq. 128, 122
cc
      do 2 i=1,L
         do 3 j=1,L
            a(i,j)=-Green0up(i,j)*(exp(Xlambda*real(Is(j)))-One)
            b(i,j)=-Green0do(i,j)*(exp(-Xlambda*real(Is(j)))-One)
3        continue
         a(i,i)=1-(Green0up(i,i)-One)*(exp(Xlambda*real(Is(i)))-One)
         b(i,i)=1-(Green0do(i,i)-One)*(exp(-Xlambda*real(Is(i)))-One)
2     continue
      call inverse(a,ainv)
      call inverse(b,binv)
      do 4 i=1,L
         do 4 j=1,L
            suma=Zero
            sumb=Zero
            do 6 k=1,L
               suma=suma+ainv(i,k)*Green0up(k,j)
               sumb=sumb+binv(i,k)*Green0do(k,j)
6           continue
         Gnewup(i,j)=suma
         Gnewdo(i,j)=sumb
4     continue
      end
c======================================================================
C PROGRAM: inverse.f
C TYPE   : subroutine
C PURPOSE: calculate inverse of matrix
C I/O    :
C VERSION: 30-Sep-95
C COMMENT: Here we use the (public domain) SLATEC routines dgeco and
C          dgedi (which perform the Gaussian elimination) + dependencies.
C          These routines can be obtained, e.g., from netlib
C          (to learn about netlib, send an otherwise empty
C          message to netlib@research.att.com
C          containing 'send index' in the subject header, 
C          on WWW, look under the address 
C          http://netlib.att.com/netlib/master/readme.html).
C========+=========+=========+=========+=========+=========+=========+=$
      Subroutine inverse(a,y)
      include 'lisaqmc.dat'
      dimension a(L,L),y(L,L)
      dimension z(L),ipvt(L)
      do 1 i=1,L
         do 2 j=1,L
              y(i,j)=a(i,j)
2        continue
1     continue
      call dgeco(y,L,L,ipvt,rcond,z)
cc
cc    we only want the Inverse, so set job = 01
cc
      job=01
      call dgedi(y,L,L,ipvt,det,z,job)
      end
c======================================================================
C       PROGRAM: rheader.f  
C       TYPE   : subroutine
C       PURPOSE: read the program's header from unit k 
C       I/O    :
C       VERSION: 30-Sep-95
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
C       VERSION: 30-Sep-95
C       COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine wheader(k) 
        include 'lisaqmc.dat'
        character *80 xyz
        write(k,'(a55)')'========================================'
        write(k,'(a55)')' Lisaqmc  simulation: rev. 06/30/95'
        write(k,'(a55)')'========================================'
        rewind(1)
        write(k,'(4(a6,I4))') 'L=',L
        read(1,*)
        read(1,*)
        read(1,*)
        do 3 i=1,6
        read(1,'(a60)')xyz
        write(k,'(a60)')xyz
3       continue 
        rewind(1)
        end
C========+=========+=========+=========+=========+=========+=========+=$
C       PROGRAM: dasum.f daxpy.f  ddot.f dgeco.f dgedi.f dgefa.f 
C                dscal.f dswap.f idamax.f
C       TYPE   : collection of subroutines 
C       PURPOSE: calculate inverse and determinant (look at 
C                subroutine dgedi.f) 
C       I/O    :
C       VERSION: 30-Sep-95
C       COMMENT: the following subroutines are a bunch of
C                functions obtained from the slatec library
C                at Netlib, which allow the calculation of 
C                inverse and determinant.
C                You can replace these programs by the 
C                corresponding routines of your favorite library,
C                e.g. Numerical Recipes (which is not in the
C                public domain).
C                Notice that we are using the double precision 
C                versions of the programs.
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
*DECK DASUM
      DOUBLE PRECISION FUNCTION DASUM (N, DX, INCX)
C***BEGIN PROLOGUE  DASUM
C***PURPOSE  Compute the sum of the magnitudes of the elements of a
C            vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A3A
C***TYPE      DOUBLE PRECISION (SASUM-S, DASUM-D, SCASUM-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SUM OF MAGNITUDES OF A VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    DASUM  double precision result (zero if N .LE. 0)
C
C     Returns sum of magnitudes of double precision DX.
C     DASUM = sum from 0 to N-1 of ABS(DX(IX+I*INCX)),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DASUM
      DOUBLE PRECISION DX(*)
      INTEGER I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  DASUM
      DASUM = 0.0D0
      IF (N .LE. 0) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DASUM = DASUM + ABS(DX(IX))
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 6.
C
   20 M = MOD(N,6)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DASUM = DASUM + ABS(DX(I))
   30 CONTINUE
      IF (N .LT. 6) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DASUM = DASUM + ABS(DX(I)) + ABS(DX(I+1)) + ABS(DX(I+2)) +
     1          ABS(DX(I+3)) + ABS(DX(I+4)) + ABS(DX(I+5))
   50 CONTINUE
      RETURN
      END
*DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A7
C***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DAXPY
      DOUBLE PRECISION DX(*), DY(*), DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DDOT
C***PURPOSE  Compute the inner product of two vectors.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A4
C***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DDOT
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +
     1              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DGECO
      SUBROUTINE DGECO (A, LDA, N, IPVT, RCOND, Z)
C***BEGIN PROLOGUE  DGECO
C***PURPOSE  Factor a matrix using Gaussian elimination and estimate
C            the condition number of the matrix.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGECO-S, DGECO-D, CGECO-C)
C***KEYWORDS  CONDITION NUMBER, GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGECO factors a double precision matrix by Gaussian elimination
C     and estimates the condition of the matrix.
C
C     If  RCOND  is not needed, DGEFA is slightly faster.
C     To solve  A*X = B , follow DGECO by DGESL.
C     To compute  INVERSE(A)*C , follow DGECO by DGESL.
C     To compute  DETERMINANT(A) , follow DGECO by DGEDI.
C     To compute  INVERSE(A) , follow DGECO by DGEDI.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an INTEGER vector of pivot indices.
C
C        RCOND   DOUBLE PRECISION
C                an estimate of the reciprocal condition of  A .
C                For the system  A*X = B , relative perturbations
C                in  A  and  B  of size  EPSILON  may cause
C                relative perturbations in  X  of size  EPSILON/RCOND .
C                If  RCOND  is so small that the logical expression
C                           1.0 + RCOND .EQ. 1.0
C                is true, then  A  may be singular to working
C                precision.  In particular,  RCOND  is zero  if
C                exact singularity is detected or the estimate
C                underflows.
C
C        Z       DOUBLE PRECISION(N)
C                a work vector whose contents are usually unimportant.
C                If  A  is close to a singular matrix, then  Z  is
C                an approximate null vector in the sense that
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DASUM, DAXPY, DDOT, DGEFA, DSCAL
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGECO
      INTEGER LDA,N,IPVT(*)
      DOUBLE PRECISION A(LDA,*),Z(*)
      DOUBLE PRECISION RCOND
C
      DOUBLE PRECISION DDOT,EK,T,WK,WKM
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  DGECO
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = MAX(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(A(K,K))) GO TO 30
            S = ABS(A(K,K))/ABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/ABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/ABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 150
            S = ABS(A(K,K))/ABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
*DECK DGEDI
      SUBROUTINE DGEDI (A, LDA, N, IPVT, DET, WORK, JOB)
C***BEGIN PROLOGUE  DGEDI
C***PURPOSE  Compute the determinant and inverse of a matrix using the
C            factors computed by DGECO or DGEFA.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D3A1, D2A1
C***TYPE      DOUBLE PRECISION (SGEDI-S, DGEDI-D, CGEDI-C)
C***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGEDI computes the determinant and inverse of a matrix
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        WORK    DOUBLE PRECISION(N)
C                work vector.  Contents destroyed.
C
C        JOB     INTEGER
C                = 11   both determinant and inverse.
C                = 01   inverse only.
C                = 10   determinant only.
C
C     On Return
C
C        A       inverse of original matrix if requested.
C                Otherwise unchanged.
C
C        DET     DOUBLE PRECISION(2)
C                determinant of original matrix if requested.
C                Otherwise not referenced.
C                Determinant = DET(1) * 10.0**DET(2)
C                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
C                or  DET(1) .EQ. 0.0 .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains
C        a zero on the diagonal and the inverse is requested.
C        It will not occur if the subroutines are called correctly
C        and if DGECO has set RCOND .GT. 0.0 or DGEFA has set
C        INFO .EQ. 0 .
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, DSWAP
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGEDI
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),DET(2),WORK(*)
C
      DOUBLE PRECISION T
      DOUBLE PRECISION TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C***FIRST EXECUTABLE STATEMENT  DGEDI
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (ABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (ABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
*DECK DGEFA
      SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
C***BEGIN PROLOGUE  DGEFA
C***PURPOSE  Factor a matrix using Gaussian elimination.
C***LIBRARY   SLATEC (LINPACK)
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
C***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGEFA factors a double precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DSCAL
      SUBROUTINE DSCAL (N, DA, DX, INCX)
C***BEGIN PROLOGUE  DSCAL
C***PURPOSE  Multiply a vector by a constant.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A6
C***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSCAL
      DOUBLE PRECISION DA, DX(*)
      INTEGER I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  DSCAL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END
*DECK DSWAP
      SUBROUTINE DSWAP (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DSWAP
C***PURPOSE  Interchange two vectors.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A5
C***TYPE      DOUBLE PRECISION (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I)
C***KEYWORDS  BLAS, INTERCHANGE, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DX  input vector DY (unchanged if N .LE. 0)
C       DY  input vector DX (unchanged if N .LE. 0)
C
C     Interchange double precision DX and double precision DY.
C     For I = 0 to N-1, interchange  DX(LX+I*INCX) and DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSWAP
      DOUBLE PRECISION DX(*), DY(*), DTEMP1, DTEMP2, DTEMP3
C***FIRST EXECUTABLE STATEMENT  DSWAP
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP1 = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP1
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 3.
C
   20 M = MOD(N,3)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   30 CONTINUE
      IF (N .LT. 3) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP1 = DX(I)
        DTEMP2 = DX(I+1)
        DTEMP3 = DX(I+2)
        DX(I) = DY(I)
        DX(I+1) = DY(I+1)
        DX(I+2) = DY(I+2)
        DY(I) = DTEMP1
        DY(I+1) = DTEMP2
        DY(I+2) = DTEMP3
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DTEMP1 = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP1
   70 CONTINUE
      RETURN
      END
*DECK IDAMAX
      INTEGER FUNCTION IDAMAX (N, DX, INCX)
C***BEGIN PROLOGUE  IDAMAX
C***PURPOSE  Find the smallest index of that component of a vector
C            having the maximum magnitude.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A2
C***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  IDAMAX
      DOUBLE PRECISION DX(*), DMAX, XMAG
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF (N .LE. 0) RETURN
      IDAMAX = 1
      IF (N .EQ. 1) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increments not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DMAX = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
        XMAG = ABS(DX(IX))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increments equal to 1.
C
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
        XMAG = ABS(DX(I))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
   30 CONTINUE
      RETURN
      END
