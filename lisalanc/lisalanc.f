C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: lisalanc.f  
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
C          You are kindly asked to cite the paper (and, if applicable,
C          the original works) if you use this program.
C
C          The present algorithm was first described in:
C          M. Caffarel and W. Krauth, Phys. Rev. Let. 72, 1545 (1994)
C
C          The programs have  been thoroughly tested on SUN Sparks,
C          HP 900, IBM RS6000 stations.
C
C TYPE   : main 
C PURPOSE: self-consistency using Anderson diag. Lanczos program
C I/O    : cf file README_lisalanc for details on code and on
C          input-output
C VERSION: 30-Sep-95
C AUTHOR : W. Krauth (krauth@physique.ens.fr). Important parts of the 
C          code were originally written by M. Caffarel.
C COMMENT: 
C          The Lanczos algorithm presented here is as simple as
C          possible (even though all Lanczos methods which exclusively
C          compute the groundstate properties are simple)
C          cf. Lin and Gubernatis (1993). The present implementation 
C          is particularly simple in that the non-zero elements 
C          of the hamiltonian matrix are explicitly written down
C          in a vector (and dumped onto external files). 
C
C          You are invited to check the present code with the 
C          finite temperature version (lisadiag.f). A convenient way
C          to do that is to go to around line 271, and put id=2.
C          The result should agree to about machine precision.
C
C          It is impossible to present a code which contains all
C          the sophistications of the versions in actual use.
C          In this version, we only calculate thermodynamic Green's
C          functions for the single-band Hubbard model, in the 
C          paramagnetic phase. Generalization of the program to more
C          complicated cases was actually performed. The necessary 
C          modifications are straightforward.
C
C          Even though FORTRAN is case-insensitive, I have capitalized
C          exactly the global variables, i. e. the variables appearing
C          in the common/global/ block (cf file lisalanc.dat).
C
C          At the end of the program, two SLATEC routines
C          have been appended. You may not want to  print these.
C========+=========+=========+=========+=========+=========+=========+=$
      program lisalanc
      include 'lisalanc.dat' 
      character*80 xyz 
      complex*16 omega, adummy, cdummy1
      Zero=0
      One=1
      Two=2
      Pi=acos(-One)
      Xi=cmplx(Zero,One)
      Half=.false.
cc
cc    open parameter file
cc
      open (unit=1,file='lisalanc.input',form='formatted',status='old')
C========+=========+=========+=========+=========+=========+=========+=$
C     initial set-up
C========+=========+=========+=========+=========+=========+=========+=$
      read (1,'(a80)') xyz 
      read (1,'(a80)') xyz 
      read (1,'(a80)') xyz 
      read (1,'(a80)') xyz 
c
cc    xmut = mu - epsd0 -U/Two 
cc
      read (1,*)Beta,U,xmu,Cutoff
      xmut=xmu-U/Two
      open (15,file='lisalanc.andpar',form='formatted',status='old')
      read (1,'(a80)') xyz 
cc
cc    the following line without importance for the Lanczos calculation
cc
      read (1,*) itaucalc,jcut
      open (unit=17,file='lisalanc.green',form='formatted',
     &  status='new')
      open (unit=18,file='lisalanc.diff',form='formatted',
     &  status='new')
cc    open andpar file
      call initial()
      print*,' initial termine'
cc
cc
cc    calculate GwAnderson from the parameters of the Anderson model
cc 
      do i=0,iwmax
         om=(Two*i+One)*Pi/Beta
         call calcg0(om,cdummy1)
         g0wand(i)=cdummy1
      end do
      call diag
cc
cc    calculate G0wnew_Hubbard
cc
      do 818 i=0,Iwmax
         omega=(Two*i+One)*Pi/Beta*Xi
cc
cc       update matrix Gp^-1 : xmu is the true chemical potential
cc                            halffilled: xmu=U/Two
cc
          
         adummy=omega-Epsk(1)-U/Two-Gw(i)/Two
cc
cc       one of the possible ways to write the self-consistency:
cc
         G0w(i)=(One/adummy+G0wand(i))/Two
818   continue
cc
cc    fit an Anderson model G0_Anderson from G0w_Hubbard
cc
      Nitermax=4000
c     call search(chi2,Nitermax)
cc
cc    replace G0_Hubbard by G0_Anderson
cc
      rewind(15)
      call wheader(15)
      write(15,'(9a)')' Eps(k)'
      do i=2,Ns
         write(15,*)Epsk(i)
      end do
      write(15,'(9a)')' V(k)'
      do i=1,Ns-1
         write(15,*)V(i)
      end do
C========+=========+=========+=========+=========+=========+=========+=$
C HERE IS SOME SAMPLE OUTPUT: THERE IS MUCH MORE YOU MAY WANT TO 
C PLOT (for example time-dependent Green's functions, densities of
C states, correlation functions ...
C========+=========+=========+=========+=========+=========+=========+=$
      write(17,'(a18,I4,3f7.2)')' TitleText: G(w)  ',
     %Ns,Beta,U,xmu
      do i=0,50
         om=(2*i+1)*Pi/Beta
         write(17,'(2f17.10)')om,real(Gw(i))
      end do
      write(17,*) 
      do i=0,50
         om=(2*i+1)*Pi/Beta
         write(17,'(2f17.10)')om,imag(Gw(i))
      end do
      write(17,*) 
      write(18,'(a26,I4,3f7.2)')' TitleText: G0(w) G0wand  ',
     %Ns,Beta,U,xmu
      write(18,'(a11)')'"Re(G0w)"'
      do i=0,50
         om=(2*i+1)*Pi/Beta
         write(18,'(2f17.10)')om,real(G0w(i))
      end do
      write(18,*) 
      write(18,'(a11)')'"Im(G0w)"'
      do i=0,50
         om=(2*i+1)*Pi/Beta
         write(18,'(2f17.10)')om,imag(G0w(i))
      end do
      write(18,*) 
      write(18,'(a11)')'"Re(G0wand)"'
      do i=0,50
         om=(2*i+1)*Pi/Beta
         write(18,'(2f17.10)')om,real(G0wand(i))
      end do
      write(18,*) 
      write(18,'(a11)')'"Im(G0wand)"'
      do i=0,50
         om=(2*i+1)*Pi/Beta
         write(18,'(2f17.10)')om,imag(G0wand(i))
      end do
      write(18,*) 
C========+=========+=========+=========+=========+=========+=========+=$
c     END OF Iteration 
C========+=========+=========+=========+=========+=========+=========+=$
      write(6,'(a20)')'     End lisalanc '
      write(6,'(a60)')'========================================'
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: initial
C TYPE   : subroutine
C PURPOSE: initialize the search 
C I/O    :
C VERSION: 30-Sep-95
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine initial()
      include 'lisalanc.dat'
      rewind(15)
      call rheader(15)
cc
cc       Epsk(1) = - xmut - U/2
cc
      Epsk(1)=-xmu
      read(15,*)
      do i=2,Ns
         read(15,*)Epsk(i)
      end do
      read(15,*)
      do i=1,Ns-1
         read(15,*)V(i)
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: rheader.f  
C TYPE   : subroutine
C PURPOSE: read the program's header from unit k 
C I/O    :
C VERSION: 30-Sep-95
C COMMENT: 
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine rheader(k) 
      do 1 i=1,8
         read(k,*)
1     continue
      end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: wheader.f  
C TYPE   : subroutine
C PURPOSE: write the program's header onto unit k
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: 
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine wheader(k) 
      include 'lisalanc.dat'
      character *80 xyz
      character *8 xxx
      xxx=' 1-band '
      write(k,'(a55)')'========================================'
      write(k,'(a25,a30)')xxx,'30-Sep-95 LANCZOS  '
      write(k,'(a55)')'========================================'
      rewind(1)
      write(k,'(4(a6,I4))') 'NSITE ',Ns,'IWMAX',iwmax
      read(1,'(a60)')xyz
      read(1,'(a60)')xyz
      read(1,'(a60)')xyz
      do 3 i=1,4
      read(1,'(a60)')xyz
      write(k,'(a60)')xyz
3     continue 
      rewind(1)
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: diag.f  
C TYPE   : main 
C PURPOSE: diagonalization of Anderson problem
C          in standard (Landau - Lifschitz) notation by the Lanczos
C          Algorithm
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: To understand this subroutine, one may consult the correspon-
C          ding subroutine in the program lisadiag.f. The essential
C          points at which the present program differs from lisadiag.f
C          are indicated in the comments. To understand the method,
C          you may consult Gollub and van Loan (1983). 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine diag
      parameter (numzeromax=2)
      include 'lisalanc.dat'
      dimension vvinit(Nslmax),vectorout(Nslmax)
      dimension alfalanc(Nl), betalanc(Nl)
      dimension groundstate(Nslmax,numzeromax)
      dimension nupzero(numzeromax),ndozero(numzeromax)
      oldzero=1000
C 
C  Lanczos diagonalization 
C
      do nup=0,Ns
         do ndo=0,Ns
            call setup(nup,ndo)
cc
cc          let's compute the groundstate eigenvalue of the 
cc          present sector and compute the eigenvector
cc          (first application of the  Lanczos procedure).
cc
            dummy=nstates(nup,ndo)
            do i=1,Nstates(nup,ndo)
               vvinit(i)=One/sqrt(dummy)
            enddo
C
            call findgroundstatel (nstates(nup,ndo),vvinit,enemin,
     &            enemax,nlanc)
cc
cc          check the Lanczos computation of the eigenvector 
cc          by a simple vector iteration: The Groundstate eigen-
cc          vector should be invariant under the operation of 
cc          the Hamiltonian. Note that here we are very
cc          careful...we might in fact skip the following routine.
cc
            call findgroundstate (nstates(nup,ndo),vvinit,enemin,
     &            enemax)
            if (enemin.lt.oldzero-tiny*10.) then
               print*,enemin,'   enemin, new', numzero,'  numzero'
               numzero=1
               nupzero(1)=nup
               ndozero(1)=ndo
               oldzero=enemin
cc
cc             save groundstate vector
cc
               do i=1,nstates(nup,ndo)
                  groundstate(i,numzero)=vvinit(i)
               end do
cc
cc             we have found a new groundstate
cc 
            else if (abs(enemin-oldzero).le.tiny*10.) then
               print*,enemin,'   enemin, new',numzero,'  numzero'
               numzero=numzero+1
               if (numzero.gt.numzeromax) pause ' too many groundstates'
               nupzero(numzero)=nup
               ndozero(numzero)=ndo
               oldzero=min(oldzero,enemin)
cc
cc             save groundstate vector
cc
               do i=1,nstates(nup,ndo)
                  groundstate(i,numzero)=vvinit(i)
               end do
            end if
         end do
      end do
      print*,' here is the  list of groundstates:'
      do i=1,numzero
         print*,nupzero(i),ndozero(i)
      end do
cc
cc    for the  groundstates we've found...
cc
      do iw=0,iwmax
         Gw(iw) =cmplx(Zero,Zero)
      end do
cc
cc    it is evident that we have to do two calculations per
cc    groundstate, one in which we create a particle, and
cc    one in which a particle is destroyed (cf. Landau Lifshitz,
cc    Statistical Physics II, chapter 8). In order to simplify
cc    comprehension of the code, I put these two calculations
cc    one behind the other. Some computer time can be saved
cc    in regrouping the calculations.
cc    The loop over izero visits all the possible ground states
cc    of the system 
cc
      do izero=1, numzero
cc
cc       we first calculate the vector a+ |gs>. (which we use as new 
cc       vvinit
cc
         nup=nupzero(izero)
         ndo=ndozero(izero)
cc
cc       redo the setup of the Hilbert space of sector (nup, ndo)
cc       (we don't need the hamiltonian matrix)
cc
         call setup(nup,ndo)
         do i=1,nstates(nup+1,ndo)
            vvinit(i)=Zero
         end do
         vecdnorm=Zero
         do i=1,nstates(nup,ndo)
            if (Ivec(1,i).eq.0) then
               call adag(1,Nstock(i),k,isign1)
               k=Ninv(k)
               vvinit(k)=isign1*groundstate(i,izero)
               vecdnorm=vecdnorm+vvinit(k)**2
cc
cc             by keeping track of vecdnorm, we can compute the density
cc             <0|a^+ a|0>
            endif
         end do
cc
cc       now, we have produced a vector of the sector (nup+1,ndo), 
cc       which will serve as a initial vector of our Lanczos routine.
cc       First normalize the vector.
cc
         nup=nup+1
         do i=1,nstates(nup,ndo)
            vvinit(i)=vvinit(i)/sqrt(vecdnorm)
         end do
cc
cc       read in the relevant matrix elements
cc
         call setup(nup,ndo)
cc
cc       ... and redo the Lanczos iteration
cc
         do i=1,nstates(nup,ndo)
            vectorout(i)= Zero
         end do
         itermax=Nl
         nlanc=0
         do iter=1,itermax
            nlanc=nlanc+1
            call lanczos(vvinit,vectorout,alf, bet,iter,itermax,
     &        nstates(nup,ndo))
            alfalanc(iter)=alf
            if (iter.ne.Nl) betalanc(iter+1)=bet
            if (abs(bet).lt.tiny) goto 102
         end do
102      continue
cc 
cc       in this case, we are not at all interested in the groundstate
cc       (which could be calculated from the alfa's and beta's), but in the 
cc       values of the alfalanc's and betalanc's.
cc       From these two sets of parameters we can now calculate the 
cc       contribution to the Green's function, using a continued-
cc       fraction representation (cf. Lin and Gubernatis, 1993).
cc       In this example, we only compute thermodynamic Green's
cc       functions for the up-spin, but we could use exactly the same 
cc       method for calculating the real-frequency Green's functions
cc       (to get the density of states) etc.
cc
         isign=1
         call computegreen
     &     (vecdnorm,oldzero,nlanc,alfalanc,betalanc,isign)


cc
cc       redo the setup of the Hilbert space of sector (nup, ndo)
cc       (we don't need the hamiltonian matrix)
cc
         nup=nupzero(izero)
         ndo=ndozero(izero)
         call setup(nup,ndo)
cc
cc       we now calculate the vector [a |gs>] (which we use as new 
cc       vvinit)
cc
         do i=1,nstates(nup-1,ndo)
            vvinit(i)=Zero
         end do
         vecdnorm=Zero
         do i=1,nstates(nup,ndo)
            if (Ivec(1,i).eq.1) then
               call a(1,Nstock(i),k,isign1)
               k=Ninv(k)
               vvinit(k)=isign1*groundstate(i,izero)
               vecdnorm=vecdnorm+vvinit(k)**2
            endif
         end do
cc
cc       now, we have produced a vector of the sector (nup-1,ndo), 
cc       which will serve as a initial vector of our Lanczos routine.
cc       First normalize the vector.
cc
         nup=nup-1
         do i=1,nstates(nup,ndo)
            vvinit(i)=vvinit(i)/sqrt(vecdnorm)
         end do
cc
cc       read in the relevant matrix elements
cc
         call setup(nup,ndo)
cc
cc       ... and redo the Lanczos iteration
cc
         do i=1,nstates(nup,ndo)
            vectorout(i)= Zero
         end do
         itermax=Nl
         nlanc=0
         do iter=1,itermax
            nlanc=nlanc+1
            call lanczos(vvinit,vectorout,alf, bet,iter,itermax,
     &        nstates(nup,ndo))
            alfalanc(iter)=alf
            if (iter.ne.Nl) betalanc(iter+1)=bet
            if (abs(bet).lt.tiny) goto 10002
         end do
10002      continue
cc 
         isign=-1
         call computegreen
     &   (vecdnorm,oldzero,nlanc,alfalanc,betalanc,isign)
      end do
      factor=numzero
      do iomega=0,Iwmax
         gw(iomega)=gw(iomega)/factor
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: adaggger
C TYPE   : subroutine
C PURPOSE: operation of adagger (creation operator)
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: operation of adag(i) on state # j produces 
C          state # k with sign isign
C          k=0 === vacuum
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine adag(i,jold,k,isign)
      include 'lisalanc.dat'
      integer n(Ip),npov2(0:Ip),npovm1(0:Ip)
      logical logic
      data iswitch/1/
      save
      if (iswitch.eq.1) then
         iswitch=0
         npov2(0)=1
         npovm1(0)=1
         do it=1,Ip
            npov2(it)=2*npov2(it-1)
            npovm1(it)=-npovm1(it-1)
         end do
      end if 
cc
cc       calculate binary representation of number j
cc
      do ll=1,Ip
         logic=btest(jold,ll-1)
         n(Ip+1-ll)=0
         if(logic)n(Ip+1-ll)=1
      end do
cc
cc    calculate sign of new state
cc
      isign=0
      do 2 ll=1,i-1
         isign=isign+n(ll)
2     continue
      isign=npovm1(isign)
      k=jold+npov2(Ip-i)
      end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: a
C TYPE   : subroutine
C PURPOSE: operation of a  (destruction operator)
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: operation of a(i) on state # j produces 
C          state # k with sign isign
C          k=0 === vacuum
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine a(i,jold,k,isign)
      include 'lisalanc.dat'
      integer n(Ip)
      logical logic
cc
cc    calculate binary representation of number j
cc
      do ll=1,Ip
         logic=btest(jold,ll-1)
         n(Ip+1-ll)=0
         if(logic)n(Ip+1-ll)=1
      end do
c     print*,(n(kk),kk=1,6), 'a',jold
      if (n(i).eq.0) pause ' a'
cc
cc    calculate sign of new state
cc
      isign=0
      do 2 ll=1,i-1
         isign=isign+n(ll)
2     continue
      isign=(-1)**isign
      k=jold-2**(Ip-i)
      end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: buildbasis
C TYPE   : subroutine
C PURPOSE: construct the Hilbert space in the sector (nup,ndo)
C VERSION: 30-Sep-95
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine buildbasis(nup,ndo,Ns,nstates,Nstock,nsmax)
      implicit double precision(a-h,o-z)
      dimension Nstock(nsmax)
      parameter(nbmax=13000)
      integer*4 number(0:nbmax,0:16),point(0:16),i,nombre,poids
      integer*4 decal
      dimension nconfg(32)
      character*1 conf(32)
      logical permut
      npart=nup+ndo
      do i=1,nsmax
      Nstock(i)=0
      end do
      npartinit=npart
cc
cc generation of possible configurations
cc
      nbb=Ns
      permut=.false.
      if(2*npart.gt.nbb) then
         npart=nbb-npart
         permut=.true.
      endif
      poids=2**((nbb+1)/2)
      decal=(nbb+1)/2
      numbermax=2**nbb
      do i=0,16
         point(i)=0
      end do
      do i=0,poids-1
         call verif(i,nbit,nbb)
         point(nbit)=point(nbit)+1
C
         if(point(nbit).gt.nbmax)then
            write(*,*)'nbmax too small'
            stop
         endif
C
         number(point(nbit),nbit)=i
      end do
      indice=1
      kcp=0
      do i=0,npart
         j=npart-i
         do k=1,point(i)
            do l=1,point(j)
c           nombre=lshift(number(k,i),decal)+number(l,j)
            nombre=number(k,i)*2**decal+number(l,j)
*           nombre=poids*number(k,i)+number(l,j)
            if(nombre.lt.numbermax) then
               if(permut) nombre=not(nombre)
               call convert(nombre,conf,nbb,indice)
c              print*,' nconf'
c              write(*,192)(conf(jj),jj=Ns,1,-1)
c192            format(32a1)
               nsup=0
               nsdo=0
              do jj=1,Ns
                 if(conf(jj).eq.'0')then
                 nconfg(jj)=0
            else
             if(conf(jj).eq.'1')nconfg(jj)=1
                    if (jj.gt.Ns/2) then
                        nsup=nsup+1
            else
                        nsdo=nsdo+1
            endif
                end if
               end do
c              print*,nsup,nsdo,'   nsup, nsdo'
               if (nsup.eq.nup.and.nsdo.eq.ndo) then
                  kcp=kcp+1
                  do jj=1,Ns
                     if(nconfg(jj).eq.1)
     &                     Nstock(kcp)=ibset(Nstock(kcp),jj-1)
                  end do
c                 write(*,*)Nstock(kcp)
               end if
               indice=indice+1
            endif
            end do
         end do
      end do
      npart=npartinit
      nstates=kcp
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: convert
C TYPE   : subroutine
C VERSION: 30-Sep-95
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine convert(n,conf,dim,index)
      integer*4 n,dim,index
      character*1 conf(32)
      logical logic
      do i=0,dim-1
         logic=btest(n,i)
         if(logic) then
            conf(i+1)='1'
         else
            conf(i+1)='0'
         endif
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: verif
C TYPE   : subroutine
C VERSION: 30-Sep-95
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine verif(nombre,nbit,nboite)
      integer nombre,nbit,nboite
      logical logic
      nbit=0
      do i=0,15
         logic=btest(nombre,i)
         if(logic) then
            nbit=nbit+1
         endif
      end do
      return
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: search
C TYPE   : subroutine
C PURPOSE: 
C I/O    : parameters in lisalanc.dat 
C VERSION: 30-Sep-95
C COMMENT: programs kondo and moulin
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine search(fmin,Nitermax)
      include 'lisalanc.dat'
      dimension hess(nmpara**2),g(nmpara),xtemp(nmpara),
     &          w(nmpara**2),xprmt(nmpara)
      external energy
      data iexit/0/
      data iprint/0/
      data hh/1.e-5/
c number of parameters to be optimized (per spin):
c  eps(2)...... eps(Ns) --->  Ns-1
c  v(1)........ v(Ns-1) --->  Ns-1
      nbparm=2*Ns-2
      if (nbparm.gt.nmpara) stop ' nmpara too small'
C
         icount=0
         do i=2,Ns
            icount=icount+1
            xtemp(icount)=Epsk(icount+1)
         end do
         do i=1,Ns-1
            icount=icount+1
            xtemp(icount)=V(i)
         end do
         do i=nbparm+1,nmpara
            xtemp(i)=Zero
         end do
      
         do i=1,nbparm
            xprmt(i)=dabs(xtemp(i))+1.d-15
         end do
C
         mode=1
c use va10 for search
         dfn=-.5
         deps=.00001
         call minimize(energy,nbparm,xtemp,fmin,g,hess,w
     +   ,dfn,xprmt,hh,deps,mode,Nitermax,iprint,iexit)
         write (6,30) iexit,fmin
30       format(' iexit fmin ',i5,e14.6)
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: computegreen
C TYPE   : subroutine
C PURPOSE:
C I/O    : parameters in lisalanc.dat
C COMMENT: This subroutine calculates the CONTRIBUTION of the current
C          groundstate to the Green's function! Notice that 
C          Gw(iomega) may not be zero on entering this routine, since
C          several groundstates may be contributing to the zero-temperature
C          Green's function.
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine computegreen
     &          (anorm,enemin,nlanc,alfalanc,betalanc,isign)
      include 'lisalanc.dat'
      double complex ener,dl,det,omega
      dimension alfalanc(Nl),betalanc(Nl),dl(Nl),det(Nl)
cc
cc calculation of G(i omega) for imaginary omega
cc (the routine is trivially modified to calculate the Green's function
cc elsewhere).
c
      do iomega=0,Iwmax
         omega=(Two*iomega+One)*Pi/Beta*Xi
         ener=omega+enemin*isign
         do i=1,nlanc
            det(i)=cmplx(zero,zero)
         end do
         do i=1,nlanc
            dl(i)=ener-alfalanc(i)*isign
         enddo
         det(nlanc)=dl(nlanc)
         det(nlanc-1)=dl(nlanc-1)*dl(nlanc)-betalanc(nlanc)**2
         do i=nlanc-2,1,-1
            det(i)=dl(i)*det(i+1)-betalanc(i+1)**2*det(i+2)
         enddo
         gw(iomega)=gw(iomega)+anorm*det(2)/det(1)
      enddo
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: setup
C TYPE   : subroutine
C PURPOSE: compute the tables Hstate, Nconnect, Nhstate in the 
C          sector (nup,ndo)
C I/O    : 
C VERSION:
C COMMENT: Hstate, Nconnect, Nhstate are global variables, transmitted
C          in the common/general block
C          
C
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine setup(nup,ndo)
      include 'lisalanc.dat'
      dimension kcp(Nslmax)
      logical logic
cc
cc    The program starts as in the exact diagonalization case:
cc
cc    set up the Hilbert space: Nstock(i) and,
cc    equivalently, Ivec(*,i) code the  
cc    i'th vector of the present sector (nup,ndo). As the 
cc    only difference with the program lisadiag.f, we no
cc    longer stock the matrix for different sectors
cc
      call buildbasis(nup,ndo,Ip,
     &         Nstates(nup,ndo),Nstock,nslmax)
      if (Nstates(nup,ndo).gt.nslmax) stop ' nslmax too small' 
      do i=1,nstates(nup,ndo)
         Ninv(Nstock(i))=i
         do ll=1,Ip
            logic=btest(Nstock(i),ll-1)
            Ivec(Ip+1-ll,i)=0
            if(logic)Ivec(Ip+1-ll,i)=1
         end do
      end do
      do i=1,nstates(nup,ndo)
         kcp(i)=0
         Nhstate(i)=0
         do k=1,Ncontx
            Hstate(i,k)=0
         end do
      enddo
      do i=1,Nstates(nup,ndo)
cc
cc       diagonal part
cc
         xtemp= Epsk(1)*Ivec(1,i)
         xtemp= xtemp+Epsk(1)*Ivec(Ns+1,i)
         do k= 2,Ns
            xtemp=xtemp + Epsk(k)*Ivec(k,i)
     %        + Epsk(k)*Ivec(Ns+k,i)
         end do
         xtemp=xtemp + 
     %   U*Ivec(1,i)*(Ivec(Ns+1,i))
         kcp(i)=kcp(i)+1
         Nconnect(i,kcp(i))=i
         Hstate(i,kcp(i))=Hstate(i,kcp(i))+xtemp
         Nhstate(i)=Nhstate(i)+1
cc      
cc      non-diag - just to understand what's going on, realize that
cc      we do the a^+ and a operation on state i, to see that 
cc      it is connected to state k2. This means that we've found
cc      one more element of the row vector xmat(i,*) (i.e., we
cc      increment kcp(i)), but also (for reasons of symmetry of
cc      the hamiltonian), of the column vector xmat(*,k2). We
cc      therefore have to increment kcp(k2)
cc
         if (Ivec(1,i).eq.0) then
            do k=1,Ns-1
               if (Ivec(k+1,i).eq.1) then
                  call a(k+1,Nstock(i),k1,isign1)
                  call adag(1,k1,k2,isign2)
                  k2=Ninv(k2)
                  kcp(i)=kcp(i)+1
                  Nconnect(i,kcp(i))=k2
                  Hstate(i,kcp(i))=V(k)*isign1*isign2
                  Nhstate(i)=Nhstate(i)+1
                  kcp(k2)=kcp(k2)+1
                  Nconnect(k2,kcp(k2))=i
                  Hstate(k2,kcp(k2))=V(k)*isign1*isign2
                  Nhstate(k2)=Nhstate(k2)+1
               end if
            end do
         end if
         if (Ivec(Ns+1,i).eq.0) then
            do k=1,Ns-1
               if (Ivec(Ns+k+1,i).eq.1) then
                  call a(k+1+Ns,Nstock(i),k1,isign1)
                  call adag(1+Ns,k1,k2,isign2)
                  k2=Ninv(k2)
                  kcp(i)=kcp(i)+1
                  Nconnect(i,kcp(i))=k2
                  Hstate(i,kcp(i))=V(k)*isign1*isign2
                  Nhstate(i)=Nhstate(i)+1
                  kcp(k2)=kcp(k2)+1
                  Nconnect(k2,kcp(k2))=i
                  Hstate(k2,kcp(k2))=V(k)*isign1*isign2
                  Nhstate(k2)=Nhstate(k2)+1
               end if
            end do
         end if
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: minimize
C TYPE   : subroutine
C PURPOSE: conjugent gradient search
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: This is a most reliable conjugent gradient routine! It has
C          served us well for many years, and is capable to cope with
C          a very large number of variables. Unfortunately, we don't
C          know who wrote this routine (original name: 'va10a'), and 
C          we find it very obscure.
C          Don't worry, it works just fine.
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine minimize (funct, n, x, f, g, h, w, dfn, xm,
     $  hh, eps, mode, maxfn, iprint, iexit)
      implicit double precision (a-h,o-z)
      dimension x(*), g(*), h(*), w(*), xm(*)
      external funct
      data zero, half, one, two /0.0d0, 0.5d0, 1.0d0, 2.0d0/
      if (iprint .ne. 0) write (6,1000)
 1000 format (' entry into minimize')
      np = n + 1
      n1 = n - 1
      nn=(n*np)/2
      is = n
      iu = n
      iv = n + n
      ib = iv + n
      idiff = 1
      iexit = 0
      if (mode .eq. 3) go to 15
      if (mode .eq. 2) go to 10
      ij = nn + 1
      do 5 i = 1, n
      do 6 j = 1, i
      ij = ij - 1
   6  h(ij) = zero
   5  h(ij) = one
      go to 15
  10  continue
      ij = 1
      do 11 i = 2, n
      z = h(ij)
      if (z .le. zero) return
      ij = ij + 1
      i1 = ij
      do 11 j = i, n
      zz = h(ij)
      h(ij) = h(ij) / z
      jk = ij
      ik = i1
      do 12 k = i, j
      jk = jk + np - k
      h(jk) = h(jk) - h(ik) * zz
      ik = ik + 1
  12  continue
      ij = ij + 1
  11  continue
      if (h(ij) .le. zero) return
  15  continue
      ij = np
      dmin = h(1)
      do 16 i = 2, n
      if (h(ij) .ge. dmin) go to 16
      dmin = h(ij)
  16  ij = ij + np - i
      if (dmin .le. zero) return
      z = f
      itn = 0
      call funct (n, x, f)
      ifn = 1
      df = dfn
      if (dfn .eq. zero) df = f - z
      if (dfn .lt. zero) df = abs (df * f)
      if (df .le. zero) df = one
  17  continue
      do 19 i = 1, n
      w(i) = x(i)
  19  continue
      link = 1
      if (idiff - 1) 100, 100, 110
  18  continue
      if (ifn .ge. maxfn) go to 90
  20  continue
      if (iprint .eq. 0) go to 21
      if (mod (itn, iprint) .ne. 0) go to 21
       write (6,1001) itn, ifn
1001  format (1x,'itn = ',i5,' ifn = ',i5)
      write (6,1002) f
1002  format (1x,'f = ',e15.7)
      if (iprint .lt. 0) go to 21
      write (6,1003) (x(i), i = 1, n)
***
***
1003  format (1x,'x = ',4e15.7 / (5x, 4e15.7))
      write (6,1004) (g(i), i = 1, n)
1004  format (1x,'g = ',4e15.7 / (5x, 4e15.7))
**
***
  21  continue
      itn = itn + 1
      w(1) = -g(1)
      do 22 i = 2, n
      ij = i
      i1 = i - 1
      z = -g(i)
      do 23 j = 1, i1
      z = z - h(ij) * w(j)
      ij = ij + n - j
  23  continue
  22  w(i) = z
      w(is+n) = w(n) / h(nn)
      ij = nn
      do 25 i = 1, n1
      ij = ij - 1
      z = zero
      do 26 j = 1, i
      z = z + h(ij) * w(is+np-j)
      ij = ij - 1
  26  continue
  25  w(is+n-i) = w(n-i) / h(ij) - z
      z = zero
      gs0 = zero
      do 29 i = 1, n
      if (z * xm(i) .ge. abs (w(is+i))) go to 28
      z = abs (w(is+i)) / xm(i)
  28  gs0 = gs0 + g(i) * w(is+i)
  29  continue
      aeps = eps / z
      iexit = 2
      if (gs0 .ge. zero) go to 92
      alpha = -two * df / gs0
      if (alpha .gt. one) alpha = one
      ff = f
      tot = zero
      int = 0
      iexit = 1
  30  continue
      if (ifn .ge. maxfn) go to 90
      do 31 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  31  continue
      call funct (n, w, f1)
      ifn = ifn + 1
      if (f1 .ge. f) go to 40
      f2 = f
      tot = tot + alpha
  32  continue
      do 33 i = 1, n
      x(i) = w(i)
  33  continue
      f = f1
      if (int - 1) 35, 49, 50
  35  continue
      if (ifn .ge. maxfn) go to 90
      do 34 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  34  continue
      call funct (n, w, f1)
      ifn = ifn + 1
      if (f1 .ge. f) go to 50
      if ((f1 + f2 .ge. f + f) .and.
     $  (7.0d0 * f1 + 5.0d0 * f2 .gt. 12.0d0 * f)) int = 2
      tot = tot + alpha
      alpha = two * alpha
      go to 32
  40  continue
      if (alpha .lt. aeps) go to 92
      if (ifn .ge. maxfn) go to 90
      alpha = half * alpha
      do 41 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  41  continue
      call funct (n, w, f2)
      ifn = ifn + 1
      if (f2 .ge. f) go to 45
      tot = tot + alpha
      f = f2
      do 42 i = 1, n
      x(i) = w(i)
  42  continue
      go to 49
  45  continue
      z = 0.1d0
      if (f1 + f .gt. f2 + f2)
     $  z = one + half * (f - f1) / (f + f1 - f2 - f2)
      if (z .lt. 0.1d0) z = 0.1d0
      alpha = z * alpha
      int = 1
      go to 30
  49  continue
      if (tot .lt. aeps) go to 92
  50  continue
      alpha = tot
      do 56 i = 1, n
      w(i) = x(i)
      w(ib+i) = g(i)
  56  continue
      link = 2
      if (idiff - 1) 100, 100, 110
  54  continue
      if (ifn .ge. maxfn) go to 90
      gys = zero
      do 55 i = 1, n
      w(i) = w(ib+i)
      gys = gys + g(i) * w(is+i)
  55  continue
      df = ff - f
      dgs = gys - gs0
      if (dgs .le. zero) go to 20
      link = 1
      if (dgs + alpha * gs0 .gt. zero) go to 52
      do 51 i = 1, n
      w(iu + i) = g(i) - w(i)
  51  continue
      sig = one / (alpha * dgs)
      go to 70
  52  continue
      zz = alpha / (dgs - alpha * gs0)
      z = dgs * zz - one
      do 53 i = 1, n
      w(iu+i) = z * w(i) + g(i)
  53  continue
      sig = one / (zz * dgs * dgs)
      go to 70
  60  continue
      link = 2
      do 61 i = 1, n
      w(iu+i) = w(i)
  61  continue
      if (dgs + alpha * gs0 .gt. zero) go to 62
      sig = one / gs0
      go to 70
  62  continue
      sig = -zz
  70  continue
      w(iv+1) = w(iu+1)
      do 71 i = 2, n
      ij = i
      i1 = i - 1
      z = w(iu+i)
      do 72 j = 1, i1
      z = z - h(ij) * w(iv+j)
      ij = ij + n - j
  72  continue
      w(iv+i) = z
  71  continue
      ij = 1
      do 75 i = 1, n
      z = h(ij) + sig * w(iv+i) * w(iv+i)
      if (z .le. zero) z = dmin
      if (z .lt. dmin) dmin = z
      h(ij) = z
      w(ib+i) = w(iv+i) * sig / z
      sig = sig - w(ib+i) * w(ib+i) * z
      ij = ij + np - i
  75  continue
      ij = 1
      do 80 i = 1, n1
      ij = ij + 1
      i1 = i + 1
      do 80 j = i1, n
      w(iu+j) = w(iu+j) - h(ij) * w(iv+i)
      h(ij) = h(ij) + w(ib+i) * w(iu+j)
      ij = ij + 1
  80  continue
      go to (60, 20), link
  90  continue
      iexit = 3
      go to 94
  92  continue
      if (idiff .eq. 2) go to 94
      idiff = 2
      go to 17
  94  continue
      if (iprint .eq. 0) return
      write (6,1005) itn, ifn, iexit
1005  format (1x,'itn = ',i5, ' ifn = ',i5,' iexit = ',i5)
      write (6,1002) f
      write (6,1003) (x(i), i = 1, n)
      write (6,1004) (g(i), i = 1, n)
      return
 100  continue
      do 101 i = 1, n
      z = hh * xm(i)
      w(i) = w(i) + z
      call funct (n, w, f1)
      g(i) = (f1 - f) / z
      w(i) = w(i) - z
 101  continue
      ifn = ifn + n
      go to (18, 54), link
 110  continue
      do 111 i = 1, n
      z = hh * xm(i)
      w(i) = w(i) + z
      call funct (n, w, f1)
      w(i) = w(i) - z - z
      call funct (n, w, f2)
      g(i) = (f1 - f2) / (two * z)
      w(i) = w(i) + z
 111  continue
      ifn = ifn + n + n
      go to (18, 54), link
      end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: energy
C TYPE   : subroutine
C PURPOSE: calculates the energy 
C          (mismatch between G0_and and G_0new)
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine energy(nbparm,x,f)
        include 'lisalanc.dat'
        complex*16  cdummy1
        dimension x(nmpara)
        icount=0
        do i=2,Ns
           icount=icount+1
           Epsk(icount+1)=x(icount)
        end do
        do i=1,Ns-1
           icount=icount+1
           v(i)=x(icount)
        end do
        diff=Zero
        do i=0,Iwmax
           om=(Two*i+One)*Pi/Beta
           call calcg0(om,cdummy1)
           g0wand(i)=cdummy1
           diff = diff + abs(g0w(i)-g0wand(i))
        end do
        f=diff/dfloat(Iwmax+1)
        end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: calcg0
C TYPE   : subroutine
C PURPOSE: computes G0and F0and
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine calcg0(om,g0and)
        include 'lisalanc.dat'
        complex*16 g0and
cc
cc      use simple formula for the G_0 function
cc
        g0and=Xi*om-Epsk(1)-U/Two
        do i=1,Ns-1
           g0and=g0and-V(i)**2/(Xi*om-Epsk(i+1))
        end do
cc
cc      compiler bug on HP 735 calculates incorrectly 1/g0and
cc
c       g0and=One/g0and
        g0and=conjg(g0and)/(conjg(g0and)*g0and)
        end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: amultv.f
C TYPE   : subroutine
C PURPOSE: calculate  uu = uu + H vv
C          in the Lanczos procedure
C I/O    : parameters in lisalanc.dat
C COMMENT: Note that we suppose that the spare-matrix representation
C          of xmat has already been calculated (vector Hstate)
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine amultv(uu,vv,nsttot)
      include 'lisalanc.dat'
      dimension uu(Nslmax),vv(Nslmax)
      do i=1,nsttot
         do j=1,Nhstate(i)
            uu(i)=uu(i)+Hstate(i,j)*vv(Nconnect(i,j))
         end do
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: findgroundstatel
C TYPE   : subroutine
C PURPOSE: computes groundstate by Lanczos iteration +
C          QL diagonalization.
C I/O    :
C VERSION: 08/NOV/93
C COMMENT: if you are unhappy with this simple program, you 
C          may look into the SLATEC routines which provide
C          full-fledged Lanczos codes. However, these codes
C          are unnecessary, since we are only interested in
C          groundstate properties (cf. Gollub and van Loan (1983)
C          for a very thorough discussion).
C
C          on input, vvinit is an arbitrary start vector 
C          on output, vvinit is the groundstate vector, as computed by
C          the Lanczos algorithm.
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine findgroundstatel (nsttot,vvinit,enemin,enemax,nlanc)
      include 'lisalanc.dat'
      dimension z(Nl,Nl)
      dimension vvinit(Nslmax),vectorin(Nslmax),vectorout(Nslmax)
     *,alfalanc(Nl),betalanc(Nl),diag(Nl),subdiag(Nl)
cc
cc
cc    Here, we specify the starting vector of the iteration
cc
      do i=1,nsttot
         vectorin(i) = vvinit(i)
         vectorout(i)= Zero
      end do
cc
cc    Here, we call the proper Lanczos routine, i. e., we compute
cc    the coefficients alfalanc and betalanc.
cc
      itermax=min(Nl,nsttot)
      nlanc=0
      do iter=1,itermax
         nlanc=nlanc+1
         call lanczos(vectorin,vectorout,alf, bet,iter,itermax,nsttot)
         alfalanc(iter)=alf
         if (iter.ne.Nl) betalanc(iter+1)=bet
         if (abs(bet).lt.tiny) goto 1002
      end do
1002  continue
cc
cc    The little Lanczos routine results in a tridiagonal 
cc    (nlanc x nlanc) matrix
cc    (coded in the alfalanc and betalanc), which is simply
cc    diagonalized using the appropriate SLATEC routine tql2. 
cc    To check your
cc    particular application, you may want to verify that the 
cc    result depends VERY little on the arbitrary parameter Nl
cc    (as long as it is relatively large, cf the discussion of the 
cc    parameter Nl in appendix C of GKKR).
cc    
      do i=1,nlanc
         do j=1,nlanc
            Z(j,i)=Zero
         end do
         Z(i,i)=One
      end do
      do i=1,nlanc
         diag(i)=alfalanc(i)
      end do
      do i=2,nlanc
         subdiag(i)=betalanc(i)
      end do
cc
      call tql2(Nl,nlanc,diag,subdiag,Z,ierr)
cc
cc    eigenvalues are now in vector  diag. We keep the smallest
cc    one diag (1), and, for reference, the largest one d(Nl)
      enemin=diag(1)
      enemax=diag(nlanc)
cc
cc
cc    now we do the same thing  once over, in order to calculate the 
cc    groundstate eigenvector. Notice that, in the previous application
cc    of the subroutine lanczos, we did not save the vectors vectorin
cc    and vectorout (for lack of space). Here, we simply run once more
cc    through the subroutine lanczos, in order to 
cc    recreate the  different basis vectors (this procedure is standard).
cc    The vector vvinit now accumulates the contri-
cc    butions to the groundstate eigenvector.
cc     
      do i=1,nsttot
          vectorin(i)=vvinit(i)
          vectorout(i)=Zero
          vvinit(i)=Zero
      end do
      do iter=1,nlanc
         call lanczos(vectorin,vectorout,alf, bet,iter,itermax,nsttot)
         do i=1,nsttot
            vvinit(i)=vvinit(i)+vectorin(i)*Z(iter,1)
         end do
      end do
cc
cc    normalize the ground state vector
cc
      xnorm=Zero
      do i=1,nsttot
         xnorm=xnorm+vvinit(i)**2
      end do
      xnorm=sqrt(xnorm)
      do i=1,nsttot
         vvinit(i)=vvinit(i)/xnorm
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: lanczos
C TYPE   : subroutine
C PURPOSE: perform a no-frills  iteration of the Lanczos algorithm
C I/O    : parameters in lisalanc.dat
C COMMENT: We adopt a notation in which the (NxN) tridiagonal Lanczos matrix
C          is described by its diagonal elements alfalanc(1)...alfalanc(N)
C          and by the       subdiagonal elements beta(2)...beta(N)
C          the iteration i  determines alfalanc(i),betalanc(i+1)
C
C          Except for the first iteration (iter>1),
C          on input of iteration iter
C          vecold contains the basis vector w_(iter-1), while 
C          vecnew contains the basis vector w_iter*beta(iter) 
C          betalanc is the parameter beta(iter), and alfalanc is arbitrary
C
C          on output,   
C          vecold contains the basis vector w_(iter), while 
C          vecnew contains the basis vector w_(iter+1)*betalanc(iter+1) 
C          betalanc is the parameter betalanc(iter+1), and 
C          alfalanc is alfalanc(iter)
C
C          For iter=1
C          vecold contains the basis vector w_(1), the unit norm starting
C                 vector
C          vecnew contains the 0
C          betalanc and alfalanc can be arbitrary
C
C          if (iter.eq.itermax), the parameter betalanc(iter+1) is not
C          calculated
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine lanczos(vecold,vecnew,alfalanc,
     &         betalanc,iter,itermax,idim)
      include 'lisalanc.dat'
      dimension vecold(idim),vecnew(idim)
      idim1=idim
      if (iter.ne.1) then
cc
cc       transport vecnew into vecold,... use vecnew as a workvector,
cc       which now contains  -beta(i) w_i
cc
         do i=1,idim
            dummy=vecold(i)
            vecold(i)=vecnew(i)/betalanc
            vecnew(i)=-betalanc*dummy
         end do
      else
cc
cc       two basic checks which you can take out, if you're careful
cc
         xnorm=0.
         do i=1,idim
            if (abs(vecnew(i)).gt.tiny) pause 'error Lanczos 1'
            xnorm=xnorm+vecold(i)**2
         end do
         if (abs(xnorm-One).gt.tiny) then
            print*, ' error', xnorm, ' xnorm'
            print*,nup,ndo,' nup ndo;'
            pause 'error Lanczos 2'
         end if
      end if
      call amultv(vecnew,vecold,idim1)
c     print*,'after amultv vecold vecnew'
c     do i=1,idim
c        write(*,'(I4,2f15.4)')i,vecold(i),vecnew(i)
c     end do
      alfalanc=Zero
      do i=1,idim
         alfalanc=alfalanc+vecnew(i)*vecold(i)
      end do
c     print*,alfalanc, ' alfalanc'
      do i=1,idim
         vecnew(i)=vecnew(i)-alfalanc*vecold(i)
      end do
      if (iter.eq.itermax) goto 100
      betalanc=zero
      do i=1,idim
         betalanc=betalanc+vecnew(i)**2
      end do
      betalanc=sqrt(betalanc)
100   continue
      end  
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: findgroundstate
C TYPE   : subroutine
C PURPOSE: Determination of the groundstate by a simple projection
C          method. The subroutine serves as check, and as clean-up
C          of the lanczos procedure, in case the eigenvector 
C          has not very very well converged. In principle, this 
C          routine could be used instead of findgroundstatel. In 
C          practice, however, it is very often too slow in 
C          separating nearby eigenvalues and eigenvectors.
C I/O    : parameters in lisalanc.dat
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine findgroundstate(nsttot,gstate,enemin,enemax)
      include 'lisalanc.dat'
      dimension ctemp(Nslmax),gstate(Nslmax)
      relprec=1.e-10
      itermax=100
cc
cc    This procedure is extremely simple: 
cc    Matrix    Eigenvalues
cc    xmat      e_1< e_2 , ... < e_Nl   
cc
cc    since there is very fast convergence of both the extremal eigen
cc    values, we know that enemax, as calculated by Lanczos,
cc    is very close to the largest eigenvalue of the matrix. 
cc    The matrix -H + e_Nl   has therefore eigenvalues 
cc      -e_1 + enemax > ....... > -e_Nl + enemax = 0
cc    To avoid a problem if the dimension of the Hilbert space is 1,
cc    we calculate with enemaxx=enemax + 0.1
cc
      eneminold=enemin
      enemaxx=enemax+0.1
      do iter=1,100
cc
cc       Here we multiply (in fact) (-H + enemaxx) * vectorin
cc       we also calculate the Rayleigh quotient
cc       <gs|-H + enemaxx|gs>/<gs|gs> = sum2/sum1, which should be very 
cc       close to  -enemin+enemaxx
cc
         sum1=Zero
         do i=1,nsttot
            ctemp(i)=Zero
            sum1=sum1+gstate(i)**2
         end do
         call amultv(ctemp,gstate,nsttot)
cc
cc       now, we have ctemp = H * gstate
cc
         sum2=Zero
         sum3=Zero
         do i=1,nsttot
            gstatenew=gstate(i)*enemaxx -ctemp(i)
            sum2=sum2 + gstate(i)*gstatenew
            gstate(i)=gstatenew
            sum3=sum3+gstate(i)**2
         end do
         rayleigh=sum2/sum1
cc       calculation of the energy
cc
         eneminnew = enemaxx - rayleigh
         do i=1,nsttot
            gstate(i)=gstate(i)/sqrt(sum3)
         end do
         if (abs(eneminnew-eneminold).lt.relprec) goto 100
         eneminold=eneminnew
         if (mod(iter,50).eq.0)
     +        write(*,*)'iter= ',iter,' e0= ',esavenew
      end do
      print*,'niter too small to find  e0'
100   continue
cc
cc    get the normed eigenvector 
cc
      enemin=eneminnew
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: tql2
C TYPE   : subroutine
C PURPOSE: Compute the eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix
C I/O    : 
C VERSION:
C COMMENT: Minimally modified SLATEC routine from netlib.
C          To learn about netlib, send an otherwise empty
C          message to netlib@research.att.com
C          containing 'send index' in the subject header)
C          The WWW address of netlib is
C          http://netlib.att.com/netlib/search.html
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
*DECK TQL2
      SUBROUTINE TQL2 (NM, N, D, E, Z, IERR)
C***BEGIN PROLOGUE  TQL2
C***PURPOSE  Compute the eigenvalues and eigenvectors of symmetric
C            tridiagonal matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A5, D4C2A
C***TYPE      SINGLE PRECISION (TQL2-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TQL2,
C     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
C     Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     This subroutine finds the eigenvalues and eigenvectors
C     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
C     The eigenvectors of a FULL SYMMETRIC matrix can also
C     be found if  TRED2  has been used to reduce this
C     full matrix to tridiagonal form.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameter, Z, as declared in the calling program
C          dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        D contains the diagonal elements of the symmetric tridiagonal
C          matrix.  D is a one-dimensional DOUBLE PRECISION array, 
C          dimensioned D(N).
C
C        E contains the subdiagonal elements of the symmetric
C          tridiagonal matrix in its last N-1 positions.  E(1) is
C          arbitrary.  E is a one-dimensional DOUBLE PRECISION array, 
C          dimensioned E(N).
C
C        Z contains the transformation matrix produced in the
C          reduction by  TRED2, if performed.  If the eigenvectors
C          of the tridiagonal matrix are desired, Z must contain
C          the identity matrix.  Z is a two-dimensional DOUBLE PRECISION array,
C          dimensioned Z(NM,N).
C
C      On Output
C
C        D contains the eigenvalues in ascending order.  If an
C          error exit is made, the eigenvalues are correct but
C          unordered for indices 1, 2, ..., IERR-1.
C
C        E has been destroyed.
C
C        Z contains orthonormal eigenvectors of the symmetric
C          tridiagonal (or full) matrix.  If an error exit is made,
C          Z contains the eigenvectors associated with the stored
C          eigenvalues.
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C
C     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  PYTHAG
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  TQL2
C
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
      DOUBLE PRECISION D(*),E(*),Z(NM,*)
      DOUBLE PRECISION B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
      DOUBLE PRECISION PYTHAG
C
C***FIRST EXECUTABLE STATEMENT  TQL2
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F = 0.0D0
      B = 0.0D0
      E(N) = 0.0D0
C
      DO 240 L = 1, N
         J = 0
         H = ABS(D(L)) + ABS(E(L))
         IF (B .LT. H) B = H
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
         DO 110 M = L, N
            IF (B + ABS(E(M)) .EQ. B) GO TO 120
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
C
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
C     .......... FORM SHIFT ..........
         L1 = L + 1
         L2 = L1 + 1
         G = D(L)
         P = (D(L1) - G) / (2.0D0 * E(L))
         R = PYTHAG(P,1.0D0)
         D(L) = E(L) / (P + SIGN(R,P))
         D(L1) = E(L) * (P + SIGN(R,P))
         DL1 = D(L1)
         H = G - D(L)
         IF (L2 .GT. N) GO TO 145
C
         DO 140 I = L2, N
  140    D(I) = D(I) - H
C
  145    F = F + H
C     .......... QL TRANSFORMATION ..........
         P = D(M)
         C = 1.0D0
         C2 = C
         EL1 = E(L1)
         S = 0.0D0
         MML = M - L
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
         DO 200 II = 1, MML
            C3 = C2
            C2 = C
            S2 = S
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0D0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0D0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0D0)
            E(I+1) = S * E(I) * R
            S = 1.0D0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     .......... FORM VECTOR ..........
            DO 180 K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
C
  200    CONTINUE
C
         P = -S * S2 * C3 * EL1 * E(L) / DL1
         E(L) = S * P
         D(L) = C * P
         IF (B + ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
  240 CONTINUE
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
      END
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: pythag 
C TYPE   : function
C PURPOSE: compute dsqrt(a**2+b**2) without overflow or
C          destructive underflow
C I/O    : 
C VERSION:
C COMMENT: 
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      double precision function pythag(a,b)
      double precision a,b
      double precision p,r,s,t,u
      p = dmax1(dabs(a),dabs(b))
      if (p .eq. 0.0d0) go to 20
      r = (dmin1(dabs(a),dabs(b))/p)**2
   10 continue
         t = 4.0d0 + r
         if (t .eq. 4.0d0) go to 20
         s = r/t
         u = 1.0d0 + 2.0d0*s
         p = u*p
         r = (s/u)**2 * r
      go to 10
   20 pythag = p
      return
      end
