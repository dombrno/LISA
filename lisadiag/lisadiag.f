C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: lisadiag.f  
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
C           in: Reviews of Modern Physics 68, 13 (1996)
C
C          (the paper will be referred to as ``GKKR'').
C          you are kindly asked to cite the paper (and, if applicable,
C              the original works) if you use this program.
C
C          The present algorithm was first described in:
C          M. Caffarel and W. Krauth, Phys. Rev. Let. 72, 1545 (1994)
C
C          The programs have  been thoroughly tested on SUN Sparks,
C          HP 900, IBM RS6000 stations.
C
C TYPE   : main 
C PURPOSE: self-consistency using Anderson diag.
C I/O    : cf file README_lisadiag
C VERSION: 30-Sep-95
C VERSION: 07-Mar-03 improved compatibility with gnu f77 compiler
C AUTHOR : W. Krauth (krauth@physique.ens.fr). Parts of the code were 
C          originally written by M. Caffarel.
C COMMENT: It is impossible to present a code which contains all
C          the sophistications of the versions in actual use.
C          In this version, we only calculate thermodynamic Green's
C          functions for the single-band Hubbard model, in the 
C          paramagnetic phase. Generalization of the program to more
C          complicated cases was actually performed. The necessary 
C          modifications are straightforward.
C          Even though FORTRAN is case-insensitive, I have capitalized
C          all the global variables, i. e. the variables appearing
C          in the common/global/ block (cf file lisadiag.dat).
C          At the end of the program, a number of SLATEC routines
C          have been appended. You may not want to  print these.

C========+=========+=========+=========+=========+=========+=========+=$
      program lisadiag
      include 'lisadiag.dat' 
      character*80 xyz 
      complex*16 omega, adummy, cdummy1
      Zero=0
      One=1
      Two=2
      Pi=acos(-One)
      Xi=cmplx(Zero,One)
      Half=.true.
cc
cc    open parameter file
cc
      open (unit=1,file='lisadiag.input',form='formatted',status='old')
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
      open (15,file='lisadiag.andpar',form='formatted',status='old')
      read (1,'(a80)') xyz 
      read (1,*) itaucalc,jcut
      Taucalc=.false.
      if (itaucalc.eq.1) Taucalc=.true.
      open (unit=17,file='lisadiag.green',form='formatted',
     &  status='new')
      open (unit=18,file='lisadiag.diff',form='formatted',
     &  status='new')
cc    open andpar file
      call initial()
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
      call search(chi2,Nitermax)
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
C PLOT (for example time-dependent Green's functions, which were
C computed). 
C========+=========+=========+=========+=========+=========+=========+=$
      write(17,'(a18,I4,3f7.2)')' TitleText: G(w)  ',
     %Ns,Beta,U,xmu
      do i=0,50
         om=(2*i+1)*Pi/Beta
         write(17,'(2f17.10)')om,dble(Gw(i))
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
         write(18,'(2f17.10)')om,dble(G0w(i))
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
         write(18,'(2f17.10)')om,dble(G0wand(i))
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
      write(6,'(a20)')'         Summary - lisadiag'
      write(6,'(a60)')'========================================'
      write(*,'(a20,f15.7)')' DENSITY         =',Two*Ddens
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: initial
C TYPE   : subroutine
C PURPOSE: initialize the search 
C I/O    :
C VERSION: 30-Sep-95
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine initial()
      include 'lisadiag.dat'
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
Cnoprint===============================================================
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
Cnoprint===============================================================
      subroutine wheader(k) 
      include 'lisadiag.dat'
      character *80 xyz
      character *8 xxx
      xxx=' 1-band '
      write(k,'(a55)')'========================================'
      write(k,'(a25,a30)')xxx,'30-Sep-95 FULL DIAG'
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
C          in standard (Landau - Lifschitz) notation
C          put onto external memory 
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine diag
      include 'lisadiag.dat'
      dimension xmat(nsmax,nsmax)
      dimension evec(nsmax,nsmax),evecp(nsmax,nsmax)
      dimension nstock(nsmax,0:Ns,0:Ns),ninv(0:2**ip-1)
      dimension e0(0:Ns,0:Ns),fv1(nsmax),fv2(nsmax)
      dimension dup(nsmax,nsmax),ddo(nsmax,nsmax)
      character*10 filename
      complex*16 fac3,omega
      logical logic
      e0ref=10.
      betstep=Beta/L
C 
C  full diagonalization by QR algorithm
C
      Xnorm=Zero
      id=1
      if (id.eq.1) then
          Beta1=Beta
      else
         Beta1=10000
         print*,'ATTENTION BETA1.ne.BETA know what you_re doing?',Beta1
      end if
      do nup=0,Ns
         do ndo=0,Ns
            call buildbasis(nup,ndo,ip,
     &               nstates(nup,ndo),nstock(1,nup,ndo),nsmax)
            if (nstates(nup,ndo).gt.nsmax) stop ' nsmax too small' 
c           print*,nstates(nup,ndo),' <=== number of states'
            do i=1,nstates(nup,ndo)
               ninv(nstock(i,nup,ndo))=i
               do ll=1,ip
                  logic=btest(nstock(i,nup,ndo),ll-1)
                  Ivec(ip+1-ll,i,nup,ndo)=0
                  if(logic)Ivec(ip+1-ll,i,nup,ndo)=1
               enddo
c              write(*,'(3i4,10i1)')i,nstock(i,nup,ndo),
c    &            (Ivec(ii,i,nup,ndo),ii=1,ip)
            end do
         end do
      end do
      do nup = 0,Ns
         do ndo = 0,Ns 
            do i=1,nstates(nup,ndo)
               do j=1,nstates(nup,ndo)
                  xmat(i,j)=Zero
               end do
            end do
            do i=1,nstates(nup,ndo)
cc
cc             diagonal part
cc
               xmat(i,i)= Epsk(1)*Ivec(1,i,nup,ndo)
               xmat(i,i)= xmat(i,i)+Epsk(1)*Ivec(Ns+1,i,nup,ndo)
               do k= 2,Ns
                  xmat(i,i)=xmat(i,i) + Epsk(k)*Ivec(k,i,nup,ndo)
     %              + Epsk(k)*Ivec(Ns+k,i,nup,ndo)
               end do
               xmat(i,i)=xmat(i,i) + 
     %         U*Ivec(1,i,nup,ndo)*(Ivec(Ns+1,i,nup,ndo))
cc      
cc      non-diag 
cc
               if (Ivec(1,i,nup,ndo).eq.0) then
                  do k=1,Ns-1
                     if (Ivec(k+1,i,nup,ndo).eq.1) then
                        call a(k+1,nstock(i,nup,ndo),k1,isign1)
                        call adag(1,k1,k2,isign2)
                        k2=ninv(k2)
                        xmat(i,k2)=V(k)*isign1*isign2
                        xmat(k2,i)=xmat(i,k2)
                     end if
                  end do
               end if
               if (Ivec(Ns+1,i,nup,ndo).eq.0) then
                  do k=1,Ns-1
                     if (Ivec(Ns+k+1,i,nup,ndo).eq.1) then
                        call a(k+1+Ns,nstock(i,nup,ndo),k1,isign1)
                        call adag(1+Ns,k1,k2,isign2)
                        k2=ninv(k2)
                        xmat(i,k2)=V(k)*isign1*isign2
                        xmat(k2,i)=xmat(i,k2)
                     end if
                  end do
               end if
            end do
            matz=1
            call  rs (nsmax,nstates(nup,ndo),xmat,Eval(1,nup,ndo), 
     &      matz, evec, fv1, fv2, ierr)
c           print*,' ground state vector'
c           do i=1,nstates(nup,ndo)
c              write(*,'(f20.10)')evec(i,1)
c           end do
            if (ierr.gt.0) print*, ' problem in RS', ierr
cc
cc          write eigenvectors onto file
cc 
            write(filename,'(a2,2I1,a6)')'zz',nup,ndo,'evfile'
            open(99,FORM='UNFORMATTED',file=filename,status='unknown')
            write(99)((evec(j,i),j=1,nstates(nup,ndo)),
     &      i=1,nstates(nup,ndo))
            close(99)
            emint=1.e20
            do i=1,nstates(nup,ndo)
               if(Eval(i,nup,ndo).lt.emint)emint=Eval(i,nup,ndo)
            end do
            e0(nup,ndo)=emint
            print*,'e0(',nup,ndo,')= ',e0(nup,ndo)
         end do
      end do
cc
cc
cc
      if(jcut.eq.1)then
      emint=1.e20
      do nup=0,Ns
         do ndo=0,Ns
            do i=1,nstates(nup,ndo)
              if(Eval(i,nup,ndo).lt.emint)emint=Eval(i,nup,ndo)
            enddo
         enddo
      enddo
c     print*,'emint= ',emint
      endif
cc
cc
cc
      do nup=0,Ns
         do ndo=0,Ns
            do i=1,nstates(nup,ndo)
cc
cc          renormalization of energies
cc
               if(jcut.eq.1)Eval(i,nup,ndo)=Eval(i,nup,ndo)-emint
c              write(21,'(f20.10)')Eval(i,nup,ndo)
               Xnorm=Xnorm+exp(-Beta1*Eval(i,nup,ndo))
            end do 
         end do
      end do
cc
cc    compute G(tau)
cc
      if (Taucalc) then
         iu=0
         do tau=0,Beta-.01,betstep
            iu=iu+1
            Gt(iu)=Zero
         end do
         do nup=0,Ns-1
            do ndo=0,Ns
cc
cc          read in files
cc
               write(filename,'(a2,2I1,a6)')'zz',nup,ndo,'evfile'
               open(99,FORM='UNFORMATTED',file=filename,status='old')
               read(99)((evec(j,i),j=1,nstates(nup,ndo)),
     &         i=1,nstates(nup,ndo))
               close(99)
               write(filename,'(a2,2I1,a6)')'zz',nup+1,ndo,'evfile'
               open(99,FORM='UNFORMATTED',file=filename,status='old')
               read(99)((evecp(j,i),j=1,nstates(nup+1,ndo)),
     &         i=1,nstates(nup+1,ndo))
               close(99)

               do i=1,nstates(nup+1,ndo)
                  do j=1,nstates(nup,ndo)
cc
cc      dupij is the matrix element <evec(i) | adag_up | evec(j)>
cc
                     dup(i,j)=Zero
                     do ll=1,nstates(nup,ndo)
                        if (Ivec(1,ll,nup,ndo).eq.0) then
                           call adag(1,nstock(ll,nup,ndo),k,isign1)
                           k=ninv(k)
                           dup(i,j)=dup(i,j)+
     %                       isign1*evec(ll,j)*evecp(k,i)
                        end if
                     end do
                     iu=0
                     do tau=0,Beta-.01,betstep
                        iu=iu+1
                        gt(iu)=gt(iu)-exp(-tau*Eval(i,nup+1,ndo))*
     %                     exp(-(Beta1-tau)*Eval(j,nup,ndo))*dup(i,j)**2
                     end do
                  end do
               end do
cc
cc          dump dup(i,j) onto file (green's function also calculated)
cc
               write(filename,'(a2,2I1,a6)')'zz',nup,ndo,'dupfil'
               open(99,file=filename,FORM='UNFORMATTED',
     &                      status='unknown')
               write(99)((dup(i,j),i=1,nstates(nup+1,ndo)),
     &         j=1,nstates(nup,ndo))
               close(99) 
            end do 
         end do 
         do nup=0,Ns
            do ndo=0,Ns-1
cc
cc          read in files - calculation of ddo(i,j)
cc
               write(filename,'(a2,2I1,a6)')'zz',nup,ndo,'evfile'
               open(99,FORM='UNFORMATTED',file=filename,status='old')
               read(99)((evec(j,i),j=1,nstates(nup,ndo)),
     &         i=1,nstates(nup,ndo))
               close(99)
               write(filename,'(a2,2I1,a6)')'zz',nup,ndo+1,'evfile'
               open(99,FORM='UNFORMATTED',file=filename,status='old')
               read(99)((evecp(j,i),j=1,nstates(nup,ndo+1)),
     &         i=1,nstates(nup,ndo+1))
               close(99)
               do i=1,nstates(nup,ndo+1)
                  do j=1,nstates(nup,ndo)
cc
cc      ddoij is the matrix element <evec(i) | adag_do | evec(j)>
cc
                     ddo(i,j)=Zero
                     do ll=1,nstates(nup,ndo)
                        if (Ivec(1+Ns,ll,nup,ndo).eq.0) then
                           call adag(1+Ns,nstock(ll,nup,ndo),k,isign1)
                           k=ninv(k)
                           ddo(i,j)=ddo(i,j)+
     %                       isign1*evec(ll,j)*evecp(k,i)
                        end if
                     end do
                  end do
               end do
cc
cc             dump ddo(i,j) onto file
cc
               write(filename,'(a2,2I1,a6)')'zz',nup,ndo,'ddofil'
               open(99,file=filename,FORM='UNFORMATTED',
     &                 status='unknown')
               write(99)((ddo(i,j),i=1,nstates(nup,ndo+1)),
     &         j=1,nstates(nup,ndo))
               close(99) 
            end do 
         end do 
         iu=0
         do tau=0,Beta-.01,betstep
            iu=iu+1
            Gt(iu)=Gt(iu)/Xnorm
         end do
         Ddens=1+Gt(1)
      end if
cc
cc    compute local spin-spin correlation and its integral
cc
      if (locspin) call chichi
cc    compute Gw(omega)
cc
      do i=0,iwmax
         Gw(i)=Zero
cc
cc       add the following line if you are calculating Greal
c        Greal(i)=0
cc
      end do
      do nup=0,Ns-1
         do ndo=0,Ns
cc
cc          read in dup(i,j)
cc
            call fetchup(nup,ndo,dup)

            do i=1,nstates(nup+1,ndo)
               do j=1,nstates(nup,ndo)

                  term=exp(-Beta1*Eval(i,nup+1,ndo)) + 
     &                   exp(-Beta1*Eval(j,nup,ndo))
                  if(jcut.eq.1)then
                     if(term.lt.Cutoff)go to 111
                  endif
                  fac=term/Xnorm
                  fac2=(Eval(j,nup,ndo) -Eval(i,nup+1,ndo))**2
cc
cc                Imaginary frequency Green's function
cc
                  do iw=0,iwmax
                     omega=(Two*iw+1)*Pi/Beta*Xi
                     fac3=Eval(j,nup,ndo) -Eval(i,nup+1,ndo) - omega
                     if (Half) fac3= - omega
cc
cc                   symmetrize Gw for half-filling in h=0
cc
                     gw(iw)=gw(iw) + 
     &               dup(i,j)**2*fac/(-omega**2+fac2)*fac3
                  end do
cc
cc                It is very simple to compute Real frequecy Green's 
CC                functions
cc                Here is how this is done (commented out)
cc
cc                do iw=0,iwmax
cc                   omega= -omemax + Two*omemax/real(iwmax)*iw
cc   &                 + Xi * epsilon
cc                   fac3=Eval(j,nup,ndo) -Eval(i,nup+1,ndo) - omega
cc                   Greal(iw)=Greal(iw) + 
cc   &                 dup(i,j)**2*fac/(-omega**2+fac2)*fac3
cc                end do
111               continue
               end do
            end do
         end do
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: adaggger
C TYPE   : subroutine
C PURPOSE: operation of adagger
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: operation of adag(i) on state # j produces 
C          state # k with sign isign
C          k=0 === vacuum
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine adag(i,jold,k,isign)
      include 'lisadiag.dat'
      integer n(ip),npov2(0:ip),npovm1(0:ip)
      logical logic
      data iswitch/1/
      save
      if (iswitch.eq.1) then
         iswitch=0
         npov2(0)=1
         npovm1(0)=1
         do it=1,ip
            npov2(it)=2*npov2(it-1)
            npovm1(it)=-npovm1(it-1)
         end do
      end if 
cc
cc       calculate binary representation of number j
      do ll=1,ip
         logic=btest(jold,ll-1)
         n(ip+1-ll)=0
         if(logic)n(ip+1-ll)=1
      enddo
cc
cc    calculate sign of new state
cc
      isign=0
      do 2 ll=1,i-1
         isign=isign+n(ll)
2     continue
      isign=npovm1(isign)
      k=jold+npov2(ip-i)
      end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: a
C TYPE   : subroutine
C PURPOSE: operation of c
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: operation of a(i) on state # j produces 
C          state # k with sign isign
C          k=0 === vacuum
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine a(i,jold,k,isign)
      include 'lisadiag.dat'
      integer n(ip)
      logical logic
cc
cc    calculate binary representation of number j
      do ll=1,ip
         logic=btest(jold,ll-1)
         n(ip+1-ll)=0
         if(logic)n(ip+1-ll)=1
      enddo
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
      k=jold-2**(ip-i)
      end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: chichi
C TYPE   : subroutine
C PURPOSE: calculate local spin properties
C I/O    : 
C VERSION: 30-Sep-95
C COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine chichi
      include 'lisadiag.dat'
      dimension evec(nsmax,nsmax)
      character*10 filename
      betstep=Beta/L
      iu=0
      do tau=Zero,Beta-.01,betstep
         iu=iu+1
         Chis(iu)=Zero
      enddo
      chilocd=Zero
      do i=0,iwmax
         Chiw(i)=Zero
      end do

      do nup=0,Ns
         do ndo=0,Ns
cc
cc       read in files
cc
            write(filename,'(a2,2I1,a6)')'zz',nup,ndo,'evfile'
            open(99,FORM='UNFORMATTED',file=filename,status='old')
            read(99)((evec(j,i),j=1,nstates(nup,ndo)),
     &      i=1,nstates(nup,ndo))
            do i=1,nstates(nup,ndo)
               do j=1,nstates(nup,ndo)
                  term=exp(-Beta1*Eval(i,nup,ndo))-
     &                 exp(-Beta1*Eval(j,nup,ndo))
                  spinij=Zero
                  do k=1,nstates(nup,ndo)
                     spinij=spinij+evec(k,i)*evec(k,j)*
     &               (dfloat(Ivec(1,k,nup,ndo)+
     &               Ivec(1+Ns,k,nup,ndo))-One)
                  enddo
                  fac=term/Xnorm
                  fac2=(Eval(j,nup,ndo)-Eval(i,nup,ndo))
                  if(abs(Beta1*fac2).lt.1.e-10)then
                     chilocd=chilocd+
     &               spinij**2*exp(-Beta1*Eval(i,nup,ndo))*
     &                       Beta1/Xnorm
                  else   
                     chilocd=chilocd+spinij**2*fac/fac2
                  endif

cc    Compute spin-spin correlation function in imaginary frequency
cc    iwn is bosonic for spin-spin!!!!(2 point corr. function)
cc
                  do iw=1,200
                     omega=real(iw)/200.
                     fac3=Eval(j,nup,ndo)-Eval(i,nup,ndo)+Xi*omega
                     Chiw(iw)=Chiw(iw)+
     &                 spinij**2*fac/(fac2**2+omega**2)*fac3
                  end do

cc   Treat iw=0 separately
cc
                  if(abs(fac2).lt.1.e-12)then
                     Chiw(0)=Chiw(0)+
     &               spinij**2*exp(-Beta1*Eval(i,nup,ndo))*Beta1/Xnorm
                  else
                     Chiw(0)=Chiw(0)+spinij**2*fac/fac2
                  end if

                  iu=0
                  do tau=0,Beta1-.01,betstep
                     iu=iu+1
                     Chis(iu)=Chis(iu)+exp(-tau*Eval(i,nup,ndo))*
     &               exp(-(Beta1-tau)*Eval(j,nup,ndo))*spinij**2
                  enddo        
               end do
            end do
         end do
      end do
      iu=0
      do tau=0,Beta1-.01,betstep
         iu=iu+1
         Chis(iu)=Chis(iu)/Xnorm
      enddo
cc
cc    end of local chi calculation
cc
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: fetchup
C TYPE   : subroutine
C PURPOSE: read in the ddagger matrix for a given set of 
C          quantum numbers
C VERSION: 30-Sep-95
C========+=========+=========+=========+=========+=========+=========+=$
       subroutine fetchup(nup,ndo,darray)
       include 'lisadiag.dat'
       dimension darray(nsmax,nsmax)
       character*10 filename
       write(filename,'(a2,2I1,a6)')'zz',nup,ndo,'dupfil'
       open(99,FORM='UNFORMATTED',file=filename,status='old')
       read(99)((darray(i,j),i=1,nstates(nup+1,ndo)),
     & j=1,nstates(nup,ndo))
       close(99)
       end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: fetchdo
C TYPE   : subroutine
C PURPOSE: read in the ddagger matrix for a given set of 
C          quantum numbers
C VERSION: 30-Sep-95
C========+=========+=========+=========+=========+=========+=========+=$
       subroutine fetchdo(nup,ndo,darray)
       include 'lisadiag.dat'
       dimension darray(nsmax,nsmax)
       character*10 filename
       write(filename,'(a2,2I1,a6)')'zz',nup,ndo,'ddofil'
       open(99,FORM='UNFORMATTED',file=filename,status='old')
       read(99)((darray(i,j),i=1,nstates(nup,ndo+1)),
     & j=1,nstates(nup,ndo))
       close(99)
       end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: buildbasis
C TYPE   : subroutine
C PURPOSE: construct the Hilbert space in the sector (nup,ndo)
C VERSION: 30-Sep-95
Cnoprint================================================================
      subroutine buildbasis(nup,ndo,Ns,nstates,nstock,nsmax)
      implicit double precision(a-h,o-z)
      dimension nstock(nsmax)
      parameter(nbmax=13000)
      integer*4 number(0:nbmax,0:16),point(0:16),i,nombre,poids
      integer*4 decal
      dimension nconfg(32)
      character*1 conf(32)
      logical permut
      npart=nup+ndo
      do i=1,nsmax
      nstock(i)=0
      enddo
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
      enddo
      do i=0,poids-1
         call verif(i,nbit,nbb)
         point(nbit)=point(nbit)+1

         if(point(nbit).gt.nbmax)then
            write(*,*)'nbmax too small'
            stop
         endif

         number(point(nbit),nbit)=i
      enddo
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
     &                     nstock(kcp)=ibset(nstock(kcp),jj-1)
                  enddo
c                 write(*,*)nstock(kcp)
               end if
               indice=indice+1
            endif
            enddo
         enddo
      enddo
      npart=npartinit
      nstates=kcp
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: convert
C TYPE   : subroutine
C VERSION: 30-Sep-95
Cnoprint================================================================
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
      enddo
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: verif
C TYPE   : subroutine
C VERSION: 30-Sep-95
Cnoprint================================================================
      subroutine verif(nombre,nbit,nboite)
      integer nombre,nbit,nboite
      logical logic
      nbit=0
      do i=0,15
         logic=btest(nombre,i)
         if(logic) then
            nbit=nbit+1
         endif
      enddo
      return
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: search
C TYPE   : subroutine
C PURPOSE: 
C I/O    : parameters in lisadiag.dat 
C VERSION: 30-Sep-95
C COMMENT: programs kondo and moulin
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine search(fmin,Nitermax)
      include 'lisadiag.dat'
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

         icount=0
         do i=2,Ns
            icount=icount+1
            xtemp(icount)=Epsk(icount+1)
         enddo
         do i=1,Ns-1
            icount=icount+1
            xtemp(icount)=V(i)
         enddo
         do i=nbparm+1,nmpara
            xtemp(i)=Zero
         end do
      
         do i=1,nbparm
            xprmt(i)=dabs(xtemp(i))+1.d-15
         enddo

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
Cnoprint================================================================
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
        include 'lisadiag.dat'
        complex*16  cdummy1
        dimension x(nmpara)
        icount=0
        do i=2,Ns
           icount=icount+1
           Epsk(icount+1)=x(icount)
        enddo
        do i=1,Ns-1
           icount=icount+1
           v(i)=x(icount)
        enddo
        diff=Zero
        do i=0,Iwmax
           om=(Two*i+One)*Pi/Beta
           call calcg0(om,cdummy1)
           g0wand(i)=cdummy1
           diff = diff + abs(g0w(i)-g0wand(i))
        enddo
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
        include 'lisadiag.dat'

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
C PROGRAM: a few subroutines from slatec
C          which were minimally modified to compute 
C          double precision eigensystems.
C TYPE   : main
C PURPOSE: 
C I/O    :
C VERSION: 
C COMMENT: to learn about netlib, send an otherwise empty
C          message to netlib@research.att.com
C          containing 'send index' in the subject header)
C          The WWW address of netlib is
C          http://netlib.att.com/netlib/search.html
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
*DECK RS
      SUBROUTINE RS (NM, N, A, W, MATZ, Z, FV1, FV2, IERR)
C***BEGIN PROLOGUE  RS
C***PURPOSE  Compute the eigenvalues and, optionally, the eigenvectors
C            of a real symmetric matrix.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4A1
C***TYPE      SINGLE PRECISION (RS-S, CH-C)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine calls the recommended sequence of
C     subroutines from the eigensystem subroutine package (EISPACK)
C     to find the eigenvalues and eigenvectors (if desired)
C     of a DOUBLE PRECISION SYMMETRIC matrix.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        A contains the real symmetric matrix.  A is a two-dimensional
C          DOUBLE PRECISION array, dimensioned A(NM,N).
C
C        MATZ is an INTEGER variable set equal to zero if only
C          eigenvalues are desired.  Otherwise, it is set to any
C          non-zero integer for both eigenvalues and eigenvectors.
C
C     On Output
C
C        A is unaltered.
C
C        W contains the eigenvalues in ascending order.  W is a one-
C          dimensional DOUBLE PRECISION array, dimensioned W(N).
C
C        Z contains the eigenvectors if MATZ is not zero.  The
C          eigenvectors are orthonormal.  Z is a two-dimensional
C          DOUBLE PRECISION array, dimensioned Z(NM,N).
C
C        IERR is an INTEGER flag set to
C          Zero       for normal return,
C          10*N       if N is greater than NM,
C          J          if the J-th eigenvalue has not been
C                     determined after 30 iterations.
C                     The eigenvalues, and eigenvectors if requested,
C                     should be correct for indices 1, 2, ..., IERR-1.
C
C        FV1 and FV2 are one-dimensional DOUBLE PRECISION arrays used for temporary
C          storage, dimensioned FV1(N) and FV2(N).
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  TQL2, TQLRAT,  TRED2
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  RS
C
      INTEGER N,NM,IERR,MATZ
      DOUBLE PRECISION A(NM,*),W(*),Z(NM,*),FV1(*),FV2(*)
C
C***FIRST EXECUTABLE STATEMENT  RS
      IERR = 10 * N
C     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
      CALL  TRED2(NM,N,A,W,FV1,Z)
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      END
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
C          matrix.  D is a one-dimensional DOUBLE PRECISION array, dimensioned D(N).
C
C        E contains the subdiagonal elements of the symmetric
C          tridiagonal matrix in its last N-1 positions.  E(1) is
C          arbitrary.  E is a one-dimensional DOUBLE PRECISION array, dimensioned
C          E(N).
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
*DECK TRED2
      SUBROUTINE TRED2 (NM, N, A, D, E, Z)
C***BEGIN PROLOGUE  TRED2
C***PURPOSE  Reduce a real symmetric matrix to a symmetric tridiagonal
C            matrix using and accumulating orthogonal transformations.
C***LIBRARY   SLATEC (EISPACK)
C***CATEGORY  D4C1B1
C***TYPE      SINGLE PRECISION (TRED2-S)
C***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
C***AUTHOR  Smith, B. T., et al.
C***DESCRIPTION
C
C     This subroutine is a translation of the ALGOL procedure TRED2,
C     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     This subroutine reduces a DOUBLE PRECISION SYMMETRIC matrix to a
C     symmetric tridiagonal matrix using and accumulating
C     orthogonal similarity transformations.
C
C     On Input
C
C        NM must be set to the row dimension of the two-dimensional
C          array parameters, A and Z, as declared in the calling
C          program dimension statement.  NM is an INTEGER variable.
C
C        N is the order of the matrix A.  N is an INTEGER variable.
C          N must be less than or equal to NM.
C
C        A contains the real symmetric input matrix.  Only the lower
C          triangle of the matrix need be supplied.  A is a two-
C          dimensional DOUBLE PRECISION array, dimensioned A(NM,N).
C
C     On Output
C
C        D contains the diagonal elements of the symmetric tridiagonal
C          matrix.  D is a one-dimensional DOUBLE PRECISION array, dimensioned D(N).
C
C        E contains the subdiagonal elements of the symmetric
C          tridiagonal matrix in its last N-1 positions.  E(1) is set
C          to zero.  E is a one-dimensional DOUBLE PRECISION array, dimensioned
C          E(N).
C
C        Z contains the orthogonal transformation matrix produced in
C          the reduction.  Z is a two-dimensional DOUBLE PRECISION array,
C          dimensioned Z(NM,N).
C
C        A and Z may coincide.  If distinct, A is unaltered.
C
C     Questions and comments should be directed to B. S. Garbow,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C     ------------------------------------------------------------------
C
C***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
C                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
C                 system Routines - EISPACK Guide, Springer-Verlag,
C                 1976.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   760101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  TRED2
C
      INTEGER I,J,K,L,N,II,NM,JP1
      DOUBLE PRECISION A(NM,*),D(*),E(*),Z(NM,*)
      DOUBLE PRECISION F,G,H,HH,SCALE
C
C***FIRST EXECUTABLE STATEMENT  TRED2
      DO 100 I = 1, N
C
         DO 100 J = 1, I
            Z(I,J) = A(I,J)
  100 CONTINUE
C
      IF (N .EQ. 1) GO TO 320
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
      DO 300 II = 2, N
         I = N + 2 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 2) GO TO 130
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(Z(I,K))
C
         IF (SCALE .NE. 0.0D0) GO TO 140
  130    E(I) = Z(I,L)
         GO TO 290
C
  140    DO 150 K = 1, L
            Z(I,K) = Z(I,K) / SCALE
            H = H + Z(I,K) * Z(I,K)
  150    CONTINUE
C
         F = Z(I,L)
         G = -SIGN(SQRT(H),F)
         E(I) = SCALE * G
         H = H - F * G
         Z(I,L) = F - G
         F = 0.0D0
C
         DO 240 J = 1, L
            Z(J,I) = Z(I,J) / H
            G = 0.0D0
C     .......... FORM ELEMENT OF A*U ..........
            DO 180 K = 1, J
  180       G = G + Z(J,K) * Z(I,K)
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
  200       G = G + Z(K,J) * Z(I,K)
C     .......... FORM ELEMENT OF P ..........
  220       E(J) = G / H
            F = F + E(J) * Z(I,J)
  240    CONTINUE
C
         HH = F / (H + H)
C     .......... FORM REDUCED A ..........
         DO 260 J = 1, L
            F = Z(I,J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 260 K = 1, J
               Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
  260    CONTINUE
C
  290    D(I) = H
  300 CONTINUE
C
  320 D(1) = 0.0D0
      E(1) = 0.0D0
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
      DO 500 I = 1, N
         L = I - 1
         IF (D(I) .EQ. 0.0D0) GO TO 380
C
         DO 360 J = 1, L
            G = 0.0D0
C
            DO 340 K = 1, L
  340       G = G + Z(I,K) * Z(K,J)
C
            DO 360 K = 1, L
               Z(K,J) = Z(K,J) - G * Z(K,I)
  360    CONTINUE
C
  380    D(I) = Z(I,I)
         Z(I,I) = 1.0D0
         IF (L .LT. 1) GO TO 500
C
         DO 400 J = 1, L
            Z(I,J) = 0.0D0
            Z(J,I) = 0.0D0
  400    CONTINUE
C
  500 CONTINUE
C
      RETURN
      END

      double precision function pythag(a,b)
      double precision a,b
c
c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c
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
