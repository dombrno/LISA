C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: qmcexample.f  
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
C          Unfortunately, we cannot guarantee correctness of the programs,
C          but they have been thoroughly tested on SUN Sparks, HP 9000,
C          IBM RS6000 stations.
C
C TYPE   : main 
C PURPOSE: An Anderson impurity problem (cf eq. (\ref{HAM}) with a total of 
C          N_s
C          orbitals is considered (1 impurity orbital , N_s-1 conduction
C          bath orbitals). 
C          The imaginary time interval Beta is devided into L time-slices. 
C
C          This PURELY PEDAGOGICAL example program should be used
C          for 'reasonable' values of L < 8 - 16, and N_s < 3 - 5
C          The purpose of this program is 
C          a) to illustrate the Section \ref{method} of  GKKR
C          b) to provide a test example for the QMC code lisaqmc.f
C
C          The program is made up of several sections:
C
C          Section 1
C          a) The Hamiltonian H_0 of the model is constructed
C             (cf. eq. (\ref{HAM}). 
C          b) exp(Tau*H_0) and exp(+/-Xlambda*V) are constructed
C             (a simple algorithm is used exp(B) = 1 + B - B^2/2 + ...)
C          c) A sequence of Ising spins is chosen: s(1),...,s(L)
C          d) The matrices B_1,...,B_L are constructed 
C
C          Section 2
C          application of the BSS algorithm (cf eqs (\ref{detsmall}) and
C          eq. (\ref{Gdiscret}). 
C          a) The matrices  B_1, B_1*B_2 , B_1*B_2* ... * B_L 
C             are computed. The spectrum of their eigenvalues are
C             put into file BSS_SPECTRUM. 
C             The determinant of the matrix B is calculated as the product 
C             of eigenvalues
C             (note that there are faster ways for computing the
C             determinant, but that it's precisely the eigenvalue spectrum 
C             which
C             governs the possible numerical instabilities). 
C          b) the Green's function G_{s1,...,sL} is computed 
C             (cf eq. (\ref{Gdiscret})
C
C
C          Section 3
C          a) The matrix O, as defined in eq. \ref{bigO}) is computed
C          b) the spectrum of O is computed and put onto file
C             O_SPECTRUM.
C          c) the determinant of the matrix O is calculated
C          c) The Green's function  G=O^{-1} is computed by explicit
C             matrix inversion.
C
C          Section 4
C          a) the Green's function G_0 = G(s1 = s2 = ... = sL = 0) is computed
C          b) the matrix A, as defined in eq (\ref{inversion}) is 
C             computed.
C          c) the Green's function G_{s1,...,sL} is computed according
C             to (cf eq. (\ref{}) (G_{s1,...,sL} = A^{-1} G_0 
C
C
C VERSION: 17-Sep-96 (three minor bugs removed).
C AUTHOR : W. Krauth (krauth@physique.ens.fr)
C COMMENT: 
C NOTAT. : Even though FORTRAN is case insensitive, I've exclusively
C          capitalized the variables also appearing in GKKR.
C========+=========+=========+=========+=========+=========+=========+=$
      implicit double precision (a-h,o-z)
      parameter (N_s=4,L=4,Iwmax=2**13)
      dimension Epsk(N_s), Vk(2:N_s),
     &B_up(N_s,N_s,L), B_do(N_s,N_s,L),
     &Bfin_up(N_s,N_s), Bfin_do(N_s,N_s),
     &H_0(N_s,N_s),Is(L),wr(L*N_s),wi(L*N_s),
     &O_up(L*N_s,L*N_s), O_do(L*N_s,L*N_s),
     &O_upinv(L*N_s,L*N_s), O_doinv(L*N_s,L*N_s),
     &Expmv(N_s,N_s), Exppv(N_s,N_s),
     &Green_up(L,L),Green_do(L,L),Exp_mH_0(N_s,N_s),
     &G0_up(L,L),G0_do(L,L),
     &dummymat_up(N_s,N_s),dummymat_do(N_s,N_s),
     &prodmat(N_s,N_s),
     &Green0_up(L),gtemp_up(-L+1:L-1)
      double complex G0w(0:iwmax), xi,det_up,det_do
      real dran
      open(1,file='BSS_SPECTRUM',status='new',access='sequential')
      open(2,file='O_SPECTRUM',status='new',access='sequential')
      idum=-3491811
      pi=acos(-1.)
      xi =(0.,1.)
      zero=0
      one=1
      two=2
      half=one/two
cc
cc    physical parameters (you may choose your own)
cc
      Beta=4
      U=4
      Deltau=Beta/L
cc
cc    the following two lines determine Xlambda such 
cc    that cosh(Xlambda)= exp(Deltau*U/2.)   (U>0)
cc
      dummy=exp(Deltau*U/two)
      Xlambda=log(dummy+sqrt(dummy**2-one))
cc
cc    Section 1 a)
cc 
      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print*,'Section 1'
      print*
      call null(H_0,N_s)
      call null(Exppv,N_s)
      call null(Expmv,N_s)
      do i=1,N_s
         epsk(i)=dran(idum)-half
         H_0(i,i)=epsk(i)
      end do
      do i=2,N_s
         Vk(i)=dran(idum)-half
         H_0(1,i)=Vk(i)
         H_0(i,1)=H_0(1,i)
      end do

cc
cc    Section 1 b): compute exp(-Deltau*H_0)
cc                  and expmv(1,1)= exp(-Xlambda V); as well as  
cc                  exppv(1,1)= exp(Xlambda V)
cc
      print*,'Ho',Deltau
      call matrixform(H_0,N_s)
      call exponentiate(H_0,-Deltau,Exp_mH_0,N_s)
      Exppv(1,1)=exp(Xlambda)  
      Expmv(1,1)=exp(-Xlambda)
      do i=2,N_s
         Exppv(i,i)=one
         Expmv(i,i)=one
      end do      
      print*,'Matrix Exp(v) (cf eq. (\ref{ztrott}))'
      call matrixform(Exppv,N_s)
cc
cc    Section 1 c): choose sequence of Ising spins
cc                  (you may choose your own)
cc
      do i=1,L
         Is(i)=1
         if (dran(idum).gt.0.5) Is(i)=-1
      end do
      print*,' Choice of Ising spins:'
      write(*,'(20I2)')(Is(j),j=1,L)
cc
cc    Section 1 d): construct the matrices B
cc                  notice that they depend on the physical spin
cc                  and are called B_up and B_do
cc                  B_i (=B(k,l,i);k=1...N_s;l=1,...,N_s)
cc                  Notice that there are only two different
cc                  matrices involved
cc
      do islice=1,L
         if (Is(islice).eq.1) then
            call matmult(Exp_mH_0,ExppV,B_up(1,1,islice),N_s)
            call matmult(Exp_mH_0,ExpmV,B_do(1,1,islice),N_s)
         else
            call matmult(Exp_mH_0,ExpmV,B_up(1,1,islice),N_s)
            call matmult(Exp_mH_0,ExppV,B_do(1,1,islice),N_s)
         end if
         print*,' B_up', islice
         call matrixform(B_up(1,1,islice),N_s)
      end do
cc
cc    Section 2 a) calculate B_1*B_2*...*B_L
cc
      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print*,'Section 2'
      print*
      call null(Bfin_up,N_s)
      call null(Bfin_do,N_s)
      do i=1,N_s
         Bfin_up(i,i)=one
         Bfin_do(i,i)=one
      end do
      do iprod=1,L
         call matmult(Bfin_up,B_up(1,1,iprod),dummymat_up,N_s)
         call matmult(Bfin_do,B_do(1,1,iprod),dummymat_do,N_s)
         call transfer (dummymat_up,Bfin_up,N_s)
         call transfer (dummymat_do,Bfin_do,N_s)
cc
cc   ...it's not totally trivial to calculate the eigenvalue spectrum
cc      of a non-diagonal matrix 
cc
      end do
      call eigenvalues(Bfin_up,wr,wi,ierr,N_s)
      det_up=one
      do i=1,N_s
         det_up=det_up*(wr(i)+one)
         write(1,*) wr(i),i 
         if (abs(wi(i)).gt.1.e-10) print*, 
     &      'check for imaginary part of eigenvalue Bfin_up'
      end do
      print*,det_up,'det_up'
      call eigenvalues(Bfin_do,wr,wi,ierr,N_s)
      det_do=one
      do i=1,N_s
         det_do=det_do*(wr(i)+one)
         write(1,*) wr(i),i 
         if (abs(wi(i)).gt.1.e-10) print*, 
     &      'check for imaginary part of eigenvalue Bfin_do'
      end do
      print*,det_do,'det_do'
cc
cc    Section 2 b) The Green's function is computed according to the
cc                 BSS algorithm. We only calculate the up-spin
cc                 Green's functions for positive times
cc
      call null(Green_up,L)
cc
cc    compute the product B_l2 ... B_1 B_L ... B_{l2+1}
cc
      do i = 1, L
         do j=1,i
cc
cc          one independent calculation per element of G
cc
            call null(dummymat_up,N_s)
            do ii=1,N_s
               dummymat_up(ii,ii)=one
            end do
            do k=j,L
               call matmult (B_up(1,1,k),dummymat_up,prodmat,N_s)
               call transfer(prodmat,dummymat_up,N_s)
            end do
            do k=1,j-1
               call matmult (B_up(1,1,k),dummymat_up,prodmat,N_s)
               call transfer(prodmat,dummymat_up,N_s)
            end do
            do ii=1,N_s
               prodmat(ii,ii)=prodmat(ii,ii)+one
            end do
            call inverse(prodmat,N_s,N_s,dummymat_up)
cc
cc          do the numerator
cc
            do k=j,i-1
               call matmult (B_up(1,1,k),dummymat_up,prodmat,N_s)
               call transfer(prodmat,dummymat_up,N_s)
            end do
cc
cc          this settles the computation of the formula (\ref{Gdiscret})
cc          the d-Green's function is given by the 1-1 element
cc          of the final matrix dummymat_up
cc
            Green_up(i,j)=dummymat_up(1,1)
         end do
      end do
      print*,'Green_up(i,j) acc to the BSS algorithm (for (i.ge.j))'
      call matrixform(Green_up,L)
cc
cc    Section 3 a) Notice that there are two matrices
cc    O_up and O_do for convenience only (One could also use
cc    a single matrix of size      (2*L*N_s) x (2*L*N_s)).
cc
      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print*,'Section 3'
      print*
      call null (O_up,L*N_s)
      call null (O_do,L*N_s)
      do i=1,L*N_s
         O_up(i,i)=one
         O_do(i,i)=one
      end do
cc
cc    here, the matrix O is filled as indicated in eq. (\ref{bigO})
cc
      do i=1,L-1
         do j=1,N_s
            do k=1,N_s
               ioff=(i-1)*N_s
               O_up(ioff+j+N_s,ioff+k)=-B_up(j,k,i)
               O_do(ioff+j+N_s,ioff+k)=-B_do(j,k,i)
            end do
         end do
      end do
cc
cc    wrap-around for block O(1,l)
cc
      do j=1,N_s
         do k=1,N_s
            ioff=(L-1)*N_s
            O_up(j,ioff+k)=B_up(j,k,L)
            O_do(j,ioff+k)=B_do(j,k,L)
         end do
      end do
cc
cc    compute eigenvalues and determinant of O_up and O_do
cc    (again : there are much faster ways to compute det O)
cc
      write(2,'(a40)')'Eigenvalues O matrix'
      det_up=one
      ieigenflag=0
      call eigenvalues(O_up,wr,wi,ierr,L*N_s)
      do i=1,L*N_s
         write(2,*) wr(i),i 
         if (abs(wi(i)).gt.1.e-10.and.ieigenflag.eq.0) then
            print*, 
     &      'check for imaginary part of eigenvalue'
            ieigenflag=1
         end if
         det_up=det_up*(wr(i)+xi*wi(i))
      end do
      write(2,*)
      det_do=one
      call eigenvalues(O_do,wr,wi,ierr,L*N_s)
      do i=1,L*N_s
         write(2,*) wr(i),i 
         if (abs(wi(i)).gt.1.e-10.and.ieigenflag.eq.0) then
            print*, 
     &      'check for imaginary part of eigenvalue'
            ieigenflag=1
         end if
         if (abs(wi(i)).gt.1.e-10) print*, 
     &      'check for imaginary part of eigenvalue'
         det_do=det_do*(wr(i)+xi*wi(i))
      end do
      print*,real(det_up),real(det_do),' detup,detdo'
cc
cc    Section 3 c): compute the Green's function by an explicit
cc                  inversion of O (notice that this is done only 
cc                  for purposes of illustration). 
cc
      call inverse(O_up,L*N_s,L*N_s,O_upinv)
      call inverse(O_do,L*N_s,L*N_s,O_doinv)
cc
cc    the d-d Greens function is part of O^(-1)
cc
      do i=1,L
         do j=1,L
            Green_up(i,j)=O_upinv((i-1)*N_s+1,(j-1)*N_s+1)
            Green_do(i,j)=O_doinv((i-1)*N_s+1,(j-1)*N_s+1)
         end do
      end do
      print*,'Green_up according O^{-1}'
      call matrixform(Green_up,L)
cc
cc    
cc
cc
cc    Section 4 a): compute the Green's function by an explicit
cc                  calculation of G0(omega) + fourier transform
cc
      print*,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      print*,'Section 4'
      print*
      do i=0,iwmax
         omega=pi*(2*i+one)/Beta
         G0w(i)=xi*omega-epsk(1)
         do j=2,N_s
            G0w(i)=G0w(i)-Vk(j)**2/(xi*omega-epsk(j))
         end do
         G0w(i)=-one/G0w(i)
      end do
      call invfourier(G0w,Green0_up,iwmax,L,Beta)
       print*,'Green0_up'
       print*,Green0_up
cc
cc    notice that Green0_up is just a vector, it has now to be 
cc    made into a matrix the following few lines are taken from the 
cc    subroutine initial  of the production-scale MC program,
cc    notice that gtemp(0) is the limit of times ---> 0 ^+
cc
      do i=0,L-1
         gtemp_up(i)=Green0_up(i+1)
      end do
cc
cc    reflection of G to calculate Greens function for neg. arguments
cc
      do  i=1,L-1
         gtemp_up(-i)=-gtemp_up(L-i)
      end do
      do i=1,L
         do j=1,L
            G0_up(i,j)=gtemp_up(i-j)
cc
cc          in this case (no external fields, no symmetry breaking,
cc          the G0 is independent of spin.
cc
            G0_do(i,j)=G0_up(i,j)
         end do
      end do
 
      print*,'G0_up from DYSON eq'
      call matrixform(G0_up,L)
cc
cc    the subroutine update determines A, A^(-1), and the product
cc    A^(-1)*Green0, the subroutine is almost identical to the one
cc    used in the production-size program... 
cc
      call update(G0_up,G0_do,Green_up,Green_do,Is,L,Xlambda)
      print*,'Green_up from DYSON eq'
      call matrixform(Green_up,L)
 
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: update.f
C TYPE   : subroutine
C PURPOSE:
C          calculate the Green's function
C          for a given configuration of spins
C          (in vector Is) from the Green's function
C          for spins set equal to zero.
C I/O    :
C VERSION: 30-Sep-95
C COMMENT: 
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine update(G0up,G0do,Gup,Gdo,Is,L,Xlambda)
        implicit double precision (a-h,o-z)
        parameter (nm=20)
        dimension a(nm,nm),b(nm,nm),ainv(nm,nm),binv(nm,nm)
        dimension Is(L),Gup(L,L),Gdo(L,L)
        dimension G0up(L,L),G0do(L,L)
        if (L.gt.nm) pause 'problem in update'
        one=1
        zero=0
cc
cc      calculate the matrix a=1-(g-1)(exp(v'-1))  NB: v=0
cc
        do i=1,L
           do j=1,L
              a(i,j)=-G0up(i,j)*(exp(Xlambda*real(Is(j)))-one)
              b(i,j)=-G0do(i,j)*(exp(-Xlambda*real(Is(j)))-one)
           end do
           a(i,i)=1-(G0up(i,i)-one)*(exp(Xlambda*real(Is(i)))-one)
           b(i,i)=1-(G0do(i,i)-one)*(exp(-Xlambda*real(Is(i)))-one)
        end do
        call inverse(a,L,nm,ainv)
        call inverse(b,L,nm,binv)
        do i=1,L
           do j=1,L
              suma=zero
              sumb=zero
              do k=1,L
                 suma=suma+ainv(i,k)*G0up(k,j)
                 sumb=sumb+binv(i,k)*G0do(k,j)
              end do
              Gup(i,j)=suma
              Gdo(i,j)=sumb
           end do
        end do
        end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: invfourier
C TYPE   : subroutine
C PURPOSE: inverse fourier transform
C          Greent, Greenw use physical definition
C          Greent(i)=G((i-1)*Deltau) for i=1,...,T
C          Greenw(n)=G(i w_n), for n=0,T/2-1
C                w_n=(2*n+1)pi/Beta
C          Symmetry property:
C          G(iw_(-n)=G(iw_(n-1))*
C          coupled to the impurity
C I/O    :
C VERSION: 30-Sep-95
C COMMENT:
C========+=========+=========+=========+=========+=========+=========+=$
       subroutine invfourier(G0w,G0,iwmax,L,Beta)
       implicit double precision(a-h,o-z)
       dimension G0(L)
       complex*16 G0w(0:Iwmax),cdummy
       xpi=acos(-1.)
       do i=1,L
          G0(i)=0.
          tau=(i-1)*Beta/L
          do j=0,Iwmax
             om=mod((2*j+1.)*xpi/Beta*tau,2*xpi)
             cdummy=(0.,1.)*om
             dummy=G0w(j)*exp(-cdummy)
             G0(i)=G0(i)+2./Beta*dummy
          end do
       end do
cc
cc     special treatment for tau=0
cc
       G0(1)= 1./2.+G0(1)
       end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: exponentiate
C TYPE   : subroutine
C PURPOSE: laboriously exponentiate a matrix 
C          (exp(x A) = 1. + x A +(x A)**2/2 ... )
C          use only in this example program
C I/O    :
C VERSION: 30-Sep-95
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine exponentiate(A,x,exp_xA,N_s)
        implicit double precision (a-h,o-z)
        parameter (nsmax=5)
        dimension A(N_s,N_s),exp_xA(N_s,N_s),
     &  dummy(Nsmax,Nsmax),prod_A(Nsmax,Nsmax)
        if (nsmax.lt.N_s) pause 'nsmax too small '
        call null(exp_xA,N_s)
        call null(prod_A,Nsmax)

        do i=1,N_s
           exp_xA(i,i)=1
           prod_A(i,i)=1
        end do
        fac=1
        do iter=1,20
           fac=fac*iter
cc
cc         calculate matrix product  dummy = prod_A * A
cc
           do i=1,N_s
              do j=1,N_s
                 dummy(i,j)=0
                 do k=1,N_s
                    dummy(i,j)=dummy(i,j)+prod_A(i,k)*A(k,j)
                 end do
              end do
           end do
           do i=1,N_s
              do j=1,N_s
                 prod_A(i,j)=dummy(i,j)*x
                 exp_xA(i,j)=exp_xA(i,j) + prod_A(i,j)/fac
              end do
           end do
        end do
        end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: matrixform
C TYPE   : subroutine
C PURPOSE: print a square matrix in matrix form
C VERSION: 30-Sep-95
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine matrixform(a,N_s)
        implicit double precision(a-h,o-z)
        dimension a(N_s,N_s)
        do i=1,N_s
           write(*,'(10f8.4)')(a(i,j),j=1,N_s)
        end do
        end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: dran
C TYPE   : real function
C PURPOSE: generator of (dirty) randum numbers. Use only
C          in this example program
C I/O    :
C VERSION: 30-Sep-95
C========+=========+=========+=========+=========+=========+=========+=$
        real function dran(idum)
        parameter(ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836)
        k=idum/iq 
        idum=ia*(idum - k*iq) -ir*k
        if (idum.lt.0) idum=idum+im
        dran=am*idum
        end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: matmult
C TYPE   : subroutine
C PURPOSE: multiply two square matrices
C I/O    :
C VERSION: 30-Sep-95
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine matmult(a,b,c,N)
      implicit double precision(a-h,o-z)
      dimension a(n,n),b(n,n),c(n,n)
      do i=1,n
         do j=1,n
            c(i,j)=0
            do k=1,n
               c(i,j)=c(i,j)+a(i,k)*b(k,j)
            end do
         end do
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: eigenvalues
C TYPE   : subroutine
C PURPOSE: calculate eigenvalues of general matrix
C          the routines were obtained from netlib@research.att.com
C         (to learn about netlib, send an otherwise empty
C          message to netlib@research.att.com 
C          containing 'send index' in the subject header)
C          The WWW address of netlib is
C          http://netlib.att.com/netlib/search.html
C I/O    :
C VERSION: 30-Sep-95
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine eigenvalues(xmat,wr,wi,ierr,N)
      implicit double precision (a-h,o-z)
      parameter (Nmax=300)
      dimension xmat(N,N),wr(N),wi(N)
      dimension fv1(Nmax)
      dimension a(Nmax,Nmax)
      integer iv1(Nmax)
      if (Nmax.lt.N) pause 'size of matrix too large'
      do i=1,N
         do j=1,N
            a(i,j)=xmat(i,j)
         end do
      end do
      call  balanc(nmax,n,a,is1,is2,fv1)
      call  elmhes(nmax,n,is1,is2,a,iv1)
      call  hqr(nmax,n,is1,is2,a,wr,wi,ierr)
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: balanc
C TYPE   : subroutine
C PURPOSE: balances a real matrix and isolates eigenvalues 
C          whenever possible
C I/O    :
C VERSION: 30-Sep-95
C COMMENT: obtained from netlib@research.att.com
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine balanc(nm,n,a,low,igh,scale)
      implicit double precision (a-h,o-z)
c
      integer i,j,k,l,m,n,jj,nm,igh,low,iexc
      dimension a(nm,n),scale(n)
      logical noconv
c
c     this subroutine is a translation of the algol procedure balance,
c     num. math. 13, 293-304(1969) by parlett and reinsch.
c     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).
c
c     this subroutine balances a real matrix and isolates
c     eigenvalues whenever possible.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the input matrix to be balanced.
c
c     on output
c
c        a contains the balanced matrix.
c
c        low and igh are two integers such that a(i,j)
c          is equal to zero if
c           (1) i is greater than j and
c           (2) j=1,...,low-1 or i=igh+1,...,n.
c
c        scale contains information determining the
c           permutations and scaling factors used.
c
c     suppose that the principal submatrix in rows low through igh
c     has been balanced, that p(j) denotes the index interchanged
c     with j during the permutation step, and that the elements
c     of the diagonal matrix used are denoted by d(i,j).  then
c        scale(j) = p(j),    for j = 1,...,low-1
c                 = d(j,j),      j = low,...,igh
c                 = p(j)         j = igh+1,...,n.
c     the order in which the interchanges are made is n to igh+1,
c     then 1 to low-1.
c
c     note that 1 is returned for igh if igh is zero formally.
c
c     the algol procedure exc contained in balance appears in
c     balanc  in line.  (note that the algol roles of identifiers
c     k,l have been reversed.)
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      radix = 16.0
c
      b2 = radix * radix
      k = 1
      l = n
      go to 100
c     .......... in-line procedure for row and
c                column exchange ..........
   20 scale(m) = j
      if (j .eq. m) go to 50
c
      do 30 i = 1, l
         f = a(i,j)
         a(i,j) = a(i,m)
         a(i,m) = f
   30 continue
c
      do 40 i = k, n
         f = a(j,i)
         a(j,i) = a(m,i)
         a(m,i) = f
   40 continue
c
   50 go to (80,130), iexc
c     .......... search for rows isolating an eigenvalue
c                and push them down ..........
   80 if (l .eq. 1) go to 280
      l = l - 1
c     .......... for j=l step -1 until 1 do -- ..........
  100 do 120 jj = 1, l
         j = l + 1 - jj
c
         do 110 i = 1, l
            if (i .eq. j) go to 110
            if (a(j,i) .ne. 0.) go to 120
  110    continue
c
         m = l
         iexc = 1
         go to 20
  120 continue
c
      go to 140
c     .......... search for columns isolating an eigenvalue
c                and push them left ..........
  130 k = k + 1
c
  140 do 170 j = k, l
c
         do 150 i = k, l
            if (i .eq. j) go to 150
            if (a(i,j) .ne. 0.) go to 170
  150    continue
c
         m = k
         iexc = 2
         go to 20
  170 continue
c     .......... now balance the submatrix in rows k to l ..........
      do 180 i = k, l
  180 scale(i) = 1.0
c     .......... iterative loop for norm reduction ..........
  190 noconv = .false.
c
      do 270 i = k, l
         c = 0.
         r = 0.
c
         do 200 j = k, l
            if (j .eq. i) go to 200
            c = c + abs(a(j,i))
            r = r + abs(a(i,j))
  200    continue
c     .......... guard against zero c or r due to underflow ..........
         if (c .eq. 0. .or. r .eq. 0.) go to 270
         g = r / radix
         f = 1.0
         s = c + r
  210    if (c .ge. g) go to 220
         f = f * radix
         c = c * b2
         go to 210
  220    g = r * radix
  230    if (c .lt. g) go to 240
         f = f / radix
         c = c / b2
         go to 230
c     .......... now balance ..........
  240    if ((c + r) / f .ge. 0.95 * s) go to 270
         g = 1.0 / f
         scale(i) = scale(i) * f
         noconv = .true.
c
         do 250 j = k, n
  250    a(i,j) = a(i,j) * g
c
         do 260 j = 1, l
  260    a(j,i) = a(j,i) * f
c
  270 continue
c
      if (noconv) go to 190
c
  280 low = k
      igh = l
      return
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: elmhes
C TYPE   : subroutine
C PURPOSE: reduction to hessenberg form
C I/O    :
C VERSION: 30-Sep-95
C COMMENT: netlib/handbook for auto. comp. 
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine elmhes(nm,n,low,igh,a,int)
      implicit double precision (a-h,o-z)
c
      integer i,j,m,n,la,nm,igh,kp1,low,mm1,mp1
      dimension a(nm,n)
      integer int(igh)
c
c     this subroutine is a translation of the algol procedure elmhes,
c     num. math. 12, 349-368(1968) by martin and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 339-358(1971).
c
c     given a real general matrix, this subroutine
c     reduces a submatrix situated in rows and columns
c     low through igh to upper hessenberg form by
c     stabilized elementary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        a contains the input matrix.
c
c     on output
c
c        a contains the hessenberg matrix.  the multipliers
c          which were used in the reduction are stored in the
c          remaining triangle under the hessenberg matrix.
c
c        int contains information on the rows and columns
c          interchanged in the reduction.
c          only elements low through igh are used.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      la = igh - 1
      kp1 = low + 1
      if (la .lt. kp1) go to 200
c
      do 180 m = kp1, la
         mm1 = m - 1
         x = 0.
         i = m
c
         do 100 j = m, igh
            if (abs(a(j,mm1)) .le. abs(x)) go to 100
            x = a(j,mm1)
            i = j
  100    continue
c
         int(m) = i
         if (i .eq. m) go to 130
c     .......... interchange rows and columns of a ..........
         do 110 j = mm1, n
            y = a(i,j)
            a(i,j) = a(m,j)
            a(m,j) = y
  110    continue
c
         do 120 j = 1, igh
            y = a(j,i)
            a(j,i) = a(j,m)
            a(j,m) = y
  120    continue
c     .......... end interchange ..........
  130    if (x .eq. 0.) go to 180
         mp1 = m + 1
c
         do 160 i = mp1, igh
            y = a(i,mm1)
            if (y .eq. 0.) go to 160
            y = y / x
            a(i,mm1) = y
c
            do 140 j = m, n
  140       a(i,j) = a(i,j) - y * a(m,j)
c
            do 150 j = 1, igh
  150       a(j,m) = a(j,m) + y * a(j,i)
c
  160    continue
c
  180 continue
c
  200 return
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: hqr
C TYPE   : subroutine
C PURPOSE: eigenvalues of a real
C          upper hessenberg matrix by the qr method
C I/O    :
C VERSION: 30-Sep-95
C COMMENT: netlib/handbook for auto. comp. 
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)
      implicit double precision (a-h,o-z)
C  RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG)
c
      integer i,j,k,l,m,n,en,ll,mm,na,nm,igh,itn,its,low,mp2,enm2,ierr
      dimension h(nm,n),wr(n),wi(n)
      logical notlas
c
c     this subroutine is a translation of the algol procedure hqr,
c     num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
c
c     this subroutine finds the eigenvalues of a real
c     upper hessenberg matrix by the qr method.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        low and igh are integers determined by the balancing
c          subroutine  balanc.  if  balanc  has not been used,
c          set low=1, igh=n.
c
c        h contains the upper hessenberg matrix.  information about
c          the transformations used in the reduction to hessenberg
c          form by  elmhes  or  orthes, if performed, is stored
c          in the remaining triangle under the hessenberg matrix.
c
c     on output
c
c        h has been destroyed.  therefore, it must be saved
c          before calling  hqr  if subsequent calculation and
c          back transformation of eigenvectors is to be performed.
c
c        wr and wi contain the real and imaginary parts,
c          respectively, of the eigenvalues.  the eigenvalues
c          are unordered except that complex conjugate pairs
c          of values appear consecutively with the eigenvalue
c          having the positive imaginary part first.  if an
c          error exit is made, the eigenvalues should be correct
c          for indices ierr+1,...,n.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated september 1989.
c
c     ------------------------------------------------------------------
c
      ierr = 0
      xnorm = 0.
      k = 1
c     .......... store roots isolated by balanc
c                and compute matrix norm ..........
      do 50 i = 1, n
c
         do 40 j = k, n
   40    xnorm = xnorm + abs(h(i,j))
c
         k = i
         if (i .ge. low .and. i .le. igh) go to 50
         wr(i) = h(i,i)
         wi(i) = 0.
   50 continue
c
      en = igh
      t = 0.
      itn = 30*n
c     .......... search for next eigenvalues ..........
   60 if (en .lt. low) go to 1001
      its = 0
      na = en - 1
      enm2 = na - 1
c     .......... look for single small sub-diagonal element
c                for l=en step -1 until low do -- ..........
   70 do 80 ll = low, en
         l = en + low - ll
         if (l .eq. low) go to 100
         s = abs(h(l-1,l-1)) + abs(h(l,l))
         if (s .eq. 0.) s = xnorm
         tst1 = s
         tst2 = tst1 + abs(h(l,l-1))
         if (tst2 .eq. tst1) go to 100
   80 continue
c     .......... form shift ..........
  100 x = h(en,en)
      if (l .eq. en) go to 270
      y = h(na,na)
      w = h(en,na) * h(na,en)
      if (l .eq. na) go to 280
      if (itn .eq. 0) go to 1000
      if (its .ne. 10 .and. its .ne. 20) go to 130
c     .......... form exceptional shift ..........
      t = t + x
c
      do 120 i = low, en
  120 h(i,i) = h(i,i) - x
c
      s = abs(h(en,na)) + abs(h(na,enm2))
      x = 0.75 * s
      y = x
      w = -0.4375 * s * s
  130 its = its + 1
      itn = itn - 1
c     .......... look for two consecutive small
c                sub-diagonal elements.
c                for m=en-2 step -1 until l do -- ..........
      do 140 mm = l, enm2
         m = enm2 + l - mm
         zz = h(m,m)
         r = x - zz
         s = y - zz
         p = (r * s - w) / h(m+1,m) + h(m,m+1)
         q = h(m+1,m+1) - zz - r - s
         r = h(m+2,m+1)
         s = abs(p) + abs(q) + abs(r)
         p = p / s
         q = q / s
         r = r / s
         if (m .eq. l) go to 150
         tst1 = abs(p)*(abs(h(m-1,m-1)) + abs(zz) + abs(h(m+1,m+1)))
         tst2 = tst1 + abs(h(m,m-1))*(abs(q) + abs(r))
         if (tst2 .eq. tst1) go to 150
  140 continue
c
  150 mp2 = m + 2
c
      do 160 i = mp2, en
         h(i,i-2) = 0.
         if (i .eq. mp2) go to 160
         h(i,i-3) = 0.
  160 continue
c     .......... double qr step involving rows l to en and
c                columns m to en ..........
      do 260 k = m, na
         notlas = k .ne. na
         if (k .eq. m) go to 170
         p = h(k,k-1)
         q = h(k+1,k-1)
         r = 0.
         if (notlas) r = h(k+2,k-1)
         x = abs(p) + abs(q) + abs(r)
         if (x .eq. 0.) go to 260
         p = p / x
         q = q / x
         r = r / x
  170    s = sign(sqrt(p*p+q*q+r*r),p)
         if (k .eq. m) go to 180
         h(k,k-1) = -s * x
         go to 190
  180    if (l .ne. m) h(k,k-1) = -h(k,k-1)
  190    p = p + s
         x = p / s
         y = q / s
         zz = r / s
         q = q / p
         r = r / p
         if (notlas) go to 225
c     .......... row modification ..........
         do 200 j = k, EN
            p = h(k,j) + q * h(k+1,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
  200    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 210 i = L, j
            p = x * h(i,k) + y * h(i,k+1)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
  210    continue
         go to 255
  225    continue
c     .......... row modification ..........
         do 230 j = k, EN
            p = h(k,j) + q * h(k+1,j) + r * h(k+2,j)
            h(k,j) = h(k,j) - p * x
            h(k+1,j) = h(k+1,j) - p * y
            h(k+2,j) = h(k+2,j) - p * zz
  230    continue
c
         j = min0(en,k+3)
c     .......... column modification ..........
         do 240 i = L, j
            p = x * h(i,k) + y * h(i,k+1) + zz * h(i,k+2)
            h(i,k) = h(i,k) - p
            h(i,k+1) = h(i,k+1) - p * q
            h(i,k+2) = h(i,k+2) - p * r
  240    continue
  255    continue
c
  260 continue
c
      go to 70
c     .......... one root found ..........
  270 wr(en) = x + t
      wi(en) = 0.
      en = na
      go to 60
c     .......... two roots found ..........
  280 p = (y - x) / 2.0
      q = p * p + w
      zz = sqrt(abs(q))
      x = x + t
      if (q .lt. 0.) go to 320
c     .......... real pair ..........
      zz = p + sign(zz,p)
      wr(na) = x + zz
      wr(en) = wr(na)
      if (zz .ne. 0.) wr(en) = x - w / zz
      wi(na) = 0.
      wi(en) = 0.
      go to 330
c     .......... complex pair ..........
  320 wr(na) = x + p
      wr(en) = x + p
      wi(na) = zz
      wi(en) = -zz
  330 en = enm2
      go to 60
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
 1001 return
      end
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: inverse.f
C TYPE   : subroutine
C PURPOSE: calculate inverse of matrix
C VERSION: 30-Sep-95
C COMMENT: from numerical recipes	, y is inverse of a
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine inverse(a,n,np,y)
      implicit double precision (a-h,o-z)
      parameter (npp=300)
      dimension a(np,np),y(np,np),indx(npp)
      dimension b(npp,npp)
      if (npp.lt.n) pause ' size of matrix too large'
      do i=1,n
         do j=1,n
            b(i,j)=a(i,j)
            y(i,j)=0.
         end do
         y(i,i)=1.
      end do
      call ludcmp(b,n,npp,indx,d)
      do 3 i=1,n
           call lubksb(b,n,npp,indx,y(1,i))
3     continue
      end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: transfer 
C TYPE   : subroutine
C PURPOSE: transfer n x n matrices A --> B  (physical size np x np) 
C VERSION: 30-Sep-95
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
        subroutine transfer (a,b,n)
        implicit double precision (a-h,o-z)
        dimension a(n,n),b(n,n)
        do i=1,n
           do j=1,n
              b(i,j)=a(i,j)
           end do
        end do
        end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: null
C TYPE   : subroutine
C PURPOSE: set the n x n matrix A to zero: A = 0
C VERSION: 30-Sep-95
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
        subroutine null (a,n)
        implicit double precision (a-h,o-z)
        dimension a(n,n)
        do i=1,n
           do j=1,n
              a(i,j)=0
           end do
        end do
        end 
C========+=========+=========+=========+=========+=========+=========+=$
C PROGRAM: lubksb.f
C TYPE   : subroutine
C PURPOSE: 
C I/O    :
C VERSION: 30-Sep-95
C COMMENT: from numerical recipes	, y is inverse of a
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine lubksb(a,n,np,indx,b)
      implicit double precision (a-h,o-z)
      dimension a(np,np),indx(n),b(n)
      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          end do
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
      end do
      do  i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do  J=I+1,N
            sum=sum-a(i,j)*b(j)
          end do
        endif
        b(i)=sum/a(i,i)
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
C       PROGRAM: ludcmp
C       TYPE   : subroutine
C       PURPOSE: 
C       I/O    :
C       VERSION: 30-Sep-95
C       COMMENT: from numerical recipes	, y is inverse of a
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine ludcmp(a,n,np,indx,d)
      implicit double precision (a-h,o-z)
      parameter (nmax=300,tiny=1.0e-20)
      dimension a(np,np),indx(n),vv(nmax)
      if (nmax.lt.n) pause ' matrix too large, see in ludcmp'
      D=1.
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
