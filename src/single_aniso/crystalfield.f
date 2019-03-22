************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine CRYSTALFIELD( ESOJ, DIPSO, S_SO, nDIMcf,
     &                         iDIM, nlanth, zmagn2, iopt, GRAD, iprint)

      Implicit None
#include "stdalloc.fh"
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)           :: iprint, iDIM
      Integer, intent(in)           :: nDIMcf, iopt, nlanth
      Real(kind=wp), intent(in)     :: ESOJ(nDIMcf)
      Real(kind=wp), intent(in)     :: ZMAGN2(3,3)
      Complex(kind=wp), intent(in)  :: DIPSO(3,nDIMcf,nDIMcf)
      Complex(kind=wp), intent(in)  ::  S_SO(3,nDIMcf,nDIMcf)
      Logical, intent(in)           :: GRAD
      ! local variables
      Integer :: info
      Real(kind=wp), allocatable :: wtmp(:)
      Complex(kind=wp), allocatable ::  DIPJ(:,:,:)
      Complex(kind=wp), allocatable ::    SJ(:,:,:), ztmp(:,:)
      Integer                       :: i,j
      Real(kind=wp), allocatable    :: gtens(:), zmagn(:,:)

      Call qEnter('SA_CF')

      Write(6,'(/)')
      Write(6,'(100A)') ('%',i=1,95)
      If(MOD(nDIMcf,2).eq.1) Then
        Write(6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD '//
     &                         'PARAMETERS OF THE GROUND ATOMIC '//
     &                         'MULTIPLET J = ', (nDIMcf-1)/2, '.'
      Else
        Write(6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD '//
     &                         'PARAMETERS OF THE GROUND ATOMIC '//
     &                         'MULTIPLET J = ', (nDIMcf-1),'/2.'
      End If
      Write(6,'(100A)') ('%',i=1,95)
      Write(6,*)


      Call mma_allocate(gtens,3,'gtens')
      Call mma_allocate(zmagn,3,3,'zmagn')
      Call dcopy_(  3,0.0_wp,0,gtens,1)
      Call dcopy_(3*3,0.0_wp,0,zmagn,1)
      If(iopt.eq.1) Then
c  coordinate system for decomposition of the CF matrix identic to the coordinate system
c  of the main magnetic axes of the ground multiplet (NDIM(1))
        CALL atens(DIPSO(1:3,1:idim,1:idim), idim, GTENS, ZMAGN, 1)
        Write(6,'(a)') 'The parameters of the Crystal Field matrix '//
     &                 'are written in the coordinate system:'
        If(MOD(iDIM,2).eq.0) Then
          Write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic '//
     &                        'axes of the ground pseuDospin S = |',
     &                        iDIM-1,'/2> multiplet.'
        Else
          Write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic '//
     &                        'axes of the ground pseuDospin S = |',
     &                       (iDIM-1)/2,'> multiplet.'
        End If

      Else If(iopt.eq.2) Then
c  coordinate system for decomposition of the CF matrix identic to the coordinate system
c  of the main magnetic axes of the ground multiplet (NDIM(1))
        CALL atens(DIPSO, nDIMCF, GTENS, ZMAGN, 1)
        Write(6,'(a)') 'The parameters of the Crystal Field matrix '//
     &                 'are written in the coordinate system:'
        If(MOD(nDIMCF,2).eq.0) Then
          Write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic '//
     &                        'axes of the ground atomic J = |',
     &                         nDIMCF-1,'/2> multiplet'
        Else
          Write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic '//
     &                        'axes of the ground atomic J = |',
     &                        (nDIMCF-1)/2,'> multiplet'
        End If

      Else If(iopt.eq.3) Then
        Write(6,'(a)') 'The parameters of the Crystal Field matrix '//
     &                 'are written in the coordinate system:'
        Write(6,'(a)') '(Xm, Ym, Zm) -- defined in the input file.'
        Call dcopy_(3*3,zmagn2,1,zmagn,1)
      Else
        Write(6,'(a)') 'The parameters of the Crystal Field matrix '//
     &                 'are written in the initial coordinate system.'
        Do i=1,3
          ZMAGN(i,i)=1.0_wp
        End Do
      End If !axisoption


      ! rotate the momentum:
      Call mma_allocate(DIPJ,3,nDIMcf,nDIMcf,'DIPJ')
      Call mma_allocate(SJ,3,nDIMcf,nDIMcf,'SJ')
      Call zcopy_(3*nDIMcf*nDIMcf,(0.0_wp,0.0_wp),0,DIPJ,1)
      Call zcopy_(3*nDIMcf*nDIMcf,(0.0_wp,0.0_wp),0,  SJ,1)
      CALL rotmom2( DIPSO, nDIMCF, ZMAGN, DIPJ )
      CALL rotmom2(  S_SO, nDIMCF, ZMAGN,   SJ )

      Write(6,'(a)') 'Rotation matrix from the initial coordinate '//
     &               'system to the employed coordinate system is:'

      If((iopt.eq.1).OR.(iopt.eq.2)) Then
        Write(6,'(70a)') ('-',i=1,67),'|'
        Write(6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
        Write(6,'(A,35x,A)') 'Xm, Ym, Zm -- main magnetic axes','|'
        Write(6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
        Write(6,'(6x,A,3F18.14,1x,A)') '| Xm |',(ZMAGN(j,1),j=1,3),'|'
        Write(6,'( A,A,3F18.14,1x,A)') ' R =  ','| Ym |',
     &                                 (ZMAGN(j,2),j=1,3),'|'
        Write(6,'(6x,A,3F18.14,1x,A)') '| Zm |',(ZMAGN(j,3),j=1,3),'|'
        Write(6,'(83a)') ('-',i=1,67),'|'
        Write(6,'(A,I3)') 'Quantization axis is Zm.'

      Else If(iopt.eq.3) Then

        Write(6,'(70a)') ('-',i=1,67),'|'
        Write(6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
        Write(6,'(A,11x,A)') 'Xm, Ym, Zm -- the coordinate system '//
     &                       'defined in the input','|'
        Write(6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
        Write(6,'(6x,A,3F18.14,1x,A)') '| Xm |',(ZMAGN(j,1),j=1,3),'|'
        Write(6,'( A,A,3F18.14,1x,A)') ' R =  ','| Ym |',
     &   (ZMAGN(j,2),j=1,3),'|'
        Write(6,'(6x,A,3F18.14,1x,A)') '| Zm |',(ZMAGN(j,3),j=1,3),'|'
        Write(6,'(83a)') ('-',i=1,67),'|'
        Write(6,'(A,I3)') 'Quantization axis is Zm.'

      Else

        Write(6,'(A)') 'IDENTITY matrix.'
        Write(6,'(A,I3)') 'Quantization axis is the initial z axis.'
      End If

      If(IPRINT.gt.2) Then
        CALL prMom('CRYSTALFIELD::   DIPJ(l,i,j)',DIPJ,nDIMcf)
        CALL prMom('CRYSTALFIELD::     SJ(l,i,j)',  SJ,nDIMcf)
      End If

      Call mma_allocate(ztmp,nDIMcf,nDIMcf,'z')
      Call mma_allocate(wtmp,nDIMcf,'w')
      wtmp(:)=0.0_wp
      ztmp(:,:)=(0.0_wp,0.0_wp)
      info=0
      Call diag_c2(DIPJ,nDIMcf,info,wtmp,ztmp)
      Do i=1,nDIMcf
        Write(6,'(A,i2,A,4F20.15)') 'energy: ',i,' : ',
     &             wtmp(i), wtmp(i)+wtmp(nDIMcf-i+1)
      End Do
      Call mma_deallocate(ztmp)
      Call mma_deallocate(wtmp)

      CALL CRYSTALFIELD_1(nDIMcf,nlanth,DIPJ,ESOJ,GRAD,iprint)


      Call mma_deallocate(DIPJ)
      Call mma_deallocate(SJ)
      Call mma_deallocate(gtens)
      Call mma_deallocate(zmagn)
      Call qExit('SA_CF')
      Return
      End



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine CRYSTALFIELD_1(nDIMcf,nlanth,MM,ESOJ,GRAD,iprint)

C This soubrutine calculates the crystal field parameters on the basis
C of the given fron RASSI - J multiplet.
c In a second step, the first largest 27 parameters will be used to
c recalculate the S-O energies, eigenfunctions, g- and D- tensors.

c  Employed parameters:


C  IPRINT = the print level of the calculation


C  IReturn = the error value.
c        0 = no error, happy landing

C================== Variable declarations =============================

      Implicit None
#include "stdalloc.fh"
      Integer, Parameter            :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)           :: nDIMcf, nlanth, iprint
      Logical, intent(in)           :: GRAD
      Real(kind=wp), intent(in)     :: ESOJ(nDIMcf)
      Complex(kind=wp), intent(in)  :: MM(3,nDIMcf,nDIMcf)
      ! local variables:
      Integer                       :: info, i, j, k, q
      Integer                       :: LuCF, IsFreeUnit
      Real(kind=wp)                 :: dznrm2
      Real(kind=wp), allocatable    :: Winit(:), Eloc(:), a(:)
      Real(kind=wp)                 :: BNC(nDIMcf,0:nDIMcf)
      Real(kind=wp)                 :: BNS(nDIMcf,0:nDIMcf)
      Real(kind=wp)                 :: Bstev(nDIMcf,-nDIMcf:nDIMcf)
      Complex(kind=wp)              :: trace
      Complex(kind=wp)              ::
     &                   Akq((nDIMcf-1),-(nDIMcf-1):(nDIMcf-1))
      Complex(kind=wp), allocatable :: Zinit(:,:), Z(:,:),
     &                                  HCF(:,:)
      External           :: trace, dznrm2, IsFreeUnit

      Call qEnter('SA_CF1')
C============== End of variable declarations ==========================
      Call mma_allocate(Winit,nDIMcf,'Winit')
      Call mma_allocate(Eloc,nDIMcf,'Eloc')
      Call mma_allocate(a,6,'anm')
      Call mma_allocate(Zinit,nDIMcf,nDIMcf,'Zinit')
      Call mma_allocate(Z,nDIMcf,nDIMcf,'Z')
      Call mma_allocate(HCF,nDIMcf,nDIMcf,'HCF')

      info=0
      Call dcopy_(nDIMcf,0.0_wp,0,Winit,1)
      Call dcopy_(nDIMcf,0.0_wp,0,Eloc,1)
      Call dcopy_(6,0.0_wp,0,A,1)
      Call zcopy_(nDIMcf*nDIMcf,(0.0_wp,0.0_wp),0,Zinit,1)
      Call zcopy_(nDIMcf*nDIMcf,(0.0_wp,0.0_wp),0,Z,1)
      Call zcopy_(nDIMcf*nDIMcf,(0.0_wp,0.0_wp),0,HCF,1)

      ! find the J-pseudospin:
      !iDir=3
      CALL pseudospin(MM,nDIMcf,Z,3,1,iprint)
      CALL rtrace(nDIMcf,ESOJ,ELOC)
      ! re-write the CF matrix in J-pseudospin basis:
      ! energy units =  cm-1
      Do i=1,nDIMcf
        Do j=1,nDIMcf
          Do k=1,nDIMcf
            HCF(i,j)=HCF(i,j) + ELOC(k)*CONJG(Z(k,i))*Z(k,j)
          End Do
        End Do
      End Do
      ! diagonalize the initial CF matrix:
      CALL DIAG_C2( HCF,nDIMcf,INFO,Winit,Zinit)
      Call print_ZFS('Ab Initio Calculated Crystal-Field Splitting '//
     &               'Matrix written in the basis of Pseudospin '//
     &               'Eigenfunctions',HCF,nDIMCF)

      If(IPRINT.gt.2) Then
        Write(6,*)
        Write(6,'(5X,A)') 'MAIN VALUES OF THE INITIAL CRYSTAL'//
     &                    '-FIELD HAMILTONIAN:'
        Write(6,*)
        If(MOD(nDIMcf,2).eq.1) Then
          Do I=1,nDIMcf
            Write(6,'(3X,A,I3,A,F25.16)') '|',
     &                     (nDIMcf-1)/2+(1-I),'> = ',Winit(I)-Winit(1)
          End Do
        Else
          Do I=1,nDIMcf
            Write(6,'(3X,A,I3,A,F25.16)') '|',
     &                  (nDIMcf-1)-2*(I-1),'/2 > = ',Winit(i)-Winit(1)
          End Do
        End If
        Write(6,*)

        Write(6,'(5X,A)') 'EIGENVECTORS OF THE INITIAL CRYSTAL'//
     &                    '-FIELD HAMILTONIAN:'
        Write(6,*)
        Call print_ZFS_naoya('J',Zinit,nDIMcf)
!     End  the checking of the main values of the initial crystal-field
      End If

C  calculating the coeficients of the crystal filed operators Bnm
C    Akq=(2k+1)/(2J+1) * 1/|< J || O || J >|^2 * Trace{HCF*O(k,-q)}
      Call NEWCF(HCF,nDIMcf,Akq,BNC,BNS,Bstev)
      !If(dbg) Call recover_CF(nDIMCF,HCF,Akq,BNC,BNS,Bstev)

      Call print_CFP_alpha(nlanth,nDIMCF,BNC,BNS)
      If(iprint>=4) Then
         Call print_CFP_LCLU(nDIMCF,BNC,BNS,.true.)
         Call print_CFP_stev(nDIMCF,Bstev,.true.)
         Call print_CFP_naoya(nDIMcf,Akq,.true.)
      Else
         Call print_CFP_LCLU(nDIMCF,BNC,BNS,.false.)
         Call print_CFP_stev(nDIMCF,Bstev,.false.)
         Call print_CFP_naoya(nDIMcf,Akq,.false.)
      End If

c-----------------------------------------------------------------------
      Write(6,'(/)')
      If(MOD(nDIMcf,2)==1) Then
         Write(6,'(A,I0)') 'DECOMPOSITION OF THE RASSI WAVE '//
     &                     'FUNCTIONS CORRESPONDING TO THE '//
     &                     'LOWEST ATOMIC MULTIPLET J =',(nDIMcf-1)/2
         Write(6,'(A,I0)') 'IN WAVE FUNCTIONS WITH DEFINITE '//
     &                     'PROJECTION OF THE TOTAL MOMENT '//
     &                     'ON THE QUANTIZATION AXIS'
      Else ! MOD(nDIMcf,2)==0
         Write(6,'(A,I0,A)') 'DECOMPOSITION OF THE RASSI WAVE '//
     &                       'FUNCTIONS CORRESPONDING TO THE '//
     &                       'LOWEST ATOMIC MULTIPLET J = ',
     &                         (nDIMcf-1),'/2'
         Write(6,'(A,I0)') 'IN WAVE FUNCTIONS WITH DEFINITE '//
     &                     'PROJECTION OF THE TOTAL MOMENT '//
     &                     'ON THE QUANTIZATION AXIS'
      End If

      Call print_ZFS_naoya('J',Zinit,nDIMcf)
      Call individual_ranks(nDIMCF,BNC,BNS,HCF,'J',iprint)
c-----------------------------------------------------------------------
C  saving some information for tests:
      CALL Add_Info('CRYS_BNMC_20',DBLE(BNC(2,0)),1,4)
      CALL Add_Info('CRYS_BNMC_40',DBLE(BNC(4,0)),1,4)
      CALL Add_Info('CRYS_BNMC_60',DBLE(BNC(6,0)),1,4)
c-----------------------------------------------------------------------
      ! for the interface related to CF gradient calculation:
      If (GRAD) Then
         LuCF=IsFreeUnit(81)
         Call molcas_open(LuCF,'CFMAT')
         Do k=2,nDIMcf-1,2
           Do q=0,k
               Write(LuCF,'(I3,I3,1x,2ES25.15)') k,q,BNC(k,q),BNS(k,q)
           End Do
         End Do
         Close(LuCF)
      End If
c-----------------------------------------------------------------------
      Call mma_deallocate(Winit)
      Call mma_deallocate(Eloc)
      Call mma_deallocate(a)
      Call mma_deallocate(Zinit)
      Call mma_deallocate(Z)
      Call mma_deallocate(HCF)
      Call qExit('SA_CF1')

      Return
      End









      Subroutine newCF(H,n,A, B,C,Bstev)
      Implicit none
      Integer, Parameter          :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)           :: n
      Complex(kind=wp),intent(in)   :: H(n,n)
      Complex(kind=wp), intent(out) :: A( (n-1), -(n-1):(n-1) )
      Real(kind=wp), intent(out)    :: B(n,0:n), C(n,0:n)
      Real(kind=wp), intent(out)    :: Bstev(n,-n:n)
      ! local variables:
      Integer                       :: ik,iq
      Real(kind=wp)                 :: rfact,cr,mfact,C0
      Complex(kind=wp)              :: trace, cfact
      Complex(kind=wp), allocatable :: Cp(:,:), Cm(:,:)
      Complex(kind=wp)              :: mf
      Real(kind=wp)                 :: knm(12,0:12)
      External                      :: trace
      Logical                       :: dbg

      Call qEnter('SA_newCF')
!-------------------------------------------
      If(n<1) Return
!-------------------------------------------
      dbg=.false.
      Call mma_allocate(Cp,n,n,'operator O')
      Call mma_allocate(Cm,n,n,'operator W')
!-------------------------------------------
!     n=2*J+1;  or   n=2*S+1
      Bstev(1:n,-n:n)=0.0_wp
      B(1:n,0:n)=0.0_wp
      C(1:n,0:n)=0.0_wp
      A(1:(n-1),-(n-1):(n-1))=(0.0_wp,0.0_wp)
      Call set_knm( knm )

      Do ik=1,n-1
         Do iq=0,ik
            cr=0.0_wp
            mfact=0.0_wp
            rfact=0.0_wp
            cfact=(0.0_wp,0.0_wp)
            C0=0.0_wp
            ! generate the operator matrix K=ik, Q=iq, dimension = n
            Call ITO(n,ik,iq,C0,Cp,Cm)
            Call coeff_redus_sub(n,ik,cr)

            mfact=dble((-1)**iq)
            rfact=C0*C0*dble(2*ik+1)/dble(n)
            cfact=cmplx(mfact*rfact,0.0_wp,wp)


            !-------------------------------------------
            ! Naoya's C/C0 operators:
            a(ik,-iq)=cfact*trace(n,H,Cp)
            a(ik, iq)=cfact*trace(n,H,Cm)


            !-------------------------------------------
            ! make real combinations of CF parameters:
            ! Liviu's ITO operators:
            If(iq==0) Then
               !b(ik, iq)=dble( (0.5_wp,0.0_wp)*(A(ik,iq)+A(ik,-iq)) )
               b(ik, iq)=dble(A(ik,iq))
            Else
               mf=cmplx((-1)**iq,0.0_wp,wp)
               b(ik, iq)=dble( A(ik,-iq)+mf*A(ik,iq))
               c(ik, iq)=dble((A(ik,-iq)-mf*A(ik,iq))*(0.0_wp,-1.0_wp))
            End If
            ! scale with the correct ratio:
            b(ik,iq)=b(ik,iq)/(cr*C0)
            c(ik,iq)=c(ik,iq)/(cr*C0)


            !-------------------------------------------
            ! parameters to be used in connection with ESO as in MATLAB
            ! EasySpin program
            ! The ESO operaors are not implemented for k>12 in EasySpin
            ! therefore we do not provide these parameters as well.
            If((ik<=12).and.(iq<=12)) Then
               ! scale with the correct ratio:
               If(iq==0) Then
                  bstev(ik, iq)=b(ik,iq)*knm(ik,iq)
               Else
                  bstev(ik, iq)=b(ik,iq)*knm(ik,iq)
                  bstev(ik,-iq)=c(ik,iq)*knm(ik,iq)
               End If
            End If

            If(dbg) Then
               Write(6,'(A,2I3,5(ES20.13,1x))') 'k,q, b(k,q), c(k,q)',
     &                                      ik,iq, b(ik,iq), c(ik,iq)
            End If
         End Do
      End Do

      Call mma_deallocate(Cp)
      Call mma_deallocate(Cm)
      Call qExit('SA_newCF')

      Return
      End subroutine newCF







      Subroutine recover_CF(N,HAM,Akq,B,C,Bstev)
      Implicit none
      Integer, Parameter           :: wp=selected_real_kind(p=15,r=307)
      Integer, intent(in)          :: n
      Complex(kind=wp), intent(in) :: HAM(n,n)
      Complex(kind=wp), intent(in) :: Akq((n-1), -(n-1):(n-1))
      Real(kind=wp), intent(in)    :: B(n,0:n), C(n,0:n), Bstev(n,-n:n)

      Integer          :: k,q,i,j,info
      Real(kind=wp)    :: tdiff
      Complex(kind=wp) :: Cp(n,n), Cm(n,n), redME
      Complex(kind=wp) :: O(n,n), W(n,n),zfact
      Complex(kind=wp) :: HCF(n,n), Z(n,n)
      Real(kind=wp)    :: w1(n), w2(n),c0,dznrm2_
      External         :: dznrm2_

      Do k=2,n-1
        Do q=-k,k
          Write(6,'(A,i2,A,i3,A,2ES20.10)') 'Akq(',k,',',q,') = ',
     &                                       Akq(k,q)
        End Do
      End Do
!==================================================================
      Write(6,'(A,ES20.10)') 'recover from Akq parameters'
      tdiff=0.0_wp
      Call zcopy_(n*n,(0.0_wp,0.0_wp),0,HCF,1)
      Do k=1,n-1
        Do q=0,k
          ! generate the operator matrix K=ik, Q=iq, dimension=na
          Call ITO(n,k,q,C0,Cp,Cm)
          If(q==0) Then
            Call zaxpy_(n*n, Akq(k, q), Cp, 1, HCF,1)
          Else
            Call zaxpy_(n*n, Akq(k, q), Cp, 1, HCF,1)
            Call zaxpy_(n*n, Akq(k,-q), Cm, 1, HCF,1)
          End If
        End Do !q
      End Do !k
      tdiff=dznrm2_(n*n,HAM-HCF,1)
      Write(6,'(A,ES20.10)') 'total difference between HAM-HCF=',tdiff
      Do i=1,n
        Do j=1,n
          Write(6, '(2(A,i2,A,i2,A,2ES20.10,A),2(2ES20.10,5x))')
     &                'HAM(',i,',',j,')=',HAM(i,j),'      ',
     &                'HCF(',i,',',j,')=',HCF(i,j),' diff=',
     &                 HAM(i,j)-HCF(i,j)
        End Do
      End Do
      w1(:)=0.0_wp
      Z(:,:)=(0.0_wp,0.0_wp)
      Call diag_c2(HAM,n,info,w1,Z)
      w2(:)=0.0_wp
      Z(:,:)=(0.0_wp,0.0_wp)
      Call diag_c2(HCF,n,info,w2,Z)
      Do i=1,n
         Write(6,'(2(A,i2,A,ES20.10,A),2(2ES20.10,5x))')
     &                'W1(',i,')=',w1(i)-w1(1),'      ',
     &                'W2(',i,')=',w2(i)-w2(1),' diff=',
     &                 w1(i)-w2(i)
      End Do

!==================================================================
      Write(6,'(A,ES20.10)') 'recover from B and C parameters'
      tdiff=0.0_wp
      Call zcopy_(n*n,(0.0_wp,0.0_wp),0,HCF,1)
      Do k=1,n-1
        Do q=0,k
          Call Liviu_ESO(n,k,q,O,W,redME)
          If(q==0) Then
            zfact=cmplx(B(k,0),0.0_wp,wp)
            Call zaxpy_(n*n, zfact, O, 1,  HCF,1)
          Else
            zfact=cmplx(B(k,q),0.0_wp,wp)
            Call zaxpy_(n*n,zfact,O,1,HCF,1)
            zfact=cmplx(C(k,q),0.0_wp,wp)
            Call zaxpy_(n*n,zfact,W,1,HCF,1)
          End If
        End Do
      End Do
      tdiff=dznrm2_(n*n,(HAM(1:n,1:n)-HCF(1:n,1:n)),1)
      Write(6,'(A,ES20.10)') 'total difference between HAM-HCF=',tdiff
      Do i=1,n
        Do j=1,n
          Write(6, '(2(A,i2,A,i2,A,2ES20.10,A),2(2ES20.10,5x))')
     &                'HAM(',i,',',j,')=',HAM(i,j),'      ',
     &                'HCF(',i,',',j,')=',HCF(i,j),' diff=',
     &                 HAM(i,j)-HCF(i,j)
        End Do
      End Do
      w1(:)=0.0_wp
      Z(:,:)=(0.0_wp,0.0_wp)
      Call diag_c2(HAM,n,info,w1,Z)
      w2(:)=0.0_wp
      Z(:,:)=(0.0_wp,0.0_wp)
      Call diag_c2(HCF,n,info,w2,Z)
      Do i=1,n
         Write(6,'(2(A,i2,A,ES20.10,A),2(2ES20.10,5x))')
     &                'W1(',i,')=',w1(i)-w1(1),'      ',
     &                'W2(',i,')=',w2(i)-w2(1),' diff=',
     &                 w1(i)-w2(i)
      End Do

!==================================================================
      Write(6,'(A,ES20.10)') 'recover from Bstev'
      tdiff=0.0_wp
      Call zcopy_(n*n,(0.0_wp,0.0_wp),0,HCF,1)
      Do k=1,n-1
        Do q=0,k
          Call ESO(n,k,q,O,W,redME)
          If(q==0) Then
            zfact=cmplx(Bstev(k,0),0.0_wp,wp)
            Call zaxpy_(n*n, zfact, O, 1,  HCF,1)
          Else
            zfact=cmplx(Bstev(k, q),0.0_wp,wp)
            Call zaxpy_(n*n,zfact,O,1,HCF,1)
            zfact=cmplx(Bstev(k,-q),0.0_wp,wp)
            Call zaxpy_(n*n,zfact,W,1,HCF,1)
          End If
        End Do
      End Do
      tdiff=dznrm2_(n*n,(HAM(1:n,1:n)-HCF(1:n,1:n)),1)
      Write(6,'(A,ES20.10)') 'total difference between HAM-HCF=',tdiff
      Do i=1,n
        Do j=1,n
          Write(6, '(2(A,i2,A,i2,A,2ES20.10,A),2(2ES20.10,5x))')
     &                'HAM(',i,',',j,')=',HAM(i,j),'      ',
     &                'HCF(',i,',',j,')=',HCF(i,j),' diff=',
     &                 HAM(i,j)-HCF(i,j)
        End Do
      End Do
      w1(:)=0.0_wp
      Z(:,:)=(0.0_wp,0.0_wp)
      Call diag_c2(HAM,n,info,w1,Z)
      w2(:)=0.0_wp
      Z(:,:)=(0.0_wp,0.0_wp)
      Call diag_c2(HCF,n,info,w2,Z)
      Do i=1,n
         Write(6,'(2(A,i2,A,ES20.10,A),2(2ES20.10,5x))')
     &                'W1(',i,')=',w1(i)-w1(1),'      ',
     &                'W2(',i,')=',w2(i)-w2(1),' diff=',
     &                 w1(i)-w2(i)
      End Do
!==================================================================

      Return
      End




