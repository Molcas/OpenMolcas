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
      subroutine termCF( ANGMOM, AMFI, ESFS, ldimcf, iDIM, maxes2, iopt,
     &                   nlanth, iprint )
c
c this Subroutine calculates the parameters of term-specific crystal field
c for lanthanides:
c
c  nstate                      -- total number of spin free states
c                                 mixed in RASSI
c  lDIMcf                      -- number of orbital states in the
c                                 considered term (LS term)
c  ANGMOM ( 3, nstate, nstate) -- a real array holding the matrix
c                                 elements of L in this basis
c                                 (read from RASSI):
c  ESFS (nstate)               -- energy of the orbital states form the LS term
c  maxes2(3,3)                 -- rotation matrix needed to choose the
c                                 main quantization axis
c                                 ( determinant(maxes)=1.0_wp, orthogonal vectors )
c  iopt                        -- option for choosing the main quanization axis
c                   iopt = 1   -- axis of the ground orbital multiplet,
c                                 iDIM specifies the size of pseudo L
c                   iopt = 2   -- axis of the entire L manifold
c                   iopt = 3   -- maxes is defined by the user
c                                 ( maxes2 is given as input )
c                   iopt = 4   -- maxes is the unity matrix ( original Z
c                                 is the quantization axis )
c
      Implicit None
      Integer, Parameter  :: wp=selected_real_kind(p=15,r=307)
#include "stdalloc.fh"
      Integer, intent(in)       :: ldimcf, iprint, iopt, nlanth, iDIM
      Real(kind=wp), intent(in) :: esfs(ldimcf)
      Real(kind=wp), intent(in) :: angmom(3,ldimcf,ldimcf)
      Real(kind=wp), intent(in) :: amfi(3,ldimcf,ldimcf)
      Real(kind=wp), intent(in) :: maxes2(3,3)

      Real(kind=wp)                 :: finddetr, dnrm2
!      Real(kind=wp)                 :: knm(12,0:12)
      Real(kind=wp), allocatable    :: maxes(:,:)
      Real(kind=wp), allocatable    :: gtens(:)
      Real(kind=wp), allocatable    :: eloc(:) ! lDIMcf
      Real(kind=wp), allocatable    :: Winit(:) ! lDIMcf

      Real(kind=wp)                 :: BNC(lDIMcf,0:lDIMcf)
      Real(kind=wp)                 :: BNS(lDIMcf,0:lDIMcf)
      Real(kind=wp)                 :: Bstev(lDIMcf,-lDIMcf:lDIMcf)
      Complex(kind=wp)              ::
     &                            Akq((lDIMcf-1),-(lDIMcf-1):(lDIMcf-1))
      Complex(kind=wp)              :: trace
      Complex(kind=wp), allocatable :: Angm(:,:,:) ! 3,ldimcf,ldimcf
      Complex(kind=wp), allocatable :: dipso(:,:,:) ! 3,ldimcf,ldimcf
      Complex(kind=wp), allocatable :: amfi_c(:,:,:) ! 3,ldimcf,ldimcf
      Complex(kind=wp), allocatable :: amfi2(:,:,:) ! 3,ldimcf,ldimcf
      Complex(kind=wp), allocatable :: amfi_l(:,:,:) ! 3,ldimcf,ldimcf
      Complex(kind=wp), allocatable :: Z(:,:) !ldimcf,ldimcf
      Complex(kind=wp), allocatable :: tmp(:,:) !ldimcf,ldimcf
      Complex(kind=wp), allocatable :: HCF(:,:) !ldimcf,ldimcf
      Complex(kind=wp), allocatable :: Zinit(:,:) !ldimcf,ldimcf

      Integer       :: i, j, l, info, i1, i2
      External      :: finddetr, trace, dnrm2
      Logical       :: debug =.false.
      Real(kind=wp) :: au2cm=2.194746313705d5

!
      Real(kind=wp) :: tS, tL, tJ, coeffCG, spinM, orbM,tJM,CF(100,100)
      Integer       :: MS, ML, MJ
      Integer       :: ij, iLS, nLS, ibasS(100), ibasL(100), ibasJ(100)
      Integer       :: irootL(100), ir, icas, k
      complex(kind=wp) :: CFC(100,100),zl(lDIMcf,lDIMcf)
      Call qEnter('TERMCF')
!============== End of variable declarations ==========================
      Call mma_allocate(maxes,3,3,'maxes')
      Call mma_allocate(gtens,3,'gtens')
      Call mma_allocate(eloc,lDIMcf,'eloc')
      Call mma_allocate(Winit,lDIMcf,'Winit')

      Call mma_allocate(Angm,3,ldimcf,ldimcf,'angm')
      Call mma_allocate(dipso,3,ldimcf,ldimcf,'dipso')
      Call mma_allocate(amfi_c,3,ldimcf,ldimcf,'amfi_c')
      Call mma_allocate(amfi2,3,ldimcf,ldimcf,'amfi2')
      Call mma_allocate(amfi_l,3,ldimcf,ldimcf,'amfi_l')
      Call mma_allocate(Z,ldimcf,ldimcf,'Z')
      Call mma_allocate(Zinit,ldimcf,ldimcf,'Zinit')
      Call mma_allocate(tmp,ldimcf,ldimcf,'tmp')
      Call mma_allocate(HCF,ldimcf,ldimcf,'HCF')

      Call dcopy_(3,[0.0_wp],0,gtens,1)
      Call dcopy_(3*3,[0.0_wp],0,maxes,1)
      Call dcopy_(lDIMcf,[0.0_wp],0,eloc,1)
      Call dcopy_(lDIMcf,[0.0_wp],0,Winit,1)

      Call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,Angm,1)
      Call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,dipso,1)
      Call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,amfi_c,1)
      Call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,amfi2,1)
      Call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,amfi_l,1)
      Call zcopy_(  ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,Z,1)
      Call zcopy_(  ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,Zinit,1)
      Call zcopy_(  ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,tmp,1)
      Call zcopy_(  ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,HCF,1)
!----------------------------------------------------------------------
      Write(6,'(/)')
      Write(6,'(100A)') ('%',i=1,95)
      If(MOD(lDIMcf,2)==1) Then
        Write(6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD '//
     &                         'PARAMETERS OF THE GROUND ATOMIC '//
     &                         'TERM, L = ', (lDIMcf-1)/2, '.'
      Else
        Write(6,'(5x,A,I2,A)') 'CALCULATION OF CRYSTAL-FIELD '//
     &                         'PARAMETERS OF THE GROUND ATOMIC '//
     &                         'TERM, L = ', (lDIMcf-1),'/2.'
      End If
      Write(6,'(100A)') ('%',i=1,95)
      Write(6,*)

      Do l=1,3
        Do i=1,ldimcf
          Do j=1,ldimcf
            dipso(l,i,j)=-CMPLX(0.0_wp,angmom(l,i,j),wp)
            amfi_c(l,i,j)=CMPLX(amfi(l,i,j),0.0_wp,wp)
          End Do
        End Do
      End Do
c  find the main anisotropy direction of the
c  diagonalize the angmom
      If ( iopt .eq. 1 ) Then
c                   iopt = 1   -- axis of the ground orbital Doublet
         Call atens( dipso(1:3,1:iDIM,1:iDIM), iDIM, gtens, maxes,
     &               iprint)
         Write(6,'(a)') 'The parameters of the Crystal Field matrix '//
     &                  'are written in the coordinate system:'
         If(MOD(iDIM,2) == 0) Then
           Write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic '//
     &                         'axes of the ground pseudo-L = |',
     &                          iDIM-1,'/2> orbital multiplet.'
         Else
           Write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic '//
     &                         'axes of the ground pseudo-L = |',
     &                         (iDIM-1)/2,'> orbital multiplet.'
         End If

      Else If ( iopt .eq. 2 ) Then
c                   iopt = 2   -- axis of the entire L manIfold
         Call atens( dipso(1:3,1:ldimcf,1:ldimcf), lDIMcf, gtens, maxes,
     &               iprint)
         If(MOD(lDIMCF,2).eq.0) Then
           Write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic '//
     &                         'axes of the ground atomic L = |',
     &                          lDIMCF-1,'/2> multiplet'
         Else
           Write(6,'(a,i2,a)') '(Xm, Ym, Zm) --  the main magnetic '//
     &                         'axes of the ground atomic L = |',
     &                         (lDIMCF-1)/2,'> multiplet'
         End If

      Else If ( iopt .eq. 3 ) Then
c                   iopt = 3   -- maxes is defined by the user ( maxes is given as input )
         Write(6,'(a)') 'The parameters of the Crystal Field matrix '//
     &                  'are written in the coordinate system:'
         Write(6,'(a)') '(Xm, Ym, Zm) -- defined in the input file.'
         If( finddetr(maxes2(1:3,1:3),3).eq.0.0_wp) Then
            Write(6,'(A)') 'TermCF:   iopt=3, while  DET(maxes)= 0.0'
            Call AbEnd()
         End If
         If( finddetr(maxes2(1:3,1:3),3).lt.0.0_wp) Then
           Do i=1,3
              maxes(i,1)=-maxes2(i,1)
           End Do
           If(iprint.gt.2) Write(6,'(a)')
     &                           'The original coordinate system '//
     &                           'was LEFT-handed. It has been '//
     &                           'changed to the RIGHT-handed'
         Else
           ! copy the maxes2 to maxes
           Do i=1,3
             Do j=1,3
               maxes(i,j)=maxes2(i,j)
             End Do
           End Do

         End If

      Else If ( iopt .eq. 4 ) Then
c                   iopt = 4   -- maxes is the unity matrix ( original Z is the quantization axis )
         Do i=1,3
            maxes(i,i)=1.0_wp
         End Do
      Else
         Call AbEnd()
      End If

      ! print out the rotation matrix:
      Write(6,'(a)') 'Rotation matrix from the initial coordinate '//
     &               'system to the employed coordinate system is:'

      If((iopt.eq.1).OR.(iopt.eq.2)) Then
         Write(6,'(70a)') ('-',i=1,67),'|'
         Write(6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
         Write(6,'(A,35x,A)') 'Xm, Ym, Zm -- main magnetic axes','|'
         Write(6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
         Write(6,'(6x,A,3F18.14,1x,A)')
     &                     '| Xm |',(maxes(j,1),j=1,3),'|'
         Write(6,'( A,A,3F18.14,1x,A)') ' R =  ',
     &                     '| Ym |',(maxes(j,2),j=1,3),'|'
         Write(6,'(6x,A,3F18.14,1x,A)')
     &                     '| Zm |',(maxes(j,3),j=1,3),'|'
         Write(6,'(83a)') ('-',i=1,67),'|'
         Write(6,'(A,I3)') 'Quantization axis is Zm.'

      Else If(iopt.eq.3) Then

         Write(6,'(70a)') ('-',i=1,67),'|'
         Write(6,'(A,31x,A)') 'x , y , z  -- initial Cartesian axes','|'
         Write(6,'(A,11x,A)') 'Xm, Ym, Zm -- the coordinate system '//
     &                        'defined in the input','|'
         Write(6,'(4x,3(17x,a),9x,a)') 'x','y','z','|'
         Write(6,'(6x,A,3F18.14,1x,A)')
     &                     '| Xm |',(maxes(j,1),j=1,3),'|'
         Write(6,'( A,A,3F18.14,1x,A)') ' R =  ',
     &                     '| Ym |',(maxes(j,2),j=1,3),'|'
         Write(6,'(6x,A,3F18.14,1x,A)')
     &                     '| Zm |',(maxes(j,3),j=1,3),'|'
         Write(6,'(83a)') ('-',i=1,67),'|'
         Write(6,'(A,I3)') 'Quantization axis is Zm.'

      Else

         Write(6,'(A)') 'IDENTITY matrix.'
         Write(6,'(A,I3)') 'Quantization axis is the initial z axis.'
      End If


!  rotate the angular momentum to the new axes, using "maxes"
      Call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,amfi2,1)
      Call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,angm,1)
      Call rotmom2( dipso(1:3,1:ldimcf,1:ldimcf), ldimcf,
     &              maxes, angm )
      Call rotmom2( amfi_c(1:3,1:ldimcf,1:ldimcf), ldimcf,
     &              maxes, amfi2 )
      If (debug) Call prmom('TERMCF:: ANGM',angm,ldimcf)
      If (debug) Call prmom('TERMCF:: AMFI2',amfi2,ldimcf)

!-----------------------------------------------------------------------
!  Find the pseudo-L basis of the LS manifold
      Call pseudospin(angm,ldimcf,Z,3,1,iprint)

      If ((iprint >= 4).or.debug) Then
        Write(6,*)
        Write(6,'(5X,A)') 'PSEUDO-L EIGENFUNCTIONS:'
        Write(6,*)
        If(MOD(lDIMcf,2).eq.1) Then
          Do I=1,lDIMcf
            Write(6,'(A,I3,A,3X,20(2F9.6,1X))') '|',
     &                 (lDIMcf-1)/2+(1-I),' > :',(Z(j,I),j=1,lDIMcf)
             do k=1,lDIMcf
                zl( (lDIMcf-1)/2+(1-I),k)=Z(k,I)
             enddo
          End Do
        Else
          Do I=1,lDIMcf
            Write(6,'(A,I3,A,3X,20(2F9.6,1X))') '|',
     &                 (lDIMcf-1)-2*(I-1),'/2 > :',(Z(j,I),j=1,lDIMcf)
          End Do
        End If
      End If
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  decompose the orbital moment AMSL in ITOs
C  transform AMFI integrals to pseudo-L basis:
      ! calculate the matrix elements of the spin and magnetic moment
      ! in the spin-orbit basis:
      Call zcopy_(3*ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,amfi_l,1)
      Do L=1,3
         Call zcopy_(  ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,tmp,1)
         ! amfi:
         Call ZGEMM_('C', 'N', lDIMcf, lDIMcf, lDIMcf, (1.0_wp,0.0_wp),
     &                     Z, lDIMcf,
     &          AMFI2(L,:,:), lDIMcf,           (0.0_wp,0.0_wp),
     &                   TMP, lDIMcf )
         Call ZGEMM_('N', 'N', lDIMcf, lDIMcf, lDIMcf, (1.0_wp,0.0_wp),
     &                    TMP, lDIMcf,
     &                      Z, lDIMcf,          (0.0_wp,0.0_wp),
     &          AMFI_L(L,:,:), lDIMcf )

      End Do !L

      If (debug) Call prMom_herm('TERMCF:: AMFI_L',amfi_l*au2cm,ldimcf)


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c  calculate the RASSI Crystal Field matrix
      Call dcopy_(lDIMcf,[0.0_wp],0,eloc,1)
      Call rtrace(lDIMcf,ESFS,ELOC)
      Call zcopy_(ldimcf*ldimcf,[(0.0_wp,0.0_wp)],0,HCF,1)
      Do i=1,lDIMcf
        Do I1=1,lDIMcf
          Do I2=1,lDIMcf
            HCF(I1,I2)=HCF(I1,I2)+ELOC(i)*conjg(Z(i,I1))*Z(i,I2)
          End Do
        End Do
      End Do

      info=0
      Call dcopy_(lDIMcf,[0.0_wp],0,Winit,1)
      Call zcopy_(lDIMcf*lDIMcf,[(0.0_wp,0.0_wp)],0,Zinit,1)
      Call DIAG_C2( HCF,lDIMcf,info,Winit,Zinit)
      Call print_ZFS('Ab Initio Calculated Crystal-Field Splitting '//
     &               'Matrix written in the basis of Pseudo-L '//
     &               'Eigenfunctions',HCF,lDIMCF)

      Call NEWCF(HCF,lDIMcf,Akq,BNC,BNS,BStev)
!ccccccccccc   print CF parameter ccccccccccccccccccccccccccccccc
      If((iprint>=4).or.debug) Then
         Call print_CFP_LCLU(lDIMCF,BNC,BNS,.true.)
         Call print_CFP_stev(lDIMCF,Bstev,.true.)
         Call print_CFP_naoya(lDIMcf,Akq,.true.)
      Else
         Call print_CFP_LCLU(lDIMCF,BNC,BNS,.false.)
         Call print_CFP_stev(lDIMCF,Bstev,.false.)
         Call print_CFP_naoya(lDIMcf,Akq,.false.)
      End If

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Write(6,'(/)')
      If(MOD(lDIMCF,2)==1) Then
         Write(6,'(A,I0)') 'DECOMPOSITION OF THE RASSCF WAVE '//
     &                     'FUNCTIONS CORRESPONDING TO THE '//
     &                     'LOWEST ATOMIC MULTIPLET L =',
     &                      (lDIMCF-1)/2
      Else
         Write(6,'(A,I0,A)') 'DECOMPOSITION OF THE RASSCF WAVE '//
     &                       'FUNCTIONS CORRESPONDING TO THE '//
     &                       'LOWEST ATOMIC MULTIPLET L = ',
     &                       (lDIMCF-1), '/2'
      End If
         Write(6,'(A,I0)') 'IN WAVE FUNCTIONS WITH DEFINITE '//
     &                     'PROJECTION OF THE TOTAL MOMENT '//
     &                     'ON THE QUANTIZATION AXIS'
      Call print_ZFS_naoya('L',Zinit,lDIMcf)
      Call individual_ranks(lDIMCF,BNC,BNS,HCF,'L',iprint)

!  saving some information for tests:
      Call Add_Info('CRYS_TERM_BNMC_20',BNC(2,0),1,4)
      Call Add_Info('CRYS_TERM_BNMC_40',BNC(4,0),1,4)
      Call Add_Info('CRYS_TERM_BNMC_60',BNC(6,0),1,4)
      goto 999
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      generate |J,MJ> states using the spin |S,MS> and |L,ML> states
c
      tS=0.d0
      tL=0.d0
      tJ=0.d0
      MS=0
      ML=0
      MJ=0
      If (nlanth== 1)      then  ! Ce  J=5/2
          tS=0.5d0
          tJ=2.5d0
      else if (nlanth== 2) then  ! Pr  J=4
          tS=1.d0
          tJ=4.d0
      else if (nlanth== 3) then ! Nd  J=9/2
          tS=1.5d0
          tJ=4.5d0
      else if (nlanth== 4) then ! Pm  J=4
          tS=2.0d0
          tJ=4.d0
      else if (nlanth== 5) then ! Sm  J=5/2
          tS=2.5d0
          tJ=2.5d0
      else if (nlanth== 6) then ! Eu  J=0
          tS=3.0d0
          tJ=0.0d0
      else if (nlanth== 7) then ! Gd  J=7/2; S=7/2
          tS=3.5d0
          tJ=3.5d0
      else if (nlanth== 8) then ! Tb  J=6
          tS=3.0d0
          tJ=6.0d0
      else if (nlanth== 9) then ! Dy  J=15/2
          tS=2.5d0
          tJ=7.5d0
      else if (nlanth==10) then ! Ho  J=8
          tS=2.0d0
          tJ=8.0d0
      else if (nlanth==11) then ! Er  J=15/2
          tS=1.5d0
          tJ=7.5d0
      else if (nlanth==12) then ! Tm  J=6
          tS=1.0d0
          tJ=6.0d0
      else if (nlanth==13) then ! Yb  J=7/2
          tS=0.5d0
          tJ=3.5d0
      else
          write (6,'(A)') 'not implemented yet'
      endif
      tL =(dble(lDIMcf)-1.d0)/2.d0
      ML=lDIMcf
      MS=nint(2.d0*tS+1.d0)
      MJ=nint(2.d0*tJ+1.d0)

c      print *,'nlanth=', nlanth, '  S,  L,  J =', tS, tL, tJ,
c     &                           ' MS, ML, MJ =', MS, ML, MJ

c      print *, 'build coupled basis |L,ML>|S,MS>'
      nLS=ML*MS
      ij=0
      ir=0
      do i=-ML+1,ML-1,2
         ir=ir+1
         do j=-MS+1,MS-1,2
            ij=ij+1
            ibasL(ij)=i
            ibasS(ij)=j
            irootL(ij)=ir
         enddo
      enddo
      do iLS=1,nLS
         write(6,'(A,3I4,2x,A,2F6.1)')
     &        'nLS,  ML,  MS =', iLS, ibasL(iLS), ibasS(iLS),
     &        'ML, MS =', dble(ibasL(iLS))/2.d0, dble(ibasS(iLS))/2.d0
      enddo

      write(6,*) 'proceed to build |J,MJ> states'
      ij=0
      do i=-MJ+1,MJ-1,2
        ij=ij+1
        ibasJ(ij)=i
      enddo

      do i=1,MJ
         write(6,'(i3,i6,F6.1)') i, ibasJ(i), dble(ibasJ(i))/2.d0
      enddo

      Cf=0.0_wp
      CfC=(0.0_wp,0.0_wp)

      do ij=1,MJ
         tJM=dble(ibasJ(ij))/2.d0
         do iLS=1,nLS
              ! set projections
              spinM=dble(ibasS(iLS))/2.d0
              orbM=dble(ibasL(iLS))/2.d0

              Call Clebsh_Gordan(tL,orbM,tS,spinM,tJ,tJM, coeffCG)
              Cf(iJ,iLS)=coeffCG

              If(abs(coeffCG)>1.d-20) write(6,*) 'ij,iLS,coeffCG',
     *                                            ij,iLS,coeffCG
         enddo
      enddo

      Write(6,*) 'MJ ->  (1,16), (2,15), (3,14), (4,13), '//
     *                   '(5,12), (6,11), (7,10), (8,9)'
      do iLS=1,nLS
        Write(6,'(A,2F6.1,16F11.8)') 'ML,MS: ',dble(ibasL(iLS))/2.d0,
     *   dble(ibasS(iLS))/2.d0,
     &   Cf(1,iLS), Cf(MJ,iLS),
     &   Cf(2,iLS), Cf(MJ-1,iLS),
     &   Cf(3,iLS), Cf(MJ-2,iLS),
     &   Cf(4,iLS), Cf(MJ-3,iLS),
     &   Cf(5,iLS), Cf(MJ-4,iLS),
     &   Cf(6,iLS), Cf(MJ-5,iLS),
     &   Cf(7,iLS), Cf(MJ-6,iLS),
     &   Cf(8,iLS), Cf(MJ-7,iLS)
      enddo


      Write(6,*) 're-write initial CASSCF states into |J,MJ>, using '
     &        //'( Z(j,I),j=1,lDIMcf)  coefficients'



      Do ij=1,MJ
        Do iLS=1,nLS

          k=ibasL(iLS)/2

           do iCAS=1,ML
              CFC(iJ,iLS) = CFC(iJ,iLS) + Z(k,iCAS)*Cf(iJ,iLS)
           enddo


        End Do
      End Do


      Write(6,*)  'MJ ->  (1,16), (2,15), (3,14), (4,13), '//
     *                   '(5,12), (6,11), (7,10), (8,9)'
      do iLS=1,nLS
         write(6,'(A,i2,F6.1,16(2F8.4,2x))') 'iCAS,MS: ',
     &           irootL(iLS), dble(ibasS(iLS))/2.d0,
     &   CfC(1,iLS), CfC(MJ,  iLS),
     &   CfC(2,iLS), CfC(MJ-1,iLS),
     &   CfC(3,iLS), CfC(MJ-2,iLS),
     &   CfC(4,iLS), CfC(MJ-3,iLS)
c     &   CfC(5,iLS), CfC(MJ-4,iLS)
c     &   CfC(6,iLS), CfC(MJ-5,iLS),
c     &   CfC(7,iLS), CfC(MJ-6,iLS),
c     &   CfC(8,iLS), CfC(MJ-7,iLS)
      enddo


 999  Continue
!-----------------------------------------------------------------------
      Call mma_deallocate(maxes)
      Call mma_deallocate(gtens)
      Call mma_deallocate(eloc)
      Call mma_deallocate(Winit)

      Call mma_deallocate(Angm)
      Call mma_deallocate(dipso)
      Call mma_deallocate(amfi_c)
      Call mma_deallocate(amfi2)
      Call mma_deallocate(amfi_l)
      Call mma_deallocate(Z)
      Call mma_deallocate(Zinit)
      Call mma_deallocate(tmp)
      Call mma_deallocate(HCF)

      Call qExit('TERMCF')
      Return
      End
