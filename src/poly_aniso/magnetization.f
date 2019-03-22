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
      Subroutine magnetization( exch, nLoc, nM, nH, nneq, neq, neqv,
     &                          nCenter, nTempMagn, nDir, nDirZee,
     &                          nDirTot, nss, nexch, iopt, LUZee,
     &                          TempMagn, hexp, mexp, hmin, hmax, em,
     &                          zJ, thrs, dirX, dirY, dirZ, dir_weight,
     &                          w, dipexch, s_exch, dipso, s_so, eso,
     &                          hinput, r_rot, XLM, ZLM, XRM, ZRM,
     &                          zeeman_energy, compute_Mdir_vector,
     &                          m_paranoid, m_accurate, smagn, mem )

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "mgrid.fh"
#include "stdalloc.fh"
c constants defining the sizes
      Integer, intent(in)       :: exch, nLoc, nM
      Integer, intent(in)       :: nH, nCenter, nTempMagn
      Integer, intent(in)       :: nDir, nDirZee, nDirTot, nneq, neqv
      Integer, intent(in)       :: neq(nneq), nss(nneq), nexch(nneq)
      Integer, intent(in)       :: iopt, mem
      Integer, intent(in)       :: LUZee(nDirZee)

      Logical, intent(in)       :: hinput
      Logical, intent(in)       :: zeeman_energy
      Logical, intent(in)       :: compute_Mdir_vector
      Logical, intent(in)       :: m_paranoid
      Logical, intent(in)       :: m_accurate
      Logical, intent(in)       :: smagn

      Real(kind=wp), intent(in) :: R_ROT(nneq,neqv,3,3)
      Real(kind=wp), intent(in) :: W(exch)
 ! exchange energies printed out in the previous part
      Real(kind=wp), intent(in) :: ESO(nneq,nLoc)
! spin-orbit energies from ANISO files
      Real(kind=wp), intent(in) :: Hexp(nH), Mexp(nH,nTempMagn)
      Real(kind=wp), intent(in) :: thrs
      Real(kind=wp), intent(in) :: XLM( nCenter,nTempMagn,3,3)
      Real(kind=wp), intent(in) :: ZLM( nCenter,nTempMagn)
      Real(kind=wp), intent(in) :: XRM( nCenter,nTempMagn,3,3)
      Real(kind=wp), intent(in) :: ZRM( nCenter,nTempMagn)
      Real(kind=wp), intent(in) :: dirX(nDir), dirY(nDir), dirZ(nDir)
      Real(kind=wp), intent(in) :: dir_weight(nDirZee,3)
      Real(kind=wp), intent(in) :: TempMagn(nTempMagn)
      Real(kind=wp), intent(in) :: zJ, hmin, hmax, em

      Complex(kind=wp), intent(in) :: DIPEXCH(3,exch,exch)
      Complex(kind=wp), intent(in) ::  S_EXCH(3,exch,exch)
      Complex(kind=wp), intent(in) :: dipso(nneq,3,nLoc,nLoc)
      Complex(kind=wp), intent(in) ::  s_so(nneq,3,nLoc,nLoc)
c exchange data:
c      Integer NM ! number of states included in the exchange Zeeman matrix, ( Nex .LE. exch)
      Real(kind=wp), allocatable :: Wex(:)
!                                   WEX(NM) ! Zeeman exchange energies
      Real(kind=wp), allocatable :: Zex(:)
!                                   ZEX(nTempMagn) ! exchange statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: Sex(:,:)
!                                   SEX(3,nTempMagn) ! spin magnetisation, from the exchange block;
      Real(kind=wp), allocatable :: Mex(:,:)
!                                   MEX(3,nTempMagn) ! magnetisation, from the exchange block
c data for individual sites (all states):
      Real(kind=wp), allocatable :: ZL(:,:)
!                                   ZL(nneq,nTempMagn) ! local statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: WL(:,:)
!                                   WL(nneq,nLoc) ! Zeeman local energis
      Real(kind=wp), allocatable :: SL(:,:,:)
!                                   SL(nneq,3,nTempMagn) ! spin magnetisation, from the local sites, using ALL states ;
      Real(kind=wp), allocatable :: ML(:,:,:)
!                                   ML(nneq,3,nTempMagn) ! magnetisation, from local sites, using ALL states;
c data for individual sites (only states that enter exchange):
      Real(kind=wp), allocatable :: ZR(:,:)
!                                   ZR(nneq,nTempMagn) ! local statistical sum, Boltzmann distribution, using only Nexch states
      Real(kind=wp), allocatable :: WR(:,:)
!                                   WR(nneq,nLoc) ! Zeeman local reduced energies, using only Nexch states;
      Real(kind=wp), allocatable :: SR(:,:,:)
!                                   SR(nneq,3,nTempMagn) ! spin magnetisation, from the local sites, using only Nexch states ;
      Real(kind=wp), allocatable :: MR(:,:,:)
!                                   MR(nneq,3,nTempMagn) ! magnetisation, from local sites, using only Nexch states;
c total vectors in general coordinate system:
      Real(kind=wp), allocatable :: ZRT(:,:)  !ZRT(nCenter,nTempMagn)
      Real(kind=wp), allocatable :: ZLT(:,:)  !ZLT(nCenter,nTempMagn)
      Real(kind=wp), allocatable :: MRT(:,:,:)!MRT(nCenter,3,nTempMagn)
      Real(kind=wp), allocatable :: MLT(:,:,:)!MLT(nCenter,3,nTempMagn)
      Real(kind=wp), allocatable :: SRT(:,:,:)!SRT(nCenter,3,nTempMagn)
      Real(kind=wp), allocatable :: SLT(:,:,:)!SLT(nCenter,3,nTempMagn)
c data for total system:
      Real(kind=wp), allocatable :: ZT(:,:)
!                                   ZT(nH,nTempMagn) ! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ST(:,:,:)
!                                   ST(3,nH,nTempMagn) ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT(:,:,:)
!                                   MT(3,nH,nTempMagn) ! total magnetisation
c magnetic field strength and orientation data:
      Real(kind=wp) :: dltH
      Real(kind=wp), allocatable :: H(:) ! H(nH)
      Real(kind=wp), allocatable :: dHX(:) !dHX(nDirTot)
      Real(kind=wp), allocatable :: dHY(:) !dHY(nDirTot)
      Real(kind=wp), allocatable :: dHZ(:) !dHZ(nDirTot)
      Real(kind=wp), allocatable :: dHW(:) !dHW(nDirTot)
c total average M and average S data:
      Real(kind=wp), allocatable :: MAV(:,:)      !MAV(nH,nTempMagn)
      Real(kind=wp), allocatable :: SAV(:,:)      !SAV(nH,nTempMagn)
      Real(kind=wp), allocatable :: MVEC(:,:,:,:)
!MVEC(nDirTot,nH,nTempMagn,3)
      Real(kind=wp), allocatable :: SVEC(:,:,:,:)
!SVEC(nDirTot,nH,nTempMagn,3)

      Integer                   :: IM,I,it,itEnd,J,iH,k,isite,l,n,nP
      Integer                   :: iDir,rtob,ibuf,mem_local
      Real(kind=wp)             :: cm3tomB
      Real(kind=wp)             :: dev, dnrm2_
      External                  :: dev, dnrm2_
      Logical          :: DBG
      Character(15)    :: lbl_X, lbl_Y, lbl_Z
c

      Call qEnter('magnetisation')

      DBG=.false.
c      Boltz_k=0.6950356000_wp   !   in cm^-1*K-1
c      mu_Bohr=0.4668643740_wp   !   in cm-1*T-1
      cm3tomB=0.5584938904_wp   !   in cm3 * mol-1 * T

      Write(6,*)
      Write(6,'(100A)') (('%'),J=1,96)
      Write(6,'(40X,A)') 'CALCULATION OF THE MOLAR MAGNETIZATION'
      Write(6,'(100A)') (('%'),J=1,96)
      Write(6,*)
c----------------------------------------------------------------------
      mem_local=0
      RtoB=8
      If(dbg) Write(6,*) 'MAGN:        nM=',nM
      If(dbg) Write(6,*) 'MAGN:      exch=',exch
      If(dbg) Write(6,*) 'MAGN:      nLoc=',nLoc
      If(dbg) Write(6,*) 'MAGN: nTempMagn=',nTempMagn

      If(nM>0) Then
        ! Zeeman exchange energy spectrum
        Call mma_allocate(Wex,nM,'Wex')
        Call dcopy_(nM,0.0_wp,0,Wex,1)
        mem_local=mem_local+nM*RtoB
      End If

      If(nTempMagn>0) Then
        ! exchange statistical sum, Boltzmann distribution
        Call mma_allocate(Zex,nTempMagn,'Zex')
        Call dcopy_(nTempMagn,0.0_wp,0,Zex,1)
        mem_local=mem_local+nTempMagn*RtoB
        ! spin magnetisation, from the exchange block
        Call mma_allocate(Sex,3,nTempMagn,'Sex')
        Call dcopy_(3*nTempMagn,0.0_wp,0,Sex,1)
        mem_local=mem_local+3*nTempMagn*RtoB
        ! magnetisation, from the exchange block
        Call mma_allocate(Mex,3,nTempMagn,'Mex')
        Call dcopy_(3*nTempMagn,0.0_wp,0,Mex,1)
        mem_local=mem_local+3*nTempMagn*RtoB

        If(nneq>0) Then
          ! local statistical sum, Boltzmann distribution
          Call mma_allocate(ZL,nneq,nTempMagn,'ZL')
          Call dcopy_(nneq*nTempMagn,0.0_wp,0,ZL,1)
          mem_local=mem_local+nneq*nTempMagn*RtoB
          ! spin magnetisation, from the local sites, using ALL states
          Call mma_allocate(SL,nneq,3,nTempMagn,'SL')
          Call dcopy_(3*nneq*nTempMagn,0.0_wp,0,SL,1)
          mem_local=mem_local+3*nneq*nTempMagn*RtoB
          !magnetisation, from local sites, using ALL states
          Call mma_allocate(ML,nneq,3,nTempMagn,'ML')
          Call dcopy_(3*nneq*nTempMagn,0.0_wp,0,ML,1)
          mem_local=mem_local+3*nneq*nTempMagn*RtoB

!         local statistical sum, Boltzmann distribution, using only Nexch states
          Call mma_allocate(ZR,nneq,nTempMagn,'ZR')
          Call dcopy_(nneq*nTempMagn,0.0_wp,0,ZR,1)
          mem_local=mem_local+nneq*nTempMagn*RtoB
!         spin magnetisation, from the local sites, using only Nexch states
          Call mma_allocate(SR,nneq,3,nTempMagn,'SR')
          Call dcopy_(3*nneq*nTempMagn,0.0_wp,0,SR,1)
          mem_local=mem_local+3*nneq*nTempMagn*RtoB
          ! magnetisation, from local sites, using only Nexch states
          Call mma_allocate(MR,nneq,3,nTempMagn,'MR')
          Call dcopy_(3*nneq*nTempMagn,0.0_wp,0,MR,1)
          mem_local=mem_local+3*nneq*nTempMagn*RtoB
        End If

        If(nLoc>0) Then
          ! Zeeman local energies
          Call mma_allocate(WL,nneq,nLoc,'WL')
          Call dcopy_(nneq*nLoc,0.0_wp,0,WL,1)
          mem_local=mem_local+nneq*nLoc*RtoB
          ! Zeeman local reduced energies, using only Nexch states
          Call mma_allocate(WR,nneq,nLoc,'WR')
          Call dcopy_(nneq*nLoc,0.0_wp,0,WR,1)
          mem_local=mem_local+nneq*nLoc*RtoB
        End If

        If(nCenter>0) Then
          ! ZRT(nCenter,nTempMagn)
          Call mma_allocate(ZRT,nCenter,nTempMagn,'ZRT')
          Call dcopy_(nCenter*nTempMagn,0.0_wp,0,ZRT,1)
          mem_local=mem_local+nCenter*nTempMagn*RtoB
          ! ZLT(nCenter,nTempMagn)
          Call mma_allocate(ZLT,nCenter,nTempMagn,'ZLT')
          Call dcopy_(nCenter*nTempMagn,0.0_wp,0,ZLT,1)
          mem_local=mem_local+nCenter*nTempMagn*RtoB
          ! MRT(nCenter,3,nTempMagn)
          Call mma_allocate(MRT,nCenter,3,nTempMagn,'MRT')
          Call dcopy_(nCenter*3*nTempMagn,0.0_wp,0,MRT,1)
          mem_local=mem_local+3*nCenter*nTempMagn*RtoB
          ! MLT(nCenter,3,nTempMagn)
          Call mma_allocate(MLT,nCenter,3,nTempMagn,'MLT')
          Call dcopy_(nCenter*3*nTempMagn,0.0_wp,0,MLT,1)
          mem_local=mem_local+3*nCenter*nTempMagn*RtoB
          ! SRT(nCenter,3,nTempMagn)
          Call mma_allocate(SRT,nCenter,3,nTempMagn,'SRT')
          Call dcopy_(nCenter*3*nTempMagn,0.0_wp,0,SRT,1)
          mem_local=mem_local+3*nCenter*nTempMagn*RtoB
          ! SLT(nCenter,3,nTempMagn)
          Call mma_allocate(SLT,nCenter,3,nTempMagn,'SLT')
          Call dcopy_(nCenter*3*nTempMagn,0.0_wp,0,SLT,1)
          mem_local=mem_local+3*nCenter*nTempMagn*RtoB
        End If

        If(nH>0) Then
          ! total statistical sum, Boltzmann distribution
          Call mma_allocate(ZT,nH,nTempMagn,'ZT')
          Call dcopy_(nH*nTempMagn,0.0_wp,0,ZT,1)
          mem_local=mem_local+nH*nTempMagn*RtoB
          ! total spin magnetisation
          Call mma_allocate(ST,3,nH,nTempMagn,'ST')
          Call dcopy_(3*nH*nTempMagn,0.0_wp,0,ST,1)
          mem_local=mem_local+3*nH*nTempMagn*RtoB
          ! total magnetisation
          Call mma_allocate(MT,3,nH,nTempMagn,'MT')
          Call dcopy_(3*nH*nTempMagn,0.0_wp,0,MT,1)
          mem_local=mem_local+3*nH*nTempMagn*RtoB
          ! total spin magnetisation
          Call mma_allocate(SAV,nH,nTempMagn,'SAV')
          Call dcopy_(nH*nTempMagn,0.0_wp,0,SAV,1)
          mem_local=mem_local+nH*nTempMagn*RtoB
          ! total magnetisation
          Call mma_allocate(MAV,nH,nTempMagn,'MAV')
          Call dcopy_(nH*nTempMagn,0.0_wp,0,MAV,1)
          mem_local=mem_local+nH*nTempMagn*RtoB
          If(nDirTot>0) Then
            ! total spin magnetisation vector
            Call mma_allocate(SVEC,nDirTot,nH,nTempMagn,3,'SVEC')
            Call dcopy_(nDirTot*nH*nTempMagn*3,0.0_wp,0,SVEC,1)
            mem_local=mem_local+nDirTot*nH*nTempMagn*3*RtoB
            ! total magnetisation vector
            Call mma_allocate(MVEC,nDirTot,nH,nTempMagn,3,'MVEC')
            Call dcopy_(nDirTot*nH*nTempMagn*3,0.0_wp,0,MVEC,1)
            mem_local=mem_local+nDirTot*nH*nTempMagn*3*RtoB
          End If

        End If
      End If

      If(nDirTot>0) Then
         ! orientation of the field
         Call mma_allocate(dHX,nDirTot,'dHX')
         Call mma_allocate(dHY,nDirTot,'dHY')
         Call mma_allocate(dHZ,nDirTot,'dHZ')
         Call mma_allocate(dHW,nDirTot,'dHW')
         Call dcopy_(nDirTot,0.0_wp,0,dHX,1)
         Call dcopy_(nDirTot,0.0_wp,0,dHY,1)
         Call dcopy_(nDirTot,0.0_wp,0,dHZ,1)
         Call dcopy_(nDirTot,0.0_wp,0,dHW,1)
         mem_local=mem_local+4*nDirTot*RtoB
      End If

      If(nH>0) Then
         Call mma_allocate(H,nH,'H  field')
         Call dcopy_(nH,0.0_wp,0,H,1)
         mem_local=mem_local+nH*RtoB
      End If
      If(dbg) Write(6,*) 'MAGN:  memory allocated (local):'
      If(dbg) Write(6,*) 'mem_local=', mem_local
      If(dbg) Write(6,*) 'MAGN:  memory allocated (total):'
      If(dbg) Write(6,*) 'mem_total=', mem+mem_local
c----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      nP=get_nP(nsymm,ngrid)

      Call hdir( nDir, nDirZee, dirX, dirY, dirZ, dir_weight,
     &                 nP, nsymm, ngrid, nDirTot, dHX, dHY, dHZ, dHW )
      Call Add_Info('MR_MAGN  dHX', dHX,nDirTot,10)
      Call Add_Info('MR_MAGN  dHY', dHY,nDirTot,10)
      Call Add_Info('MR_MAGN  dHZ', dHZ,nDirTot,10)
      Call Add_Info('MR_MAGN  dWX', dHW,nDirTot,10)
      If(DBG) Then
        Write(6,'(A,F10.6)') '        zJ = ',zJ
        Write(6,'(A,F10.6)') '      thrs = ',thrs
        Write(6,'(A,   I6)') '      iopt = ',iopt
        Write(6,'(A, 10I6)') '     nss(i)=',(  nss(i),i=1,nneq)
        Write(6,'(A, 10I6)') '   nexch(i)=',(nexch(i),i=1,nneq)
        Write(6,'(A,   I6)') '  nTempMagn= ',nTempMagn
        Write(6,'(A      )') 'TempMagn(i)='
        Write(6,'(10F10.6)') (TempMagn(i),i=1,nTempMagn)
        Write(6,*) 'm_paranoid =', m_paranoid
        Write(6,*) 'm_accurate =', m_accurate
        Write(6,*) '     smagn =', smagn
        Do k=1,nneq
          Write(6,'(A,i6)') 'site ', k
          Write(6,'(A,i2,A)') 'ESO(',k,')'
          Write(6,'(10F12.6)') (ESO(k,i),i=1,nss(k))
          Write(6,'(A,i2,A)') 'S_SO(',k,')'
          Do i=1,nss(k)
            Do j=1,nss(k)
              Write(6,'(3(2F15.10,3x))') (s_so(k,l,i,j),l=1,3)
            End Do
          End Do
          Write(6,'(A,i2,A)') 'DIPSO(',k,')'
          Do i=1,nss(k)
            Do j=1,nss(k)
              Write(6,'(3(2F15.10,3x))') (dipso(k,l,i,j),l=1,3)
            End Do
          End Do
        End Do !k
      End If

      Write(6,'(2X,A,i3,A)') 'Molar magnetization will be '//
     &                       'calculated in ',nH,
     &                       ' points, equally distributed in '//
     &                       'magnetic field range'
      Write(6,'(2X,F4.1,1x,a,1X,F4.1,a,5(F6.3,a))') HMIN,'--',HMAX,
     &                       ' T., at the following temperatures:'
      Do i=1,nTempMagn,10
        j=MIN(nTempMagn,i+9)
        Write(6,'(10(F8.4,A))') (TempMagn(k),' K.;',k=i,j)
      End Do
      Write(6,'(2X,A,I4,A)') 'Powder molar magnetization will be '//
     &                       'averaged on ',nP,
     &                       ' directions of the applied magnetic'//
     &                       ' field.'
      Write(6,'(2x,10A)') ('--------', i=1,10)
      If(nsymm.eq.1) Then
        Write(6,'(23x,A)') 'Lebedev-Laikov grid on a hemisphere:'
        Write(6,'(38x,A)') 'z >= 0;'
      Else If(nsymm.eq.2) Then
        Write(6,'(23x,A)') 'Lebedev-Laikov grid on a 4th-of-a-sphere:'
        Write(6,'(34x,A)') 'x >= 0; z >= 0;'
      Else If(nsymm.eq.3) Then
        Write(6,'(23x,A)') 'Lebedev-Laikov grid on a 8th-of-a-sphere:'
        Write(6,'(30x,A)') 'x >= 0; y >= 0; z >= 0;'
      End If
      Write(6,'(2x,10A)') ('--------', i=1,10)
      Write(6,'(2x,A,12x,A,2(18x,A),16x,A)') 'Nr.','x','y','z',
     &                                       'weight'
      Do i=1,nP
        Write(6,'(i4,2x,4(F18.12,1x))')
     &                  i,dHX(i+nDir+nDirZee),dHY(i+nDir+nDirZee),
     &                    dHZ(i+nDir+nDirZee),dHW(i+nDir+nDirZee)
      End Do
      Write(6,'(2x,10A)') ('--------', i=1,10)
      Write(6,'(2X,A)') 'The cut-off energy for the '//
     &                  'exact diagonalization of the Zeeman'//
     &                  ' Hamiltonian is:'
      Write(6,'(2x,a,F15.9,A)') 'E = ',EM ,' cm(-1).'
      If(NM.lt.10) Then
        Write(6,'(2X,A,i2,a)') 'The exact diagonalization of the '//
     &                         'Zeeman Hamiltonian included ',NM,
     &                         ' exchange states.'
      Else If( (NM.ge.10) .AND. (NM.lt.100) ) Then
        Write(6,'(2X,A,i3,a)') 'The exact diagonalization of the '//
     &                         'Zeeman Hamiltonian included ',NM,
     &                         ' exchange states.'
      Else If( (NM.ge.100) .AND. (NM.lt.1000) ) Then
        Write(6,'(2X,A,i4,a)') 'The exact diagonalization of the '//
     &                         'Zeeman Hamiltonian included ',NM,
     &                         ' exchange states.'
      Else If( (NM.ge.1000) .AND. (NM.lt.10000) ) Then
        Write(6,'(2X,A,i5,a)') 'The exact diagonalization of the '//
     &                         'Zeeman Hamiltonian included ',NM,
     &                         ' exchange states.'
      End If
      If(m_accurate) Then
        Write(6,'(2x,A)') 'The contribution of local excited states '//
     &                    ' is computed exactly.'
        If ( (zJ.ne.0.0_wp) .and. (m_paranoid) .and.
     &  (nTempMagn.gt.1) ) Then
          Write(6,'(2x,A)') 'The average spin is computed exactly '//
     &                      'for each temperature point.'
        Else If( (zJ.ne.0.0_wp) .and. (.not. m_paranoid) .and.
     &   (nTempMagn.gt.1) ) Then
          Write(6,'(2x,A,F9.3)') 'The average spin is computed '//
     &                           'exactly only for the temperature '//
     &                           'point T= ',MAXVAL(TempMagn(:))
          Write(6,'(2x,A     )') 'We consider this to be a good '//
     &                           'approximation.'
        End If
      Else
        If ( (zJ.ne.0.0_wp) .and. (m_paranoid) .and.
     &  (nTempMagn.gt.1) ) Then
          Write(6,'(2x,A)') 'The average spin is computed exactly '//
     &                      'for each temperature point.'
        Else If( (zJ.ne.0.0_wp) .and. (.not. m_paranoid) .and.
     &   (nTempMagn.gt.1) ) Then
          Write(6,'(2x,A,F9.3)') 'The average spin is computed '//
     &                           'exactly only for the temperature '//
     &                           ' point T= ',TempMagn(1)
          Write(6,'(2x,A)') 'We consider this to be a good '//
     &                      'approximation.'
          Write(6,'(2x,A)') 'Use MPAR keyword to compute average '//
     &                      'spin exactly, for each temperature point'//
     &                      ', in case higher accuracy is needed.'
        End If
      Write(6,'(2x,A)') 'The contribution of local excited states '//
     & ' is computed approximately (by using the susceptibility data). '
      End If
      If(compute_Mdir_vector) Then
      Write(6,'(2x,A,i2,a)') 'The magnetization vector for ',nDir,
     & ' directions of the applied magnetic field will be calculated.'
      Else
      Write(6,'(2X,A)') 'The magnetization vector was not calculated.'
      End If
      If(zeeman_energy) Then
      Write(6,'(2x,A,i2,a)') 'The Zeeman splitting for ',nDirZee,
     & ' directions of the applied magnetic field will be calculated.'
      Write(6,'(2x,a)') 'The Zeeman energies for each direction of th'//
     & 'e applied magnetic field are written in files "energy_XX.txt".'
         If(DBG) Then
            Do i=1,nDirZee
              Write(6,'(A,I4,A,I4)') 'LuZee(',i,' )=',LUZee(i)
            End Do
         End If
      Else
      Write(6,'(2X,A)') 'Computation of the Zeeman splitting was not'//
     & ' requested.'
      End If
c      If(maxval(TempMagn).lt.W(exch)) Then
c      Write(6,'(2x,A)') 'Contribution to molar magnetization coming '//
c     & 'from local excited states is taken into account.'
c      Else
c      Write(6,'(2x,A)') 'Contribution to molar magnetization coming '//
c     & 'from local excited states is NOT taken into account.'
c      Write(6,'(2x,A)') 'TMAG is requesting to compute magnetization'//
c     & ' at a larger temperature than the highest exchange state.'
c      Write(6,'(2x,A)') 'Please include more states into the '//
c     & 'exchange coupling.'
c      End If

C /// opening the loop over the field points
      Do iH=1,nH
c /// ---------------------------------------------------------------
        If (HINPUT) Then
           H(iH)=HEXP(iH)
           If(H(iH).eq.0.0_wp) Then
             H(iH)=0.0001_wp
           End If
        Else
          DLTH=(HMAX-HMIN)/DBLE(nH-1)
          If (iH.eq.1) Then
            H(iH)=HMIN+0.0001_wp
          Else
            H(iH)=HMIN+DLTH*DBLE(iH-1)
          End If
          If(H(iH).eq.0.0_wp) Then
            H(iH)=0.0001_wp
          End If
        End If

        If(DBG) Write(6,'(A,i0,A,F10.5)') 'MAGNETIZATION::  H(',
     &                                     iH,') = ', H(iH)

c ///  opening the loop over different directions of the magnetic field
        Do IM=1,NDIRTOT
          Call dcopy_(nM,0.0_wp,0,Wex,1)
          Call dcopy_(nneq*nLoc,0.0_wp,0,WL,1)
          Call dcopy_(nneq*nLoc,0.0_wp,0,WR,1)
          Call dcopy_(nTempMagn,0.0_wp,0,Zex,1)
          Call dcopy_(nneq*nTempMagn,0.0_wp,0,ZL,1)
          Call dcopy_(nneq*nTempMagn,0.0_wp,0,ZR,1)
          Call dcopy_(3*nTempMagn,0.0_wp,0,Sex,1)
          Call dcopy_(3*nTempMagn,0.0_wp,0,Mex,1)
          Call dcopy_(3*nneq*nTempMagn,0.0_wp,0,SL,1)
          Call dcopy_(3*nneq*nTempMagn,0.0_wp,0,ML,1)
          Call dcopy_(3*nneq*nTempMagn,0.0_wp,0,SR,1)
          Call dcopy_(3*nneq*nTempMagn,0.0_wp,0,MR,1)
          If(DBG) Write(6,'(A,F20.10)') 'zJ = ', zJ
c exchange magnetization:
          Call MAGN( exch, NM, dHX(iM),dHY(iM),dHZ(iM),H(iH),
     &               W, zJ, THRS,
     &               DIPEXCH,
     &                S_EXCH,
     &               nTempMagn,TempMagn,smagn,
     &               Wex,
     &               Zex,
     &               Sex,
     &               Mex, m_paranoid, DBG )
          If(DBG) Write(6,'(A,3F11.7)') 'MEX:',(Mex(l,1),l=1,3)
          Call Add_Info('MEX_MAGN    ',dnrm2_(3*nTempMagn,Mex,1),1,8)
          Call Add_Info('MR_MAGN  Wex',dnrm2_(nM,Wex,1)         ,1,8)
c compute local magnetizations:
          If(m_accurate) Then
            Do i=1,nneq
c all states:
              If ( NSS(i).gt.NEXCH(i) )  Then
c this check is to avoid the unnecessary computation, in cases when no local excited states are present
                Call MAGN( NSS(i), NEXCH(i), dHX(iM),dHY(iM),dHZ(iM),
     &                     H(iH),
     &                     ESO(i,1:NSS(i)), zJ, THRS,
     &                     DIPSO( i,1:3,1:NSS(i),1:NSS(i) ),
     &                      S_SO( i,1:3,1:NSS(i),1:NSS(i) ),
     &                     nTempMagn,TempMagn(1:nTempMagn),smagn,
     &                     WL(i,1:NEXCH(i) ),
     &                     ZL(i,1:nTempMagn),
     &                     SL(i,1:3,1:nTempMagn),
     &                     ML(i,1:3,1:nTempMagn), m_paranoid, DBG )
                Call Add_Info('MR_MAGN  WL',dnrm2_(nexch(i),WL,1),1,8)
                If(DBG) Write(6,'(A,I2,A,3F11.7)') 'ML: site',i,' : ',
     &                                   (ML(i,l,1),l=1,3)
c only local "exchange states":
                Call MAGN( NEXCH(i), NEXCH(i), dHX(iM),dHY(iM),dHZ(iM),
     &                     H(iH),
     &                     ESO(i,1:NEXCH(i)), zJ, THRS,
     &                     DIPSO( i,1:3,1:NEXCH(i),1:NEXCH(i) ),
     &                      S_SO( i,1:3,1:NEXCH(i),1:NEXCH(i) ),
     &                     nTempMagn,TempMagn,smagn,
     &                     WR(i,1:Nexch(i)),
     &                     ZR(i,1:nTempMagn),
     &                     SR(i,1:3,1:nTempMagn),
     &                     MR(i,1:3,1:nTempMagn), m_paranoid, DBG )
                Call Add_Info('MR_MAGN  WR',dnrm2_(nexch(i),WR,1),1,8)
                If(DBG) Write(6,'(A,I2,A,3F11.7)') 'MR: site',i,' : ',
     &                                   (MR(i,l,1),l=1,3)
              End If
            End Do

            Call Add_Info('ML_MAGN',dnrm2_(3*nTempMagn*nneq,ML,1),1,8)
            Call Add_Info('MR_MAGN',dnrm2_(3*nTempMagn*nneq,MR,1),1,8)

c expand the basis and rotate local vectors to the general
c coordinate system:
            Call dcopy_(  nCenter*nTempMagn,0.0_wp,0,ZRT,1)
            Call dcopy_(  nCenter*nTempMagn,0.0_wp,0,ZLT,1)
            Call dcopy_(3*nCenter*nTempMagn,0.0_wp,0,MRT,1)
            Call dcopy_(3*nCenter*nTempMagn,0.0_wp,0,MLT,1)
            Call dcopy_(3*nCenter*nTempMagn,0.0_wp,0,SRT,1)
            Call dcopy_(3*nCenter*nTempMagn,0.0_wp,0,SLT,1)
            isite=0
            Do i=1,NNEQ
              Do j=1,NEQ(i)
                isite=isite+1
c statistical distributions
                Do iT=1,nTempMagn
                  ZLT(isite,iT)=ZL(i,iT)
                  ZRT(isite,iT)=ZR(i,iT)
                End Do
c magnetizations:
c    use R_rot matrices, which have determinant +1.
c  note that  R_lg matrices have arbitrary determinant.
                Do iT=1,nTempMagn
                  Do l=1,3
                    Do n=1,3
                      MLT(isite,l,iT) = MLT(isite,l,iT) +
     &                                         r_rot(i,j,l,n)*ML(i,n,iT)
                      SLT(isite,l,iT) = SLT(isite,l,iT) +
     &                                         r_rot(i,j,l,n)*SL(i,n,iT)
                      MRT(isite,l,iT) = MRT(isite,l,iT) +
     &                                         r_rot(i,j,l,n)*MR(i,n,iT)
                      SRT(isite,l,iT) = SRT(isite,l,iT) +
     &                                         r_rot(i,j,l,n)*SR(i,n,iT)
                    End Do
                  End Do
                End Do
              End Do ! j, neq(i)
            End Do ! i, nneq
          End If ! m_accurate

c compute the total magnetizations according to the derived formulas:
          If (m_accurate) Then
            Do iT=1,nTempMagn
              If (smagn) Then
                Call MSUM ( nCenter, Sex(:,iT), Zex(iT),
     &                      SLT(:, :,iT), ZLT(:,iT),
     &                      SRT(:, :,iT), ZRT(:,iT), iopt,
     &                      ST(:,iH,iT),  ZT(iH,iT)  )
              End If
              Call MSUM ( nCenter, Mex(:,iT), Zex(iT),
     &                    MLT(:, :,iT), ZLT(:,iT),
     &                    MRT(:, :,iT), ZRT(:,iT), iopt,
     &                    MT(:,iH,iT),  ZT(iH,iT)  )
            End Do
          Else
c add the contribution from local excited states using the approximate
c X*H expression:
            Call dcopy_(3*nCenter*nTempMagn,0.0_wp,0,SRT,1)
            Call dcopy_(3*nCenter*nTempMagn,0.0_wp,0,SLT,1)
            Do iT=1,nTempMagn
              Do isite=1,nCenter
                Do l=1,3
                  MRT(isite,l,iT)=( XRM(isite,iT,l,1)*H(iH)*dHX(iM)
     &                             +XRM(isite,iT,l,2)*H(iH)*dHY(iM)
     &                             +XRM(isite,iT,l,3)*H(iH)*dHZ(iM)
     &                            ) /cm3tomB
                  MLT(isite,l,iT)=( XLM(isite,iT,l,1)*H(iH)*dHX(iM)
     &                             +XLM(isite,iT,l,2)*H(iH)*dHY(iM)
     &                             +XLM(isite,iT,l,3)*H(iH)*dHZ(iM)
     &                            ) /cm3tomB
                End Do
              End Do
              If (smagn) Then
                Call MSUM ( nCenter, Sex(:,iT), Zex(iT),
     &                      SLT(:, :,iT), ZLM(:,iT),
     &                      SRT(:, :,iT), ZRM(:,iT), iopt,
     &                       ST(:,iH,iT),  ZT(iH,iT)  )
              End If
              Call MSUM ( nCenter, Mex(:,iT), Zex(iT),
     &                    MLT(:, :,iT), ZLM(:,iT),
     &                    MRT(:, :,iT), ZRM(:,iT), iopt,
     &                     MT(:,iH,iT),  ZT(iH,iT)  )
            End Do
          End If
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c print out hte Zeeman eigenstates
      If (zeeman_energy) Then
         If((iH.eq.1).and.(iM.eq.nDir+1))
     &     Write(6,'(A)') 'Energies of the Zeeman Hamiltonian '//
     &                    'for the following directions of the '//
     &                    'applied field:'
         If ( (iH.eq.1).and.(iM.gt.nDir) .and.
     &                      (iM.le.(nDir+nDirZee)) ) Then
           Write(6,'(A,I3,A,3F10.6,3x,5A)')
     &             'direction Nr.',iM-nDir, ' : ',
     &              dHX(iM), dHY(iM), dHZ(iM),
     &             'written in file "zeeman_energy_',
     &              CHAR(48+mod(int((iM-nDir)/100),10)),
     &              CHAR(48+mod(int((iM-nDir)/10 ),10)),
     &              CHAR(48+mod(    (iM-nDir)     ,10)),
     &              '.txt".'
           Write(LUZee(iM-nDir),'(A,3F24.15)')
     &              '# direction of the applied magnetic field:',
     &              dHX(iM), dHY(iM), dHZ(iM)
           Write(LuZee(iM-nDir),'(A,6x,A,1000(I4,6x) )') '# H(T)',
     &              ' State =>',(i,i=1,nm)
         End If

         If ( (iM.gt.nDir) .and. (iM.le.(nDir+nDirZee)) ) Then
           Write(LUZee(iM-nDir),'(F8.4,1000F10.3)') H(IH),
     &                                     (Wex(I),I=1,NM)
         End If
      End If !zeeman_energy
c --------------------------------------------------------------------
c  computing the AVERAGE MOMENTS calculated at dIfferent temperatures
c (TempMagn(i))
      Do iT=1,nTempMagn
        Do l=1,3
          MVEC(iM,iH,iT,l)= MT(l,iH,iT)
          SVEC(iM,iH,iT,l)= ST(l,iH,iT)
        End Do !l

c      MAV(iH,iT)=MAV(iH,iT) + MT(1,iH,iT)*dHX(iM)*dHW(iM)
c     &                      + MT(2,iH,iT)*dHY(iM)*dHW(iM)
c     &                      + MT(3,iH,iT)*dHZ(iM)*dHW(iM)
c      SAV(iH,iT)=SAV(iH,iT) + ST(1,iH,iT)*dHX(iM)*dHW(iM)
c     &                      + ST(2,iH,iT)*dHY(iM)*dHW(iM)
c     &                      + ST(3,iH,iT)*dHZ(iM)*dHW(iM)
        ! accumulate contributions:
        Call daxpy_(1,dHX(iM)*dHW(iM),MT(1,iH,iT),1,MAV(iH,iT),1)
        Call daxpy_(1,dHY(iM)*dHW(iM),MT(2,iH,iT),1,MAV(iH,iT),1)
        Call daxpy_(1,dHZ(iM)*dHW(iM),MT(3,iH,iT),1,MAV(iH,iT),1)
        Call daxpy_(1,dHX(iM)*dHW(iM),ST(1,iH,iT),1,SAV(iH,iT),1)
        Call daxpy_(1,dHY(iM)*dHW(iM),ST(2,iH,iT),1,SAV(iH,iT),1)
        Call daxpy_(1,dHZ(iM)*dHW(iM),ST(3,iH,iT),1,SAV(iH,iT),1)
      End Do !iT
c ///  closing the loops over field strengths and directions
      End Do ! iM
      End Do ! iH
c Close Zeeman files, if opened
      If (Zeeman_Energy) Then
        Do i=1,nDirZee
          Close( LUZee(i) )
        End Do
      End If


c -------------------------------------------------------------------
C   WRITING SOME OF THE OUTPUT....
C -------------------------------------------------------------------
      If (smagn) Then
      Do iT=1,nTempMagn
c      Write(6,*)
        Do iDir=1,nDir+nDirZee
c        Write(6,*)
        Write(6,'(A,A,1x,A)') '--------|',
     &  '------------------------------------------------------------|',
     & '|------------------------------------------------------------|'
        Write(6,'(A,i3,26x,A,1x,A,60x,A)')
     &  'Direction of the applied magnetic field:',iDir,'|','|','|'
        Write(6,'(A,F18.14,44x,A,1x,A,60x,A)')
     &                              'proj X=',dHX(iDIR),'|','|','|'
        Write(6,'(A,F18.14,44x,A,1x,A,60x,A)')
     &                              'proj Y=',dHY(iDir),'|','|','|'
        Write(6,'(A,F18.14,44x,A,1x,A,60x,A)')
     &                              'proj Z=',dHZ(iDir),'|','|','|'
        Write(6,'(A,F7.4,A,41x,A,1x,A,60x,A)')
     &          'Temperature = ',TempMagn(iT),' Kelvin','|','|','|'
        Write(6,'(A,A,1x,A)') '--------|',
     &  '------------------------------------------------------------|',
     & '|------------------------------------------------------------|'
        Write(6,'(2x,A,12x,2A,1x,A,10x,2A)') 'Field |',
     & 'Magnetization Vector            |','   Total Magn. |','|',
     & 'Spin Magnetization Vector         |','   Total Magn. |'
        Write(6,'(5A,1x,5A)') '--------|',
     &  '--- proj X ---|','--- proj Y ---|','--- proj Z ---|',
     &  '- in this dir.-|',
     &  '|--- proj X ---|','--- proj Y ---|','--- proj Z ---|',
     &  '- in this dir.-|'
         Do iH=1,nH
         Write(6,'(F7.3,1x,A, 3(E13.6,1x,A),E14.7,1x,A,1x,A, '//
     &                       '3(E13.6,1x,A),E14.7,1x,A)')
     &   H(iH),'|',
     &     MVEC(iDir,iH,iT,1),' ',
     &     MVEC(iDir,iH,iT,2),' ',
     &     MVEC(iDir,iH,iT,3),'|',
     &   ( MVEC(iDir,iH,iT,1)*dHX(iDir)+
     &     MVEC(iDir,iH,iT,2)*dHY(iDir)+
     &     MVEC(iDir,iH,iT,3)*dHZ(iDir) ),'|','|',
     &     SVEC(iDir,iH,iT,1),' ',
     &     SVEC(iDir,iH,iT,2),' ',
     &     SVEC(iDir,iH,iT,3),'|',
     &   ( SVEC(iDir,iH,iT,1)*dHX(iDir)+
     &     SVEC(iDir,iH,iT,2)*dHY(iDir)+
     &     SVEC(iDir,iH,iT,3)*dHZ(iDir) ),'|'
           End Do
        Write(6,'(A,A,1x,A)') '--------|',
     &  '------------------------------------------------------------|',
     & '|------------------------------------------------------------|'
          End Do !iDir
        End Do !iT

        Else

      Do iT=1,nTempMagn
        Do iDir=1,nDir+nDirZee
        Write(6,'(2A)') '--------|',
     & '------------------------------------------------------------|'
        Write(6,'(A,i3,26x,A)')
     &  'Direction of the applied magnetic field:',iDir,'|'
        Write(6,'(A,F18.14,44x,A)') 'proj X=',dHX(iDIR),'|'
        Write(6,'(A,F18.14,44x,A)') 'proj Y=',dHY(iDir),'|'
        Write(6,'(A,F18.14,44x,A)') 'proj Z=',dHZ(iDir),'|'
        Write(6,'(A,F7.4,A,41x,A)') 'Temperature = ',TempMagn(iT),
     &                       ' Kelvin','|'
       Write(6,'(2A)') '--------|',
     & '------------------------------------------------------------|'
        Write(6,'(2x,A,12x,2A)') 'Field |',
     & 'Magnetization Vector            |','   Total Magn. |'
        Write(6,'(5A)') '--------|',
     &  '--- proj X ---|','--- proj Y ---|','--- proj Z ---|',
     &  '- in this dir.-|'
         Do iH=1,nH
            Write(6,'(F7.3,1x,A,3(E13.6,1x,A),E14.7,1x,A)')
     &                         H(iH),'|',
     &              MVEC(iDir,iH,iT,1),' ',
     &              MVEC(iDir,iH,iT,2),' ',
     &              MVEC(iDir,iH,iT,3),'|',
     &            ( MVEC(iDir,iH,iT,1)*dHX(iDir)+
     &              MVEC(iDir,iH,iT,2)*dHY(iDir)+
     &              MVEC(iDir,iH,iT,3)*dHZ(iDir) ),'|'
         End Do
         Write(6,'(2A)') '--------|',
     & '------------------------------------------------------------|'
          End Do !iDir

        End Do !iT
       End If !(smagn)
c ---------------------------------------------------------------------
c     COMPUTING THE STANDARD DEVIATION OF THE MAGNETIZATION
c      If(HINPUT) Then
c        Do iT=1,nTempMagn
c          STDEV(iT) =0.0_wp
c          STDEV(iT) = dev(nH,MAV(:,iT),Mexp(:,iT))
c        End Do
c      End If

      Write(6,*)
      Write(6,'(15X,A)') 'HIGH-FIELD POWDER MAGNETIZATION'
      Write(6,'(20X,A)') '(Units: Bohr magneton)'
      Write(6,*)
      Do iT=1,nTempMagn,5
         iTEnd=min(nTempMagn,iT+4)
         Write(6,'(A,11A)') '-----------|',('---------------|',i=iT,
     &                      iTEnd+1)
         Write(6,'(A,11(F10.3,A))') '    H(T)   |STATISTICAL SUM|',
     &                    (TempMagn(i),' K.  |',i=iT,iTEnd)
         Write(6,'(A,11A)') '-----------|',('---------------|',i=iT,
     &                      iTEnd+1)
         Do iH=1,nH
          Write(6,'(F10.6,1X,A,F14.7,1x,A,11(f14.10,1x,A))')
     &             H(iH),'|',ZT(iH,1),'|', (MAV(iH,i),'|',i=iT,iTEnd)
         End Do
         Write(6,'(A,11A)') '-----------|',('---------------|',i=iT,
     &                      iTEnd+1)
         If(HINPUT) Then
            Write(6,'(A,15x,A, 11(f14.10,1x,A) )')
     &                      'ST.DEV.M   |','|',
     &               ( dev(nH,MAV(:,i),Mexp(:,i)),'|',i=iT,iTEnd)
            Write(6,'(A,11A)') '-----------|',
     &                        ('---------------|',i=iT,iTEnd+1)
         End If
      End Do
      If(smagn) Then
         Write(6,*)
         Write(6,'(15X,A)') 'HIGH-FIELD POWDER SPIN MAGNETIZATION'
         Write(6,'(20X,A)') '(Units: Bohr magneton)'
         Write(6,*)
         Do iT=1,nTempMagn,5
           iTEnd=min(nTempMagn,iT+4)
           Write(6,'(A,11A)') '-----------|',('---------------|',i=iT,
     &                                                        iTEnd+1)
           Write(6,'(A,11(F10.3,A))') '   H(T)   |STATISTICAL SUM|',
     &                               (TempMagn(i),' K.  |',i=iT,iTEnd)
           Write(6,'(A,11A)') '-----------|',('---------------|',i=iT,
     &                                                        iTEnd+1)
           Do iH=1,nH
            Write(6,'(F10.6,1X,A,F14.7,1x,A,11(f14.10,1x,A))')
     &              H(iH),'|',ZT(iH,1),'|', (SAV(iH,i),'|',i=iT,iTEnd)
           End Do
           Write(6,'(A,11A)') '-----------|',('---------------|',i=iT,
     &                                                        iTEnd+1)
         End Do
      End If!smagn
      ! add some verification:
      Call Add_Info('H_MAGN       ',   H,nH,6)
      Call Add_Info('MAGN_AVERAGED', MAV,nH*nTempMagn,6)
      Call Add_Info('ZTL_MAGN',dnrm2_(nH*nTempMagn,ZT,1),1,8)
      Do iH=1,nH
         Write(lbl_X,'(A,i3)') 'MAGN_VECT X ',iH
         Write(lbl_Y,'(A,i3)') 'MAGN_VECT Y ',iH
         Write(lbl_Z,'(A,i3)') 'MAGN_VECT Z ',iH
         ibuf=nDirTot*nTempMagn
         Call Add_Info(lbl_X,dnrm2_(ibuf,MVEC(:,iH,:,1),1),1,8)
         Call Add_Info(lbl_Y,dnrm2_(ibuf,MVEC(:,iH,:,2),1),1,8)
         Call Add_Info(lbl_Z,dnrm2_(ibuf,MVEC(:,iH,:,3),1),1,8)
      End Do






c----------------------------------------------------------------------
      If(nM>0) Then
        Call mma_deallocate(Wex)
      End If

      If(nTempMagn>0) Then
        Call mma_deallocate(Zex)
        Call mma_deallocate(Sex)
        Call mma_deallocate(Mex)
        If(nneq>0) Then
          Call mma_deallocate(ZL)
          Call mma_deallocate(SL)
          Call mma_deallocate(ML)
          Call mma_deallocate(ZR)
          Call mma_deallocate(SR)
          Call mma_deallocate(MR)
        End If

        If(nLoc>0) Then
          Call mma_deallocate(WL)
          Call mma_deallocate(WR)
        End If

        If(nCenter>0) Then
          Call mma_deallocate(ZRT)
          Call mma_deallocate(ZLT)
          Call mma_deallocate(MRT)
          Call mma_deallocate(MLT)
          Call mma_deallocate(SRT)
          Call mma_deallocate(SLT)
        End If

        If(nH>0) Then
          Call mma_deallocate(ZT)
          Call mma_deallocate(ST)
          Call mma_deallocate(MT)
          Call mma_deallocate(SAV)
          Call mma_deallocate(MAV)
          If(nDirTot>0) Then
            Call mma_deallocate(SVEC)
            Call mma_deallocate(MVEC)
          End If
        End If
      End If

      If(nDirTot>0) Then
         Call mma_deallocate(dHX)
         Call mma_deallocate(dHY)
         Call mma_deallocate(dHZ)
         Call mma_deallocate(dHW)
      End If

      If(nH>0) Then
         Call mma_deallocate(H)
      End If

      Call qExit('magnetisation')
      Return
      End

