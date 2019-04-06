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
      Subroutine torque_pa( nneq, nCenter, neq, neqv, nLoc, exch,
     &                   nTempMagn, nH, nM, AngPoints, nexch,
     &                   iopt, nss, mem,
     &                   smagn, m_paranoid, m_accurate,
     &                   TempMagn, w, hmin, hmax, dltH0, EM, zJ, THRS,
     &                   hexp,
     &                   dipexch, s_exch, dipso, s_so, eso,
     &                   hinput, r_rot, XLM, ZLM, XRM, ZRM )

      Implicit None
      Integer, parameter           :: wp=SELECTED_REAL_KIND(p=15,r=307)
#include "mgrid.fh"
#include "stdalloc.fh"
c------------------------------------------------------------
c         INPUT VARIABLES:
c------------------------------------------------------------
      Integer, intent (in)         :: nneq, neq(nneq), neqv, nLoc,
     &                                nexch(nneq), nH, nTempMagn, nM,
     &                                exch, nCenter, AngPoints,
     &                                nss(nneq), iopt, mem
      Real(kind=wp), intent(in)    :: TempMagn(nTempMagn), hmin, hmax,
     &                                dltH0, EM, zJ, THRS, Hexp(nH)
      Logical, intent(in)          :: m_paranoid, m_accurate, smagn,
     &                                hinput
c correction to M from the local excited states:
      Real(kind=wp), intent(in)    :: XLM( nCenter,nTempMagn,3,3)
      Real(kind=wp), intent(in)    :: ZLM( nCenter,nTempMagn)
      Real(kind=wp), intent(in)    :: XRM( nCenter,nTempMagn,3,3)
      Real(kind=wp), intent(in)    :: ZRM( nCenter,nTempMagn)
c rotation matrices for equivalent sites:
      Real(kind=wp), intent(in)    :: R_ROT(nneq,neqv,3,3)
c exchange spectum:
!     exchange energies printed out in the previous part
      Real(kind=wp), intent(in)    :: W(exch)
c local spin-orbit spectum:
!     spin-orbit energies from ANISO files
      Real(kind=wp), intent(in)    :: ESO(nneq,nLoc)
c magnetic and spin moments (i.e. the BIG matrices):
      Complex(kind=wp), intent(in) :: DIPEXCH(3,EXCH,EXCH)
      Complex(kind=wp), intent(in) ::  S_EXCH(3,EXCH,EXCH)
      Complex(kind=wp), intent(in) :: dipso(nneq,3,nLoc,nLoc)
      Complex(kind=wp), intent(in) ::  s_so(nneq,3,nLoc,nLoc)

c------------------------------------------------------------
c         LOCAL VARIABLES:
c------------------------------------------------------------
c exchange data:
      Real(kind=wp), allocatable :: WEX(:)
!                                   WEX(NM)                ! Zeeman exchange energies
      Real(kind=wp), allocatable :: ZEX(:)
!                                   ZEX(nTempMagn)         ! exchange statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: SEX(:,:)
!                                   SEX(3,nTempMagn)       ! spin magnetisation, from the exchange block;
      Real(kind=wp), allocatable :: MEX(:,:)
!                                   MEX(3,nTempMagn)       ! magnetisation, form the exchange block
c data for individual sites (all states):
      Real(kind=wp), allocatable :: ZL(:,:)
!                                   ZL(nneq,nTempMagn)     ! local statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: WL(:,:)
!                                   WL(nneq,nLoc)          ! Zeeman local energies
      Real(kind=wp), allocatable :: SL(:,:,:)
!                                   SL(nneq,3,nTempMagn)   ! spin magnetisation, from the local sites, using ALL states ;
      Real(kind=wp), allocatable :: ML(:,:,:)
!                                   ML(nneq,3,nTempMagn)   ! magnetisation, from local sites, using ALL states;
c data for individual sites (only states that enter exchange):
      Real(kind=wp), allocatable :: ZR(:,:)
!                                   ZR(nneq,nTempMagn)     ! local statistical sum, Boltzmann distribution, using only NEXCH states
      Real(kind=wp), allocatable :: WR(:,:)
!                                   WR(nneq,nLoc)          ! Zeeman local reduced energies, using only NEXCH states;
      Real(kind=wp), allocatable :: SR(:,:,:)
!                                   SR(nneq,3,nTempMagn)   ! spin magnetisation, from the local sites, using only NEXCH states ;
      Real(kind=wp), allocatable :: MR(:,:,:)
!                                   MR(nneq,3,nTempMagn)   ! magnetisation, from local sites, using only NEXCH states;
c total vectors in general coordinate system:
      Real(kind=wp), allocatable :: ZRT(:,:)   ! ZRT(nCenter,nTempMagn)
      Real(kind=wp), allocatable :: ZLT(:,:)   ! ZLT(nCenter,nTempMagn)
      Real(kind=wp), allocatable :: MRT(:,:,:)
!                                   MRT(nCenter,3,nTempMagn)
      Real(kind=wp), allocatable :: MLT(:,:,:)
!                                   MLT(nCenter,3,nTempMagn)
      Real(kind=wp), allocatable :: SRT(:,:,:)
!                                   SRT(nCenter,3,nTempMagn)
      Real(kind=wp), allocatable :: SLT(:,:,:)
!                                   SLT(nCenter,3,nTempMagn)
c data for total system:
      Real(kind=wp), allocatable :: ZT(:)
!                                   ZT(nTempMagn)        ! total statistical sum, Boltzmann distribution
      Real(kind=wp), allocatable :: ST(:,:)
!                                   ST(3,nTempMagn)      ! total spin magnetisation,
      Real(kind=wp), allocatable :: MT(:,:)
!                                   MT(3,nTempMagn)      ! total magnetisation
c magnetic field strength and orientation data:
      Integer       :: nPlanes
      Real(kind=wp) :: dlth
      parameter (nPlanes=3)
      Real(kind=wp), allocatable :: H(:)    ! H(nH)
      Real(kind=wp), allocatable :: dX(:,:) ! dX(nPlanes,AngPoints)
      Real(kind=wp), allocatable :: dY(:,:) ! dY(nPlanes,AngPoints)
      Real(kind=wp), allocatable :: dZ(:,:) ! dZ(nPlanes,AngPoints)
      Real(kind=wp), allocatable :: Ang(:)  ! Ang(AngPoints)
c magnetic torque
      Real(kind=wp), allocatable :: tx(:,:,:,:)
!                                   tx(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque, X
      Real(kind=wp), allocatable :: ty(:,:,:,:)
!                                   ty(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque, Y
      Real(kind=wp), allocatable :: tz(:,:,:,:)
!                                   tz(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque, Z
      Real(kind=wp), allocatable :: sx(:,:,:,:)
!                                   sx(nPlanes,AngPoints,nH,nTempMagn) ! spin magnetization torque, X
      Real(kind=wp), allocatable :: sy(:,:,:,:)
!                                   sy(nPlanes,AngPoints,nH,nTempMagn) ! spin magnetization torque, Y
      Real(kind=wp), allocatable :: sz(:,:,:,:)
!                                   sz(nPlanes,AngPoints,nH,nTempMagn) ! spin magnetization torque, Z
      Real(kind=wp) :: cm3tomB
      Integer :: mem_local, RtoB
c local data:
      Integer :: IM,I,it
      Integer :: J,IH,k,isite,l,n,iPl
      Logical :: DBG


      Call qEnter('PA_torq')

c      Boltz_k=0.6950356000_wp   !   in cm^-1*K-1
c      mu_Bohr=0.4668643740_wp   !   in cm-1*T-1
      cm3tomB=0.5584938904_wp   !   in cm3 * mol-1 * T
      DBG=.false.

      Write(6,*)
      Write(6,'(100A)') (('%'),J=1,96)
      Write(6,'(20X,A)')
     & 'ANGULAR DEPENDENCE OF THE MAGNETIZATION TORQUE'
      Write(6,'(100A)') (('%'),J=1,96)
      Write(6,*)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Write(6,'(2X,A,i3,A)') 'Magnetization torque is '//
     & 'calculated for the ',NH,' field points, in the field Domain:'
      Write(6,'(2X,F4.1,1x,a,1X,F4.1,a,30(F6.3,a))') HMIN,'--',HMAX,
     & ' T., at the following temperatures:'
      Do i=11,nTempMagn,10
        j=MIN(nTempMagn,i+9)
        Write(6,'(17x,10(F9.3,A))') (TempMagn(k),' K.;',k=i,j)
      End Do
      Write(6,'(2x,A,i3,A)')
     & 'Angular dependence of the magnetization torque '//
     & 'is computed for ',AngPoints,' angular points distributed'
      Write(6,'(2x,A)') 'in the domain 0-180 deg.'
      Write(6,'(2X,A)') 'The cut-off energy for the '//
     & 'exact diagonalization of the Zeeman Hamiltonian is:'
      Write(6,'(2x,a,F15.9,A)') 'E = ',EM ,' cm(-1).'
      If(NM.lt.10) Then
        Write(6,'(2X,A,i2,a)') 'The exact diagonalization of the '//
     &       'Zeeman Hamiltonian included ',NM,' exchange states.'
      Else If( (NM.ge.10) .AND. (NM.lt.100) ) Then
        Write(6,'(2X,A,i3,a)') 'The exact diagonalization of the '//
     &       'Zeeman Hamiltonian included ',NM,' exchange states.'
      Else If( (NM.ge.100) .AND. (NM.lt.1000) ) Then
        Write(6,'(2X,A,i4,a)') 'The exact diagonalization of the '//
     &       'Zeeman Hamiltonian included ',NM,' exchange states.'
      Else If( (NM.ge.1000) .AND. (NM.lt.10000) ) Then
        Write(6,'(2X,A,i5,a)') 'The exact diagonalization of the '//
     &       'Zeeman Hamiltonian included ',NM,' exchange states.'
      End If
      If(m_accurate) Then
        Write(6,'(2x,A)') 'The contribution of local excited states'//
     &                    ' is computed exactly.'
      Else
        Write(6,'(2x,A)') 'The contribution of local excited states'//
     &                    ' is computed approximately (by using the'//
     &                    ' susceptibility data). '
      End If
!-----------------------------------------------------------------------
      If(dbg) Write(6,*) 'nM       = ', nM
      If(dbg) Write(6,*) 'nTempMagn= ', nTempMagn
      If(dbg) Write(6,*) 'nneq     = ', nneq
      If(dbg) Write(6,*) 'nCenter  = ', nCenter
      If(dbg) Write(6,*) 'AngPoints= ', AngPoints
      If(dbg) Write(6,*) 'nPlanes  = ', nPlanes
      If(dbg) Write(6,*) 'nH       = ', nH
      If(dbg) Write(6,*) 'nLoc     = ', nLoc
      If(dbg) Write(6,*) 'nexch()  = ', (nexch(i),i=1,nneq)
      If(dbg) Write(6,*) 'neq()    = ', (neq(i),i=1,nneq)
      If(dbg) Write(6,*) 'exch     = ', exch
      If(dbg) Write(6,*) 'iopt     = ', iopt
      If(dbg) Write(6,*) 'neqv     = ', neqv
      If(dbg) Write(6,*) 'nss()    = ', (nss(i),i=1,nneq)
      If(dbg) Write(6,*) 'W()      = ', (W(i),i=1,exch)
      If(dbg) Write(6,*) 'zJ       = ', zJ
      If(dbg) Write(6,*) 'EM       = ', EM
      If(dbg) Write(6,*) 'm_paranoi= ', m_paranoid
      If(dbg) Write(6,*) 'm_accurat= ', m_accurate
      If(dbg) Write(6,*) 'smagn    = ', smagn
      If(dbg) Write(6,*) 'hinput   = ', hinput


! Allocate memory for this calculation:
      mem_local=0
      RtoB=8
      If(nM>0) Then
         ! Zeeman exchange energy spectrum
         Call mma_allocate(Wex,nM,'Wex')
         Call dcopy_(nM,[0.0_wp],0,Wex,1)
         mem_local=mem_local+nM*RtoB
      End If

      If(dbg) Write(6,*) 'mem_local 1 = ', mem_local

      If(nTempMagn>0) Then
         ! exchange statistical sum, Boltzmann distribution
         Call mma_allocate(Zex,nTempMagn,'Zex')
         Call dcopy_(nTempMagn,[0.0_wp],0,Zex,1)
         mem_local=mem_local+nTempMagn*RtoB
         ! spin magnetisation, from the exchange block
         Call mma_allocate(SEX,3,nTempMagn,'SEX')
         Call dcopy_(3*nTempMagn,[0.0_wp],0,SEX,1)
         mem_local=mem_local+3*nTempMagn*RtoB
         ! magnetisation, from the exchange block
         Call mma_allocate(MEX,3,nTempMagn,'MEX')
         Call dcopy_(3*nTempMagn,[0.0_wp],0,MEX,1)
         mem_local=mem_local+3*nTempMagn*RtoB

         ! total statistical sum, Boltzmann distribution
         Call mma_allocate(ZT,nTempMagn,'ZT')
         Call dcopy_(nTempMagn,[0.0_wp],0,ZT,1)
         mem_local=mem_local+nTempMagn*RtoB
         ! total spin magnetisation
         Call mma_allocate(ST,3,nTempMagn,'ST')
         Call dcopy_(3*nTempMagn,[0.0_wp],0,ST,1)
         mem_local=mem_local+3*nTempMagn*RtoB
         ! total magnetisation
         Call mma_allocate(MT,3,nTempMagn,'MT')
         Call dcopy_(3*nTempMagn,[0.0_wp],0,MT,1)
         mem_local=mem_local+3*nTempMagn*RtoB

      if(dbg) Write(6,*) 'mem_local 2 = ', mem_local
         If(nneq>0) Then
            ! local statistical sum, Boltzmann distribution
            Call mma_allocate(ZL,nneq,nTempMagn,'ZL')
            Call dcopy_(nneq*nTempMagn,[0.0_wp],0,ZL,1)
            mem_local=mem_local+nneq*nTempMagn*RtoB
            ! spin magnetisation, from the local sites, using ALL states
            Call mma_allocate(SL,nneq,3,nTempMagn,'SL')
            Call dcopy_(nneq*3*nTempMagn,[0.0_wp],0,SL,1)
            mem_local=mem_local+nneq*3*nTempMagn*RtoB
            ! magnetisation, from local sites, using ALL states
            Call mma_allocate(ML,nneq,3,nTempMagn,'ML')
            Call dcopy_(nneq*3*nTempMagn,[0.0_wp],0,ML,1)
            mem_local=mem_local+nneq*3*nTempMagn*RtoB

            ! local statistical sum, Boltzmann distribution
            Call mma_allocate(ZR,nneq,nTempMagn,'ZR')
            Call dcopy_(nneq*nTempMagn,[0.0_wp],0,ZR,1)
            mem_local=mem_local+nneq*nTempMagn*RtoB
!           spin magnetisation, from the local sites, using only NEXCH states
            Call mma_allocate(SR,nneq,3,nTempMagn,'SR')
            Call dcopy_(nneq*3*nTempMagn,[0.0_wp],0,SR,1)
            mem_local=mem_local+nneq*3*nTempMagn*RtoB
!           magnetisation, from the local sites, using only NEXCH states
            Call mma_allocate(MR,nneq,3,nTempMagn,'MR')
            Call dcopy_(nneq*3*nTempMagn,[0.0_wp],0,MR,1)
            mem_local=mem_local+nneq*3*nTempMagn*RtoB
         End If

      if(dbg) Write(6,*) 'mem_local 3 = ', mem_local
         If(nCenter>0) Then
            ! ZRT(nCenter,nTempMagn)
            Call mma_allocate(ZRT,nCenter,nTempMagn,'ZRT')
            Call dcopy_(nCenter*nTempMagn,[0.0_wp],0,ZRT,1)
            mem_local=mem_local+nCenter*nTempMagn*RtoB
            ! ZLT(nCenter,nTempMagn)
            Call mma_allocate(ZLT,nCenter,nTempMagn,'ZLT')
            Call dcopy_(nCenter*nTempMagn,[0.0_wp],0,ZLT,1)
            mem_local=mem_local+nCenter*nTempMagn*RtoB
            ! MRT(nCenter,3,nTempMagn)
            Call mma_allocate(MRT,nCenter,3,nTempMagn,'MRT')
            Call dcopy_(nCenter*3*nTempMagn,[0.0_wp],0,MRT,1)
            mem_local=mem_local+nCenter*3*nTempMagn*RtoB
            ! MLT(nCenter,3,nTempMagn)
            Call mma_allocate(MLT,nCenter,3,nTempMagn,'MLT')
            Call dcopy_(nCenter*3*nTempMagn,[0.0_wp],0,MLT,1)
            mem_local=mem_local+nCenter*3*nTempMagn*RtoB
            ! SRT(nCenter,3,nTempMagn)
            Call mma_allocate(SRT,nCenter,3,nTempMagn,'SRT')
            Call dcopy_(nCenter*3*nTempMagn,[0.0_wp],0,SRT,1)
            mem_local=mem_local+nCenter*3*nTempMagn*RtoB
            ! SLT(nCenter,3,nTempMagn)
            Call mma_allocate(SLT,nCenter,3,nTempMagn,'SLT')
            Call dcopy_(nCenter*3*nTempMagn,[0.0_wp],0,SLT,1)
            mem_local=mem_local+nCenter*3*nTempMagn*RtoB
         End If

         If(dbg) Write(6,*) 'mem_local 4 = ', mem_local
c magnetic torque
         If( (nPlanes>0).and.(AngPoints>0).and.(nH>0)
     &                                    .and.(nTempMagn>0) ) Then
            Call mma_allocate(tx,nPlanes,AngPoints,nH,nTempMagn,'tx')
            Call mma_allocate(ty,nPlanes,AngPoints,nH,nTempMagn,'ty')
            Call mma_allocate(tz,nPlanes,AngPoints,nH,nTempMagn,'tz')
            Call mma_allocate(sx,nPlanes,AngPoints,nH,nTempMagn,'sx')
            Call mma_allocate(sy,nPlanes,AngPoints,nH,nTempMagn,'sy')
            Call mma_allocate(sz,nPlanes,AngPoints,nH,nTempMagn,'sz')
            Call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[0.0_wp],0,tx,1)
            Call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[0.0_wp],0,ty,1)
            Call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[0.0_wp],0,tz,1)
            Call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[0.0_wp],0,sx,1)
            Call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[0.0_wp],0,sy,1)
            Call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[0.0_wp],0,sz,1)
            mem_local=mem_local+6*nPlanes*AngPoints*nH*nTempMagn*RtoB
         End If
      End If
      if(dbg) Write(6,*) 'mem_local 5 = ', mem_local

      If((nLoc>0).and.(nneq>0)) Then
         ! Zeeman local energies
         Call mma_allocate(WL,nneq,nLoc,'WL')
         Call dcopy_(nneq*nLoc,[0.0_wp],0,WL,1)
         mem_local=mem_local+nneq*nLoc*RtoB
         ! Zeeman local reduced energies, using only NEXCH states
         Call mma_allocate(WR,nneq,nLoc,'WR')
         Call dcopy_(nneq*nLoc,[0.0_wp],0,WR,1)
         mem_local=mem_local+nneq*nLoc*RtoB
      End If
      if(dbg) Write(6,*) 'mem_local 6 = ', mem_local

      If(AngPoints>0) Then
         Call mma_allocate(Ang,AngPoints,'Ang')
         Call dcopy_(AngPoints,[0.0_wp],0,Ang,1)
         mem_local=mem_local+AngPoints*RtoB
         If(nPlanes>0) Then
            Call mma_allocate(dX,nPlanes,AngPoints,'dX')
            Call mma_allocate(dY,nPlanes,AngPoints,'dY')
            Call mma_allocate(dZ,nPlanes,AngPoints,'dZ')
            Call dcopy_(nPlanes*AngPoints,[0.0_wp],0,dX,1)
            Call dcopy_(nPlanes*AngPoints,[0.0_wp],0,dY,1)
            Call dcopy_(nPlanes*AngPoints,[0.0_wp],0,dZ,1)
            mem_local=mem_local+3*nPlanes*AngPoints*RtoB
         End If
      End If
      if(dbg) Write(6,*) 'mem_local 7 = ', mem_local

      If(nH>0) Then
         Call mma_allocate(H,nH,'H')
         Call dcopy_(nH,[0.0_wp],0,H,1)
         mem_local=mem_local+nH*RtoB
      End If

      If(dbg) Write(6,*) 'TORQ:  memory allocated (local):'
      If(dbg) Write(6,*) 'mem_local=', mem_local
      If(dbg) Write(6,*) 'TORQ:  memory allocated (total):'
      If(dbg) Write(6,*) 'mem_total=', mem+mem_local


!-----------------------------------------------------------------------
! set up the field points:
      If (HINPUT) Then
        Do iH=1,nH
          H(iH)=HEXP(iH)
          If(H(iH).eq.0.0_wp) Then
            H(iH)=0.0001_wp
          End If
        End Do
      Else
        DLTH=(HMAX-HMIN)/DBLE(NH-1)
        Do iH=1,nH
          If (iH.eq.1) Then
            H(IH)=HMIN+dltH0
          Else
            H(IH)=HMIN+DLTH*DBLE(IH-1)
          End If
          If(H(iH).eq.0.0_wp) Then
            H(iH)=0.0001_wp
          End If
        End Do
      End If
c-----------------------------------------

      Do IH=1,NH
c ///  opening the loop over field points:
         Call dcopy_(nPlanes*AngPoints,[0.0_wp],0,dX,1)
         Call dcopy_(nPlanes*AngPoints,[0.0_wp],0,dY,1)
         Call dcopy_(nPlanes*AngPoints,[0.0_wp],0,dZ,1)
         Do iPl=1,nPlanes
c ///  opening the loop over dIfferent planes of rotation of the applied magnetic field:
c  loop over various angular grids  (iPl = 1, 2 or 3)
c  iPl=1 (X) , angular grid in the plane YZ ( half of the plane , 0-180 degrees)
c  iPl=2 (Y) , angular grid in the plane XZ ( half of the plane , 0-180 degrees)
c  iPl=3 (Z) , angular grid in the plane XY ( half of the plane , 0-180 degrees)
            Call dcopy_(AngPoints,[0.0_wp],0,Ang,1)

            Call hdir2( AngPoints, iPl, dX(iPl,:), dY(iPl,:), dZ(iPl,:),
     &                  Ang, 4 )

            If(dbg) Write(6,*) 'iPl, dX, dY, dZ=', iPl, dX(iPl,:),
     &                          dY(iPl,:), dZ(iPl,:),  Ang
            Do IM=1,AngPoints
c ///  opening the loop over dIfferent directions of the magnetic field
               Call dcopy_(nM,              [0.0_wp],0,Wex,1)
               Call dcopy_(nneq*nLoc,       [0.0_wp],0,WL ,1)
               Call dcopy_(nneq*nLoc,       [0.0_wp],0,WR ,1)
               Call dcopy_(       nTempMagn,[0.0_wp],0,Zex,1)
               Call dcopy_(nneq*  nTempMagn,[0.0_wp],0,ZL ,1)
               Call dcopy_(nneq*  nTempMagn,[0.0_wp],0,ZR ,1)
               Call dcopy_(     3*nTempMagn,[0.0_wp],0,SEX,1)
               Call dcopy_(     3*nTempMagn,[0.0_wp],0,MEX,1)
               Call dcopy_(nneq*3*nTempMagn,[0.0_wp],0,SL ,1)
               Call dcopy_(nneq*3*nTempMagn,[0.0_wp],0,ML ,1)
               Call dcopy_(nneq*3*nTempMagn,[0.0_wp],0,SR ,1)
               Call dcopy_(nneq*3*nTempMagn,[0.0_wp],0,MR ,1)
               Call dcopy_(     3*nTempMagn,[0.0_wp],0,MT ,1)
               Call dcopy_(     3*nTempMagn,[0.0_wp],0,ST ,1)
               Call dcopy_(       nTempMagn,[0.0_wp],0,ZT ,1)

c exchange magnetization:
               Call MAGN( EXCH, NM,
     &                    dX(iPl,iM),dY(iPl,iM),dZ(iPl,iM), H(iH),
     &                    W, zJ, THRS,
     &                    DIPEXCH,
     &                     S_EXCH,
     &                    nTempMagn,TempMagn,smagn,
     &                    Wex,
     &                    Zex,
     &                    Sex,
     &                    Mex, m_paranoid , DBG )

          If(iPl==2) Then
          If(DBG) Write(6,'(A,I3,1x,F8.4,2x, 3F19.14,2x,3F19.14)')
     &            'MEX: iM,',iM, H(iH), (Mex(l,1),l=1,3),
     &                   dX(iPl,iM),dY(iPl,iM),dZ(iPl,iM)
          End If
c              iT=1

c              Write(6,'(F10.4,1x,2I3,3F18.14,3x,3F20.14,2x,F20.14)')
c     &           H(iH), iL, iM, dX(iM),dY(iM),dZ(iM),(Mex(j,iT),j=1,3),
c     &           Mex(1,iT)*dX(iM)+Mex(2,iT)*dY(iM)+Mex(3,iT)*dZ(iM)
               ! compute local magnetizations:
               If (m_accurate) Then
                  Do i=1,nneq
                     ! all states:
                     If ( NSS(i) > NEXCH(i) )  Then
                     ! this check is meant to avoid the unnecessary
                     ! computation, in cases when no local excited
                     ! states are present
                       Call MAGN( NSS(i), NEXCH(i),
     &                            dX(iPl,iM),dY(iPl,iM),dZ(iPl,iM),
     &                            H(iH), ESO(i,1:NSS(i)), zJ, THRS,
     &                            DIPSO(i,:,:,:),
     &                             S_SO(i,:,:,:),
     &                            nTempMagn,TempMagn,smagn,
     &                            WL(i,:),
     &                            ZL(i,:),
     &                            SL(i,:,:),
     &                            ML(i,:,:), m_paranoid, DBG )
                     ! only local "exchange states":
                       Call MAGN( NEXCH(i), NEXCH(i),
     &                            dX(iPl,iM),dY(iPl,iM),dZ(iPl,iM),
     &                            H(iH), ESO(i,1:NEXCH(i)), zJ, THRS,
     &                            DIPSO( i,1:3,1:NEXCH(i),1:NEXCH(i) ),
     &                             S_SO( i,1:3,1:NEXCH(i),1:NEXCH(i) ),
     &                            nTempMagn,TempMagn,smagn,
     &                            WR(i,1:Nexch(i)),
     &                            ZR(i,:),
     &                            SR(i,:,:),
     &                            MR(i,:,:), m_paranoid, DBG )
                     End If
                  End Do
                  ! expand the basis and rotate local vectors to the
                  ! general coordinate system:
                  Call dcopy_(  nCenter*nTempMagn,[0.0_wp],0,ZRT,1)
                  Call dcopy_(  nCenter*nTempMagn,[0.0_wp],0,ZLT,1)
                  Call dcopy_(3*nCenter*nTempMagn,[0.0_wp],0,MRT,1)
                  Call dcopy_(3*nCenter*nTempMagn,[0.0_wp],0,MLT,1)
                  Call dcopy_(3*nCenter*nTempMagn,[0.0_wp],0,SRT,1)
                  Call dcopy_(3*nCenter*nTempMagn,[0.0_wp],0,SLT,1)

                  isite=0
                  Do i=1,NNEQ
                     Do j=1,NEQ(i)
                        isite=isite+1
                        ! statistical distributions
                        Do iT=1,nTempMagn
                           ZLT(isite,iT)=ZL(i,iT)
                           ZRT(isite,iT)=ZR(i,iT)
                        End Do
!                       magnetizations:
!                       use R_rot matrices, which have determinant +1.
!                        >> note that  R_lg matrices may have arbitrary
!                        >> sign of the determinant.
                        Do iT=1,nTempMagn
                           Do l=1,3
                              Do n=1,3
                                 MLT(isite,l,iT) = MLT(isite,l,iT)
     &                                       + r_rot(i,j,l,n)*ML(i,n,iT)

                                 SLT(isite,l,iT) = SLT(isite,l,iT)
     &                                       + r_rot(i,j,l,n)*SL(i,n,iT)

                                 MRT(isite,l,iT) = MRT(isite,l,iT)
     &                                       + r_rot(i,j,l,n)*MR(i,n,iT)

                                 SRT(isite,l,iT) = SRT(isite,l,iT)
     &                                       + r_rot(i,j,l,n)*SR(i,n,iT)
                              End Do
                           End Do
                        End Do

                     End Do ! j, neq(i)
                  End Do ! i, nneq
               End If ! m_accurate

               ! compute the total magnetizations according
               ! to the derived formulas:
               If (m_accurate) Then
                  Do iT=1,nTempMagn
                    If (smagn) Then
                       Call MSUM ( nCenter, Sex(:,iT), Zex(iT),
     &                             SLT(:, :,iT), ZLT(:,iT),
     &                             SRT(:, :,iT), ZRT(:,iT), iopt,
     &                             ST(:,iT),  ZT(iT)  )
                    End If

                    Call MSUM ( nCenter, Mex(:,iT), Zex(iT),
     &                          MLT(:, :,iT), ZLT(:,iT),
     &                          MRT(:, :,iT), ZRT(:,iT), iopt,
     &                           MT(:,iT),  ZT(iT)  )
                  End Do

               Else ! (m_accurate)

                  ! add the contribution from local excited states
                  ! using the approximate X*H expression:
                  Do iT=1,nTempMagn
                     Do isite=1,nCenter
                        Do l=1,3

                           MRT(isite,l,iT)=(
     &                                XRM(isite,iT,l,1)*H(iH)*dX(iPl,iM)
     &                               +XRM(isite,iT,l,2)*H(iH)*dY(iPl,iM)
     &                               +XRM(isite,iT,l,3)*H(iH)*dZ(iPl,iM)
     &                                ) /cm3tomB

                           MLT(isite,l,iT)=(
     &                                XLM(isite,iT,l,1)*H(iH)*dX(iPl,iM)
     &                               +XLM(isite,iT,l,2)*H(iH)*dY(iPl,iM)
     &                               +XLM(isite,iT,l,3)*H(iH)*dZ(iPl,iM)
     &                                ) /cm3tomB

                        End Do ! l
                     End Do ! isite

                     If (smagn) Then
                        Call MSUM ( nCenter, Sex(:,iT), Zex(iT),
     &                              SLT(:, :,iT), ZLM(:,iT),
     &                              SRT(:, :,iT), ZRM(:,iT), iopt,
     &                              ST(:,iT),  ZT(iT)  )
                     End If
                     Call MSUM ( nCenter, Mex(:,iT), Zex(iT),
     &                           MLT(:, :,iT), ZLM(:,iT),
     &                           MRT(:, :,iT), ZRM(:,iT), iopt,
     &                           MT(:,iT),  ZT(iT)  )
                  End Do ! iT
               End If ! (m_accurate)


               !  at this point we have MT and ST computed
               !  in the direction of applied field
               ! compute the M and S torque
               Do iT=1,nTempMagn
                  tx(iPl,iM,iH,iT) = MT(2,iT)*dZ(iPl,iM)*H(iH)
     &                              -MT(3,iT)*dY(iPl,iM)*H(iH)

                  ty(iPl,iM,iH,iT) = MT(3,iT)*dX(iPl,iM)*H(iH)
     &                              -MT(1,iT)*dZ(iPl,iM)*H(iH)

                  tz(iPl,iM,iH,iT) = MT(1,iT)*dY(iPl,iM)*H(iH)
     &                              -MT(2,iT)*dX(iPl,iM)*H(iH)

                  If (smagn) Then
                     sx(iPl,iM,iH,iT) = ST(2,iT)*dZ(iPl,iM)*H(iH)
     &                                 -ST(3,iT)*dY(iPl,iM)*H(iH)

                     sy(iPl,iM,iH,iT) = ST(3,iT)*dX(iPl,iM)*H(iH)
     &                                 -ST(1,iT)*dZ(iPl,iM)*H(iH)

                     sz(iPl,iM,iH,iT) = ST(1,iT)*dY(iPl,iM)*H(iH)
     &                                 -ST(2,iT)*dX(iPl,iM)*H(iH)
                 End If
               End Do !iT

c ///  closing the loops over field strengths and directions
          End Do ! iL
        End Do ! iM
      End Do ! iH
c -------------------------------------------------------------------
C   WRITING SOME OF THE OUTPUT....
C -------------------------------------------------------------------
      Write(6,*)
      Write(6,'(25X,A)') 'ANGULAR DEPENDENCE OF THE MAGNETIZATION '//
     &                   'TORQUE'
      Write(6,'(30X,A)') '(Units of torque: [energy, cm-1])'
      Write(6,*)

        Write(6,'(5x,A)') 'Orientation of the applied magnetic'//
     & ' field employed:'
        Write(6,'(10A)') '--------|',
     & ('--------------------------------------------|',i=1,3)

        Write(6,'(2x,A,3(10x,A))') 'Angle |',
     & 'rotation in the YZ plane          |',
     & 'rotation in the XZ plane          |',
     & 'rotation in the XY plane          |'
        Write(6,'(10A)') '--------|',
     & ('--- proj X ---|','--- proj Y ---|','--- proj Z ---|',i=1,3)
         Do iM=1,AngPoints
         Write(6,'(F7.3,1x,A,3(F13.10,1x,A,F13.10,1x,A,F13.10,1x,A))')
     &   Ang(iM),'|',
     &   ( dX(iPl,iM),' ',dY(iPl,iM),' ',dZ(iPl,iM),'|', iPl=1,3 )
           End Do
        Write(6,'(10A)') '--------|',
     & ('--------------------------------------------|',i=1,3)

      Do iH=1,nH
        Do iT=1,nTempMagn
        Write(6,*)
        Write(6,'(10A)') '--------|',
     & ('--------------------------------------------|',i=1,3)
        Write(6,'(A,F9.4,A)') 'Magnetic field strength = ',H(iH),
     &                        ' Tesla'
        Write(6,'(12x,A,F9.4,A)') 'Temperature = ',TempMagn(iT),
     &                            ' Kelvin'
        Write(6,'(10A)') '--------|',
     & ('--------------------------------------------|',i=1,3)

        Write(6,'(2x,A,3(10x,A))') 'Angle |',
     & 'rotation in the YZ plane          |',
     & 'rotation in the XZ plane          |',
     & 'rotation in the XY plane          |'
        Write(6,'(10A)') '--------|',
     & ('-- torque X --|','-- torque Y --|','-- torque Z --|',i=1,3)
         Do iM=1,AngPoints
         Write(6,'(F7.3,1x,A,3(E13.6,1x,A,E13.6,1x,A,E13.6,1x,A))')
     &   Ang(iM),'|',
     &   ( tx(iPl,iM,iH,iT),' ',
     &     ty(iPl,iM,iH,iT),' ',
     &     tz(iPl,iM,iH,iT),'|', iPl=1,3 )
           End Do
        Write(6,'(10A)') '--------|',
     & ('--------------------------------------------|',i=1,3)
        End Do !iT
      End Do !iH

C -------------------------------------------------------------------
      If(smagn) Then
      Write(6,*)
      Write(6,'(25X,A)') 'ANGULAR DEPENDENCE OF THE SPIN '//
     &                   'MAGNETIZATION TORQUE'
      Write(6,'(30X,A)') '(Units of torque: [energy, cm-1])'
      Write(6,*)

        Write(6,'(5x,A)') 'Orientation of the applied magnetic'//
     & ' field employed:'
        Write(6,'(10A)') '--------|',
     & ('--------------------------------------------|',i=1,3)

        Write(6,'(2x,A,3(10x,A))') 'Angle |',
     & 'rotation in the YZ plane          |',
     & 'rotation in the XZ plane          |',
     & 'rotation in the XY plane          |'
        Write(6,'(10A)') '--------|',
     & ('--- proj X ---|','--- proj Y ---|','--- proj Z ---|',i=1,3)
         Do iM=1,AngPoints
         Write(6,'(F7.3,1x,A,3(F13.10,1x,A,F13.10,1x,A,F13.10,1x,A))')
     &   Ang(iM),'|',
     &   ( dX(iPl,iM),' ',dY(iPl,iM),' ',dZ(iPl,iM),'|', iPl=1,3 )
           End Do
        Write(6,'(10A)') '--------|',
     & ('--------------------------------------------|',i=1,3)

      Do iH=1,nH
        Do iT=1,nTempMagn
        Write(6,*)
        Write(6,'(10A)') '--------|',
     & ('--------------------------------------------|',i=1,3)
        Write(6,'(A,F9.4,A)') 'Magnetic field strength = ',H(iH),
     &                        ' Tesla'
        Write(6,'(12x,A,F9.4,A)') 'Temperature = ',TempMagn(iT),
     &                            ' Kelvin'
        Write(6,'(10A)') '--------|',
     & ('--------------------------------------------|',i=1,3)

        Write(6,'(2x,A,3(10x,A))') 'Angle |',
     & 'rotation in the YZ plane          |',
     & 'rotation in the XZ plane          |',
     & 'rotation in the XY plane          |'
        Write(6,'(10A)') '--------|',
     & ('Spin torque X |','Spin torque Y |','Spin torque Z |',i=1,3)
         Do iM=1,AngPoints
         Write(6,'(F7.3,1x,A,3(E13.6,1x,A,E13.6,1x,A,E13.6,1x,A))')
     &   Ang(iM),'|',
     &   ( sx(iPl,iM,iH,iT),' ',
     &     sy(iPl,iM,iH,iT),' ',
     &     sz(iPl,iM,iH,iT),'|', iPl=1,3 )
           End Do
        Write(6,'(10A)') '--------|',
     & ('--------------------------------------------|',i=1,3)
        End Do !iT
      End Do !iH

      End If !smagn

! 199  Continue
!-----------------------------------------------------------------------
! Deallocate memory for this calculation:
      If(nM>0) Then
         Call mma_deallocate(Wex)
      End If

      If(nTempMagn>0) Then
         Call mma_deallocate(Zex)
         Call mma_deallocate(SEX)
         Call mma_deallocate(MEX)
         Call mma_deallocate(ZT)
         Call mma_deallocate(ST)
         Call mma_deallocate(MT)

         If(nneq>0) Then
            Call mma_deallocate(ZL)
            Call mma_deallocate(SL)
            Call mma_deallocate(ML)
            Call mma_deallocate(ZR)
            Call mma_deallocate(SR)
            Call mma_deallocate(MR)
         End If

         If(nCenter>0) Then
            Call mma_deallocate(ZRT)
            Call mma_deallocate(ZLT)
            Call mma_deallocate(MRT)
            Call mma_deallocate(MLT)
            Call mma_deallocate(SRT)
            Call mma_deallocate(SLT)
         End If

         If( (nPlanes>0).and.(AngPoints>0).and.(nH>0) ) Then
            Call mma_deallocate(tx)
            Call mma_deallocate(ty)
            Call mma_deallocate(tz)
            Call mma_deallocate(sx)
            Call mma_deallocate(sy)
            Call mma_deallocate(sz)
         End If
      End If

      If((nLoc>0).and.(nneq>0)) Then
         Call mma_deallocate(WL)
         Call mma_deallocate(WR)
      End If

      If(AngPoints>0) Then
         Call mma_deallocate(Ang)
         If(nPlanes>0) Then
            Call mma_deallocate(dX)
            Call mma_deallocate(dY)
            Call mma_deallocate(dZ)
         End If
      End If

      If(nH>0) Then
         Call mma_deallocate(H)
      End If

      If(dbg) Write(6,*) 'TORQ: allocated memory was sucessfully '//
     &                   'deallocated'

      Call qExit('PA_torq')
      Return
      End
