!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine torque_pa(nneq,nCenter,neq,neqv,nLoc,exch,nTempMagn,nH,nM,AngPoints,nexch,iopt,nss,mem,smagn,m_paranoid,m_accurate, &
                     TempMagn,w,hmin,hmax,dltH0,EM,zJ,THRS,hexp,dipexch,s_exch,dipso,s_so,eso,hinput,r_rot,XLM,ZLM,XRM,ZRM)

use Constants, only: Zero, Ten, mBohr, rNAVO
use Definitions, only: wp, u6

implicit none
#include "mgrid.fh"
#include "stdalloc.fh"
!-----------------------------------------------------------------------
! INPUT VARIABLES:
!-----------------------------------------------------------------------
integer, intent(in) :: nneq, neq(nneq), neqv, nLoc, nexch(nneq), nH, nTempMagn, nM, exch, nCenter, AngPoints, nss(nneq), iopt, mem
real(kind=8), intent(in) :: TempMagn(nTempMagn), hmin, hmax, dltH0, EM, zJ, THRS, Hexp(nH)
logical, intent(in) :: m_paranoid, m_accurate, smagn, hinput
! correction to M from the local excited states:
real(kind=8), intent(in) :: XLM(nCenter,nTempMagn,3,3)
real(kind=8), intent(in) :: ZLM(nCenter,nTempMagn)
real(kind=8), intent(in) :: XRM(nCenter,nTempMagn,3,3)
real(kind=8), intent(in) :: ZRM(nCenter,nTempMagn)
! rotation matrices for equivalent sites:
real(kind=8), intent(in) :: R_ROT(nneq,neqv,3,3)
! exchange spectum:
! exchange energies printed out in the previous part
real(kind=8), intent(in) :: W(exch)
! local spin-orbit spectum:
! spin-orbit energies from ANISO files
real(kind=8), intent(in) :: ESO(nneq,nLoc)
! magnetic and spin moments (i.e. the BIG matrices):
complex(kind=8), intent(in) :: DIPEXCH(3,EXCH,EXCH)
complex(kind=8), intent(in) :: S_EXCH(3,EXCH,EXCH)
complex(kind=8), intent(in) :: dipso(nneq,3,nLoc,nLoc)
complex(kind=8), intent(in) :: s_so(nneq,3,nLoc,nLoc)
!-----------------------------------------------------------------------
! LOCAL VARIABLES:
!-----------------------------------------------------------------------
! exchange data:
real(kind=8), allocatable :: WEX(:)   ! WEX(NM)                ! Zeeman exchange energies
real(kind=8), allocatable :: ZEX(:)   ! ZEX(nTempMagn)         ! exchange statistical sum, Boltzmann distribution
real(kind=8), allocatable :: SEX(:,:) ! SEX(3,nTempMagn)       ! spin magnetisation, from the exchange block
real(kind=8), allocatable :: MEX(:,:) ! MEX(3,nTempMagn)       ! magnetisation, form the exchange block
! data for individual sites (all states):
real(kind=8), allocatable :: ZL(:,:)   ! ZL(nneq,nTempMagn)     ! local statistical sum, Boltzmann distribution
real(kind=8), allocatable :: WL(:,:)   ! WL(nneq,nLoc)          ! Zeeman local energies
real(kind=8), allocatable :: SL(:,:,:) ! SL(nneq,3,nTempMagn)   ! spin magnetisation, from the local sites, using ALL states
real(kind=8), allocatable :: ML(:,:,:) ! ML(nneq,3,nTempMagn)   ! magnetisation, from local sites, using ALL states
! data for individual sites (only states that enter exchange):
real(kind=8), allocatable :: ZR(:,:)   ! ZR(nneq,nTempMagn)     ! local statistical sum, Boltzmann distribution, using only NEXCH states
real(kind=8), allocatable :: WR(:,:)   ! WR(nneq,nLoc)          ! Zeeman local reduced energies, using only NEXCH states
real(kind=8), allocatable :: SR(:,:,:) ! SR(nneq,3,nTempMagn)   ! spin magnetisation, from the local sites, using only NEXCH states
real(kind=8), allocatable :: MR(:,:,:) ! MR(nneq,3,nTempMagn)   ! magnetisation, from local sites, using only NEXCH states
! total vectors in general coordinate system:
real(kind=8), allocatable :: ZRT(:,:)   ! ZRT(nCenter,nTempMagn)
real(kind=8), allocatable :: ZLT(:,:)   ! ZLT(nCenter,nTempMagn)
real(kind=8), allocatable :: MRT(:,:,:) ! MRT(nCenter,3,nTempMagn)
real(kind=8), allocatable :: MLT(:,:,:) ! MLT(nCenter,3,nTempMagn)
real(kind=8), allocatable :: SRT(:,:,:) ! SRT(nCenter,3,nTempMagn)
real(kind=8), allocatable :: SLT(:,:,:) ! SLT(nCenter,3,nTempMagn)
! data for total system:
real(kind=8), allocatable :: ZT(:)   ! ZT(nTempMagn)        ! total statistical sum, Boltzmann distribution
real(kind=8), allocatable :: ST(:,:) ! ST(3,nTempMagn)      ! total spin magnetisation,
real(kind=8), allocatable :: MT(:,:) ! MT(3,nTempMagn)      ! total magnetisation
! magnetic field strength and orientation data:
integer :: nPlanes
real(kind=8) :: dlth
parameter(nPlanes=3)
real(kind=8), allocatable :: H(:)    ! H(nH)
real(kind=8), allocatable :: dX(:,:) ! dX(nPlanes,AngPoints)
real(kind=8), allocatable :: dY(:,:) ! dY(nPlanes,AngPoints)
real(kind=8), allocatable :: dZ(:,:) ! dZ(nPlanes,AngPoints)
real(kind=8), allocatable :: Ang(:)  ! Ang(AngPoints)
! magnetic torque
real(kind=8), allocatable :: tx(:,:,:,:) ! tx(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque, X
real(kind=8), allocatable :: ty(:,:,:,:) ! ty(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque, Y
real(kind=8), allocatable :: tz(:,:,:,:) ! tz(nPlanes,AngPoints,nH,nTempMagn) ! magnetization torque, Z
real(kind=8), allocatable :: sx(:,:,:,:) ! sx(nPlanes,AngPoints,nH,nTempMagn) ! spin magnetization torque, X
real(kind=8), allocatable :: sy(:,:,:,:) ! sy(nPlanes,AngPoints,nH,nTempMagn) ! spin magnetization torque, Y
real(kind=8), allocatable :: sz(:,:,:,:) ! sz(nPlanes,AngPoints,nH,nTempMagn) ! spin magnetization torque, Z
integer :: mem_local, RtoB
! local data:
integer :: IM, I, it
integer :: J, IH, k, isite, l, n, iPl
logical :: DBG
real(kind=8), parameter :: cm3tomB = rNAVO*mBohr/Ten ! in cm3 * mol-1 * T

DBG = .false.

write(u6,*)
write(u6,'(100A)') (('%'),J=1,96)
write(u6,'(20X,A)') 'ANGULAR DEPENDENCE OF THE MAGNETIZATION TORQUE'
write(u6,'(100A)') (('%'),J=1,96)
write(u6,*)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
write(u6,'(2X,A,i3,A)') 'Magnetization torque is calculated for the ',NH,' field points, in the field Domain:'
write(u6,'(2X,F4.1,1x,a,1X,F4.1,a,30(F6.3,a))') HMIN,'--',HMAX,' T., at the following temperatures:'
do i=11,nTempMagn,10
  j = min(nTempMagn,i+9)
  write(u6,'(17x,10(F9.3,A))') (TempMagn(k),' K.;',k=i,j)
end do
write(u6,'(2x,A,i3,A)') 'Angular dependence of the magnetization torque is computed for ',AngPoints,' angular points distributed'
write(u6,'(2x,A)') 'in the domain 0-180 deg.'
write(u6,'(2X,A)') 'The cut-off energy for the exact diagonalization of the Zeeman Hamiltonian is:'
write(u6,'(2x,a,F15.9,A)') 'E = ',EM,' cm(-1).'
if (NM < 10) then
  write(u6,'(2X,A,i2,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 10) .and. (NM < 100)) then
  write(u6,'(2X,A,i3,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 100) .and. (NM < 1000)) then
  write(u6,'(2X,A,i4,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
else if ((NM >= 1000) .and. (NM < 10000)) then
  write(u6,'(2X,A,i5,a)') 'The exact diagonalization of the Zeeman Hamiltonian included ',NM,' exchange states.'
end if
if (m_accurate) then
  write(u6,'(2x,A)') 'The contribution of local excited states is computed exactly.'
else
  write(u6,'(2x,A)') 'The contribution of local excited states is computed approximately (by using the susceptibility data). '
end if
!-----------------------------------------------------------------------
if (dbg) write(u6,*) 'nM       = ',nM
if (dbg) write(u6,*) 'nTempMagn= ',nTempMagn
if (dbg) write(u6,*) 'nneq     = ',nneq
if (dbg) write(u6,*) 'nCenter  = ',nCenter
if (dbg) write(u6,*) 'AngPoints= ',AngPoints
if (dbg) write(u6,*) 'nPlanes  = ',nPlanes
if (dbg) write(u6,*) 'nH       = ',nH
if (dbg) write(u6,*) 'nLoc     = ',nLoc
if (dbg) write(u6,*) 'nexch()  = ',(nexch(i),i=1,nneq)
if (dbg) write(u6,*) 'neq()    = ',(neq(i),i=1,nneq)
if (dbg) write(u6,*) 'exch     = ',exch
if (dbg) write(u6,*) 'iopt     = ',iopt
if (dbg) write(u6,*) 'neqv     = ',neqv
if (dbg) write(u6,*) 'nss()    = ',(nss(i),i=1,nneq)
if (dbg) write(u6,*) 'W()      = ',(W(i),i=1,exch)
if (dbg) write(u6,*) 'zJ       = ',zJ
if (dbg) write(u6,*) 'EM       = ',EM
if (dbg) write(u6,*) 'm_paranoi= ',m_paranoid
if (dbg) write(u6,*) 'm_accurat= ',m_accurate
if (dbg) write(u6,*) 'smagn    = ',smagn
if (dbg) write(u6,*) 'hinput   = ',hinput

! Allocate memory for this calculation:
mem_local = 0
RtoB = 8
! Zeeman exchange energy spectrum
call mma_allocate(Wex,nM,'Wex')
call dcopy_(nM,[Zero],0,Wex,1)
mem_local = mem_local+nM*RtoB

if (dbg) write(u6,*) 'mem_local 1 = ',mem_local

! exchange statistical sum, Boltzmann distribution
call mma_allocate(Zex,nTempMagn,'Zex')
call dcopy_(nTempMagn,[Zero],0,Zex,1)
mem_local = mem_local+nTempMagn*RtoB
! spin magnetisation, from the exchange block
call mma_allocate(SEX,3,nTempMagn,'SEX')
call dcopy_(3*nTempMagn,[Zero],0,SEX,1)
mem_local = mem_local+3*nTempMagn*RtoB
! magnetisation, from the exchange block
call mma_allocate(MEX,3,nTempMagn,'MEX')
call dcopy_(3*nTempMagn,[Zero],0,MEX,1)
mem_local = mem_local+3*nTempMagn*RtoB

! total statistical sum, Boltzmann distribution
call mma_allocate(ZT,nTempMagn,'ZT')
call dcopy_(nTempMagn,[Zero],0,ZT,1)
mem_local = mem_local+nTempMagn*RtoB
! total spin magnetisation
call mma_allocate(ST,3,nTempMagn,'ST')
call dcopy_(3*nTempMagn,[Zero],0,ST,1)
mem_local = mem_local+3*nTempMagn*RtoB
! total magnetisation
call mma_allocate(MT,3,nTempMagn,'MT')
call dcopy_(3*nTempMagn,[Zero],0,MT,1)
mem_local = mem_local+3*nTempMagn*RtoB

if (dbg) write(u6,*) 'mem_local 2 = ',mem_local
! local statistical sum, Boltzmann distribution
call mma_allocate(ZL,nneq,nTempMagn,'ZL')
call dcopy_(nneq*nTempMagn,[Zero],0,ZL,1)
mem_local = mem_local+nneq*nTempMagn*RtoB
! spin magnetisation, from the local sites, using ALL states
call mma_allocate(SL,nneq,3,nTempMagn,'SL')
call dcopy_(nneq*3*nTempMagn,[Zero],0,SL,1)
mem_local = mem_local+nneq*3*nTempMagn*RtoB
! magnetisation, from local sites, using ALL states
call mma_allocate(ML,nneq,3,nTempMagn,'ML')
call dcopy_(nneq*3*nTempMagn,[Zero],0,ML,1)
mem_local = mem_local+nneq*3*nTempMagn*RtoB

! local statistical sum, Boltzmann distribution
call mma_allocate(ZR,nneq,nTempMagn,'ZR')
call dcopy_(nneq*nTempMagn,[Zero],0,ZR,1)
mem_local = mem_local+nneq*nTempMagn*RtoB
! spin magnetisation, from the local sites, using only NEXCH states
call mma_allocate(SR,nneq,3,nTempMagn,'SR')
call dcopy_(nneq*3*nTempMagn,[Zero],0,SR,1)
mem_local = mem_local+nneq*3*nTempMagn*RtoB
! magnetisation, from the local sites, using only NEXCH states
call mma_allocate(MR,nneq,3,nTempMagn,'MR')
call dcopy_(nneq*3*nTempMagn,[Zero],0,MR,1)
mem_local = mem_local+nneq*3*nTempMagn*RtoB

if (dbg) write(u6,*) 'mem_local 3 = ',mem_local
! ZRT(nCenter,nTempMagn)
call mma_allocate(ZRT,nCenter,nTempMagn,'ZRT')
call dcopy_(nCenter*nTempMagn,[Zero],0,ZRT,1)
mem_local = mem_local+nCenter*nTempMagn*RtoB
! ZLT(nCenter,nTempMagn)
call mma_allocate(ZLT,nCenter,nTempMagn,'ZLT')
call dcopy_(nCenter*nTempMagn,[Zero],0,ZLT,1)
mem_local = mem_local+nCenter*nTempMagn*RtoB
! MRT(nCenter,3,nTempMagn)
call mma_allocate(MRT,nCenter,3,nTempMagn,'MRT')
call dcopy_(nCenter*3*nTempMagn,[Zero],0,MRT,1)
mem_local = mem_local+nCenter*3*nTempMagn*RtoB
! MLT(nCenter,3,nTempMagn)
call mma_allocate(MLT,nCenter,3,nTempMagn,'MLT')
call dcopy_(nCenter*3*nTempMagn,[Zero],0,MLT,1)
mem_local = mem_local+nCenter*3*nTempMagn*RtoB
! SRT(nCenter,3,nTempMagn)
call mma_allocate(SRT,nCenter,3,nTempMagn,'SRT')
call dcopy_(nCenter*3*nTempMagn,[Zero],0,SRT,1)
mem_local = mem_local+nCenter*3*nTempMagn*RtoB
! SLT(nCenter,3,nTempMagn)
call mma_allocate(SLT,nCenter,3,nTempMagn,'SLT')
call dcopy_(nCenter*3*nTempMagn,[Zero],0,SLT,1)
mem_local = mem_local+nCenter*3*nTempMagn*RtoB

if (dbg) write(u6,*) 'mem_local 4 = ',mem_local
! magnetic torque
call mma_allocate(tx,nPlanes,AngPoints,nH,nTempMagn,'tx')
call mma_allocate(ty,nPlanes,AngPoints,nH,nTempMagn,'ty')
call mma_allocate(tz,nPlanes,AngPoints,nH,nTempMagn,'tz')
call mma_allocate(sx,nPlanes,AngPoints,nH,nTempMagn,'sx')
call mma_allocate(sy,nPlanes,AngPoints,nH,nTempMagn,'sy')
call mma_allocate(sz,nPlanes,AngPoints,nH,nTempMagn,'sz')
call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[Zero],0,tx,1)
call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[Zero],0,ty,1)
call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[Zero],0,tz,1)
call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[Zero],0,sx,1)
call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[Zero],0,sy,1)
call dcopy_(nPlanes*AngPoints*nH*nTempMagn,[Zero],0,sz,1)
mem_local = mem_local+6*nPlanes*AngPoints*nH*nTempMagn*RtoB
if (dbg) write(u6,*) 'mem_local 5 = ',mem_local

! Zeeman local energies
call mma_allocate(WL,nneq,nLoc,'WL')
call dcopy_(nneq*nLoc,[Zero],0,WL,1)
mem_local = mem_local+nneq*nLoc*RtoB
! Zeeman local reduced energies, using only NEXCH states
call mma_allocate(WR,nneq,nLoc,'WR')
call dcopy_(nneq*nLoc,[Zero],0,WR,1)
mem_local = mem_local+nneq*nLoc*RtoB
if (dbg) write(u6,*) 'mem_local 6 = ',mem_local

call mma_allocate(Ang,AngPoints,'Ang')
call dcopy_(AngPoints,[Zero],0,Ang,1)
mem_local = mem_local+AngPoints*RtoB
call mma_allocate(dX,nPlanes,AngPoints,'dX')
call mma_allocate(dY,nPlanes,AngPoints,'dY')
call mma_allocate(dZ,nPlanes,AngPoints,'dZ')
call dcopy_(nPlanes*AngPoints,[Zero],0,dX,1)
call dcopy_(nPlanes*AngPoints,[Zero],0,dY,1)
call dcopy_(nPlanes*AngPoints,[Zero],0,dZ,1)
mem_local = mem_local+3*nPlanes*AngPoints*RtoB
if (dbg) write(u6,*) 'mem_local 7 = ',mem_local

call mma_allocate(H,nH,'H')
call dcopy_(nH,[Zero],0,H,1)
mem_local = mem_local+nH*RtoB

if (dbg) write(u6,*) 'TORQ:  memory allocated (local):'
if (dbg) write(u6,*) 'mem_local=',mem_local
if (dbg) write(u6,*) 'TORQ:  memory allocated (total):'
if (dbg) write(u6,*) 'mem_total=',mem+mem_local

!-----------------------------------------------------------------------
! set up the field points:
if (HINPUT) then
  do iH=1,nH
    H(iH) = HEXP(iH)
    if (H(iH) == Zero) then
      H(iH) = 0.0001_wp
    end if
  end do
else
  DLTH = (HMAX-HMIN)/real(NH-1,kind=wp)
  do iH=1,nH
    if (iH == 1) then
      H(IH) = HMIN+dltH0
    else
      H(IH) = HMIN+DLTH*real(IH-1,kind=wp)
    end if
    if (H(iH) == Zero) then
      H(iH) = 0.0001_wp
    end if
  end do
end if
!-----------------------------------------------------------------------

do IH=1,NH
  ! ///  opening the loop over field points:
  call dcopy_(nPlanes*AngPoints,[Zero],0,dX,1)
  call dcopy_(nPlanes*AngPoints,[Zero],0,dY,1)
  call dcopy_(nPlanes*AngPoints,[Zero],0,dZ,1)
  do iPl=1,nPlanes
    ! ///  opening the loop over different planes of rotation of the applied magnetic field:
    !  loop over various angular grids  (iPl = 1, 2 or 3)
    !  iPl=1 (X) , angular grid in the plane YZ ( half of the plane , 0-180 degrees)
    !  iPl=2 (Y) , angular grid in the plane XZ ( half of the plane , 0-180 degrees)
    !  iPl=3 (Z) , angular grid in the plane XY ( half of the plane , 0-180 degrees)
    call dcopy_(AngPoints,[Zero],0,Ang,1)

    call hdir2(AngPoints,iPl,dX(iPl,:),dY(iPl,:),dZ(iPl,:),Ang,4)

    if (dbg) write(u6,*) 'iPl, dX, dY, dZ=',iPl,dX(iPl,:),dY(iPl,:),dZ(iPl,:),Ang
    do IM=1,AngPoints
      ! ///  opening the loop over different directions of the magnetic field
      call dcopy_(nM,[Zero],0,Wex,1)
      call dcopy_(nneq*nLoc,[Zero],0,WL,1)
      call dcopy_(nneq*nLoc,[Zero],0,WR,1)
      call dcopy_(nTempMagn,[Zero],0,Zex,1)
      call dcopy_(nneq*nTempMagn,[Zero],0,ZL,1)
      call dcopy_(nneq*nTempMagn,[Zero],0,ZR,1)
      call dcopy_(3*nTempMagn,[Zero],0,SEX,1)
      call dcopy_(3*nTempMagn,[Zero],0,MEX,1)
      call dcopy_(nneq*3*nTempMagn,[Zero],0,SL,1)
      call dcopy_(nneq*3*nTempMagn,[Zero],0,ML,1)
      call dcopy_(nneq*3*nTempMagn,[Zero],0,SR,1)
      call dcopy_(nneq*3*nTempMagn,[Zero],0,MR,1)
      call dcopy_(3*nTempMagn,[Zero],0,MT,1)
      call dcopy_(3*nTempMagn,[Zero],0,ST,1)
      call dcopy_(nTempMagn,[Zero],0,ZT,1)

      ! exchange magnetization:
      call MAGN(EXCH,NM,dX(iPl,iM),dY(iPl,iM),dZ(iPl,iM),H(iH),W,zJ,THRS,DIPEXCH,S_EXCH,nTempMagn,TempMagn,smagn,Wex,Zex,Sex,Mex, &
                m_paranoid,DBG)

      if (iPl == 2) then
        if (DBG) write(u6,'(A,I3,1x,F8.4,2x, 3F19.14,2x,3F19.14)') 'MEX: iM,',iM,H(iH),(Mex(l,1),l=1,3),dX(iPl,iM),dY(iPl,iM), &
                                                                   dZ(iPl,iM)
      end if
      !iT = 1

      !write(u6,'(F10.4,1x,2I3,3F18.14,3x,3F20.14,2x,F20.14)') H(iH),iL,iM,dX(iM),dY(iM),dZ(iM),(Mex(j,iT),j=1,3), &
      !                                                        Mex(1,iT)*dX(iM)+Mex(2,iT)*dY(iM)+Mex(3,iT)*dZ(iM)
      ! compute local magnetizations:
      if (m_accurate) then
        do i=1,nneq
          ! all states:
          if (NSS(i) > NEXCH(i)) then
            ! this check is meant to avoid the unnecessary
            ! computation, in cases when no local excited
            ! states are present
            call MAGN(NSS(i),NEXCH(i),dX(iPl,iM),dY(iPl,iM),dZ(iPl,iM),H(iH),ESO(i,1:NSS(i)),zJ,THRS,DIPSO(i,:,:,:),S_SO(i,:,:,:), &
                      nTempMagn,TempMagn,smagn,WL(i,:),ZL(i,:),SL(i,:,:),ML(i,:,:),m_paranoid,DBG)
            ! only local "exchange states":
            call MAGN(NEXCH(i),NEXCH(i),dX(iPl,iM),dY(iPl,iM),dZ(iPl,iM),H(iH),ESO(i,1:NEXCH(i)),zJ,THRS, &
                      DIPSO(i,1:3,1:NEXCH(i),1:NEXCH(i)),S_SO(i,1:3,1:NEXCH(i),1:NEXCH(i)),nTempMagn,TempMagn,smagn, &
                      WR(i,1:Nexch(i)),ZR(i,:),SR(i,:,:),MR(i,:,:),m_paranoid,DBG)
          end if
        end do
        ! expand the basis and rotate local vectors to the
        ! general coordinate system:
        call dcopy_(nCenter*nTempMagn,[Zero],0,ZRT,1)
        call dcopy_(nCenter*nTempMagn,[Zero],0,ZLT,1)
        call dcopy_(3*nCenter*nTempMagn,[Zero],0,MRT,1)
        call dcopy_(3*nCenter*nTempMagn,[Zero],0,MLT,1)
        call dcopy_(3*nCenter*nTempMagn,[Zero],0,SRT,1)
        call dcopy_(3*nCenter*nTempMagn,[Zero],0,SLT,1)

        isite = 0
        do i=1,NNEQ
          do j=1,NEQ(i)
            isite = isite+1
            ! statistical distributions
            do iT=1,nTempMagn
              ZLT(isite,iT) = ZL(i,iT)
              ZRT(isite,iT) = ZR(i,iT)
            end do
            ! magnetizations:
            ! use R_rot matrices, which have determinant +1.
            !  >> note that  R_lg matrices may have arbitrary
            !  >> sign of the determinant.
            do iT=1,nTempMagn
              do l=1,3
                do n=1,3
                  MLT(isite,l,iT) = MLT(isite,l,iT)+r_rot(i,j,l,n)*ML(i,n,iT)

                  SLT(isite,l,iT) = SLT(isite,l,iT)+r_rot(i,j,l,n)*SL(i,n,iT)

                  MRT(isite,l,iT) = MRT(isite,l,iT)+r_rot(i,j,l,n)*MR(i,n,iT)

                  SRT(isite,l,iT) = SRT(isite,l,iT)+r_rot(i,j,l,n)*SR(i,n,iT)
                end do
              end do
            end do

          end do ! j, neq(i)
        end do ! i, nneq
      end if ! m_accurate

      ! compute the total magnetizations according
      ! to the derived formulas:
      if (m_accurate) then
        do iT=1,nTempMagn
          if (smagn) then
            call MSUM(nCenter,Sex(:,iT),Zex(iT),SLT(:,:,iT),ZLT(:,iT),SRT(:,:,iT),ZRT(:,iT),iopt,ST(:,iT),ZT(iT))
          end if

          call MSUM(nCenter,Mex(:,iT),Zex(iT),MLT(:,:,iT),ZLT(:,iT),MRT(:,:,iT),ZRT(:,iT),iopt,MT(:,iT),ZT(iT))
        end do

      else ! (m_accurate)

        ! add the contribution from local excited states
        ! using the approximate X*H expression:
        do iT=1,nTempMagn
          do isite=1,nCenter
            do l=1,3

              MRT(isite,l,iT) = (XRM(isite,iT,l,1)*H(iH)*dX(iPl,iM)+XRM(isite,iT,l,2)*H(iH)*dY(iPl,iM)+ &
                  XRM(isite,iT,l,3)*H(iH)*dZ(iPl,iM))/cm3tomB

              MLT(isite,l,iT) = (XLM(isite,iT,l,1)*H(iH)*dX(iPl,iM)+XLM(isite,iT,l,2)*H(iH)*dY(iPl,iM)+ &
                  XLM(isite,iT,l,3)*H(iH)*dZ(iPl,iM))/cm3tomB

            end do ! l
          end do ! isite

          if (smagn) then
            call MSUM(nCenter,Sex(:,iT),Zex(iT),SLT(:,:,iT),ZLM(:,iT),SRT(:,:,iT),ZRM(:,iT),iopt,ST(:,iT),ZT(iT))
          end if
          call MSUM(nCenter,Mex(:,iT),Zex(iT),MLT(:,:,iT),ZLM(:,iT),MRT(:,:,iT),ZRM(:,iT),iopt,MT(:,iT),ZT(iT))
        end do ! iT
      end if ! (m_accurate)

      ! at this point we have MT and ST computed
      ! in the direction of applied field
      ! compute the M and S torque
      do iT=1,nTempMagn
        tx(iPl,iM,iH,iT) = MT(2,iT)*dZ(iPl,iM)*H(iH)-MT(3,iT)*dY(iPl,iM)*H(iH)

        ty(iPl,iM,iH,iT) = MT(3,iT)*dX(iPl,iM)*H(iH)-MT(1,iT)*dZ(iPl,iM)*H(iH)

        tz(iPl,iM,iH,iT) = MT(1,iT)*dY(iPl,iM)*H(iH)-MT(2,iT)*dX(iPl,iM)*H(iH)

        if (smagn) then
          sx(iPl,iM,iH,iT) = ST(2,iT)*dZ(iPl,iM)*H(iH)-ST(3,iT)*dY(iPl,iM)*H(iH)

          sy(iPl,iM,iH,iT) = ST(3,iT)*dX(iPl,iM)*H(iH)-ST(1,iT)*dZ(iPl,iM)*H(iH)

          sz(iPl,iM,iH,iT) = ST(1,iT)*dY(iPl,iM)*H(iH)-ST(2,iT)*dX(iPl,iM)*H(iH)
        end if
      end do !iT

      ! ///  closing the loops over field strengths and directions
    end do ! iL
  end do ! iM
end do ! iH
! ----------------------------------------------------------------------
! WRITING SOME OF THE OUTPUT....
! ----------------------------------------------------------------------
write(u6,*)
write(u6,'(25X,A)') 'ANGULAR DEPENDENCE OF THE MAGNETIZATION TORQUE'
write(u6,'(30X,A)') '(Units of torque: [energy, cm-1])'
write(u6,*)

write(u6,'(5x,A)') 'Orientation of the applied magnetic field employed:'
write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

write(u6,'(2x,A,3(10x,A))') 'Angle |','rotation in the YZ plane          |','rotation in the XZ plane          |', &
                            'rotation in the XY plane          |'
write(u6,'(10A)') '--------|',('--- proj X ---|','--- proj Y ---|','--- proj Z ---|',i=1,3)
do iM=1,AngPoints
  write(u6,'(F7.3,1x,A,3(F13.10,1x,A,F13.10,1x,A,F13.10,1x,A))') Ang(iM),'|',(dX(iPl,iM),' ',dY(iPl,iM),' ',dZ(iPl,iM),'|',iPl=1,3)
end do
write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

do iH=1,nH
  do iT=1,nTempMagn
    write(u6,*)
    write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)
    write(u6,'(A,F9.4,A)') 'Magnetic field strength = ',H(iH),' tesla'
    write(u6,'(12x,A,F9.4,A)') 'Temperature = ',TempMagn(iT),' kelvin'
    write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

    write(u6,'(2x,A,3(10x,A))') 'Angle |','rotation in the YZ plane          |','rotation in the XZ plane          |', &
                                'rotation in the XY plane          |'
    write(u6,'(10A)') '--------|',('-- torque X --|','-- torque Y --|','-- torque Z --|',i=1,3)
    do iM=1,AngPoints
      write(u6,'(F7.3,1x,A,3(ES13.6,1x,A,ES13.6,1x,A,ES13.6,1x,A))') Ang(iM),'|',(tx(iPl,iM,iH,iT),' ',ty(iPl,iM,iH,iT),' ', &
                                                                     tz(iPl,iM,iH,iT),'|',iPl=1,3)
    end do
    write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)
  end do !iT
end do !iH

! ----------------------------------------------------------------------
if (smagn) then
  write(u6,*)
  write(u6,'(25X,A)') 'ANGULAR DEPENDENCE OF THE SPIN MAGNETIZATION TORQUE'
  write(u6,'(30X,A)') '(Units of torque: [energy, cm-1])'
  write(u6,*)

  write(u6,'(5x,A)') 'Orientation of the applied magnetic field employed:'
  write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

  write(u6,'(2x,A,3(10x,A))') 'Angle |','rotation in the YZ plane          |','rotation in the XZ plane          |', &
                              'rotation in the XY plane          |'
  write(u6,'(10A)') '--------|',('--- proj X ---|','--- proj Y ---|','--- proj Z ---|',i=1,3)
  do iM=1,AngPoints
    write(u6,'(F7.3,1x,A,3(F13.10,1x,A,F13.10,1x,A,F13.10,1x,A))') Ang(iM),'|', &
                                                                   (dX(iPl,iM),' ',dY(iPl,iM),' ',dZ(iPl,iM),'|',iPl=1,3)
  end do
  write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

  do iH=1,nH
    do iT=1,nTempMagn
      write(u6,*)
      write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)
      write(u6,'(A,F9.4,A)') 'Magnetic field strength = ',H(iH),' tesla'
      write(u6,'(12x,A,F9.4,A)') 'Temperature = ',TempMagn(iT),' kelvin'
      write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)

      write(u6,'(2x,A,3(10x,A))') 'Angle |','rotation in the YZ plane          |','rotation in the XZ plane          |', &
                                  'rotation in the XY plane          |'
      write(u6,'(10A)') '--------|',('Spin torque X |','Spin torque Y |','Spin torque Z |',i=1,3)
      do iM=1,AngPoints
        write(u6,'(F7.3,1x,A,3(ES13.6,1x,A,ES13.6,1x,A,ES13.6,1x,A))') Ang(iM),'|',(sx(iPl,iM,iH,iT),' ',sy(iPl,iM,iH,iT),' ', &
                                                                       sz(iPl,iM,iH,iT),'|',iPl=1,3)
      end do
      write(u6,'(10A)') '--------|',('--------------------------------------------|',i=1,3)
    end do !iT
  end do !iH

end if !smagn

! 199 continue
!-----------------------------------------------------------------------
! Deallocate memory for this calculation:
call mma_deallocate(Wex)

call mma_deallocate(Zex)
call mma_deallocate(SEX)
call mma_deallocate(MEX)
call mma_deallocate(ZT)
call mma_deallocate(ST)
call mma_deallocate(MT)

call mma_deallocate(ZL)
call mma_deallocate(SL)
call mma_deallocate(ML)
call mma_deallocate(ZR)
call mma_deallocate(SR)
call mma_deallocate(MR)

call mma_deallocate(ZRT)
call mma_deallocate(ZLT)
call mma_deallocate(MRT)
call mma_deallocate(MLT)
call mma_deallocate(SRT)
call mma_deallocate(SLT)

call mma_deallocate(tx)
call mma_deallocate(ty)
call mma_deallocate(tz)
call mma_deallocate(sx)
call mma_deallocate(sy)
call mma_deallocate(sz)

call mma_deallocate(WL)
call mma_deallocate(WR)

call mma_deallocate(Ang)
call mma_deallocate(dX)
call mma_deallocate(dY)
call mma_deallocate(dZ)

call mma_deallocate(H)

if (dbg) write(u6,*) 'TORQ: allocated memory was sucessfully deallocated'

return

end subroutine torque_pa
