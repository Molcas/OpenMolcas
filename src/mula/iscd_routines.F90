!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2009, Giovanni Ghigo                                   *
!***********************************************************************

! Contains:
!   ISCD_Rate
!   ISCD_Ene
!   ISCD_LogEVec
!   ISCD_MakeGraphs
!   ISCD_MakenMat
!   ISCD_MakenIncDec
!   ISCD_ReloadNMAT
!   ISCD_FCval

!  InterSystem Crossing rate
!  Author: Giovanni Ghigo
!          Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
!          June 2009

subroutine ISCD_Rate(iPrint,nOsc,max_nOrd,iMaxYes,nYes,dMinWind,lBatch,nIndex,VibWind2,lNMAT0,lNMAT,lNINC,lNDEC,lnTabDim,nnTabDim, &
                     C1,C2,W2,det0,det1,det2,C,W,r01,r02,m_max,n_max,FC00,dRho,nMat,nInc,nDec)

use mula_global, only: hbarcm, maxMax_n, TranDip
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Pi, auTocm
use Definitions, only: wp, iwp, u6, ItoB

implicit none
integer(kind=iwp), intent(in) :: iPrint, nOsc, iMaxYes, nYes, lBatch, nIndex(3,0:maxMax_n), VibWind2(nYes), lNMAT0, lNMAT, lNINC, &
                                 lNDEC, lnTabDim, nnTabDim(0:lnTabDim), m_max, n_max
integer(kind=iwp), intent(out) :: max_nOrd, nMat(nOsc,lBatch), nInc(nOsc,lBatch), nDec(nOsc,lBatch)
real(kind=wp), intent(in) :: dMinWind, C1(nOsc,nOsc), C2(nOsc,nosc), W2(nOsc,nOsc), det0, det1, det2, C(nOsc,nOsc), W(nOsc,nOsc), &
                             r01(nOsc), r02(nOsc), dRho
real(kind=wp), intent(out) :: FC00
integer(kind=iwp) :: ii, m_max_ord, max_nInc, mx_max_ord, n_max_ord, nvTabDim, nx_max_ord
real(kind=wp) :: const, dLT, dRate, dSoc, dSum
real(kind=wp), allocatable :: FCWind2(:)
integer(kind=iwp), parameter :: MB = 1048576

call TabDim(n_max,nOsc,nvTabDim)
max_nOrd = nvTabDim-1
call TabDim(m_max,nOsc,nvTabDim)
m_max_ord = nvTabDim-1
call TabDim(min(n_max,m_max+1),nOsc,nvTabDim)
mx_max_ord = nvTabDim-1
call TabDim(min(m_max,n_max+1),nOsc,nvTabDim)
nx_max_ord = nvTabDim-1
call TabDim(n_max-1,nOsc,nvTabDim)
max_nInc = nvTabDim-1
call TabDim(n_max,nOsc,nvTabDim)
n_max_ord = nvTabDim-1

mx_max_ord = 0 ! CGGn
if (iPrint >= 3) write(u6,*) ' Memory allocated for U matrix:',(n_max_ord+1)*(mx_max_ord+1),' words,  ', &
                             (n_max_ord+1)*(mx_max_ord+1)*ItoB/MB,' MB.'
call XFlush(u6)

!GGt -------------------------------------------------------------------
!write(u6,*) '     lnTabDim+1=',lnTabDim+1,':'
!do i=0,lnTabDim
!  iIndex0 = nnTabDim(i)
!  call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
!  write(u6,*) i,' read at',nnTabDim(i),'  M:',(nMat0(j),j=1,nOsc)
!end do
!write(u6,*) '-----------------------------------------------'
!GGt -------------------------------------------------------------------
!call XFlush(u6)

call mma_allocate(FCWind2,nYes,label='FCWind2')
call ISCD_FCval(iPrint,iMaxYes,lnTabDim,nnTabDim,lNMAT0,lNMAT,lNINC,lNDEC,lBatch,nIndex,C1,det1,r01,C2,W2,det2,r02,m_max_ord, &
                n_max_ord,mx_max_ord,max_nInc,nx_max_ord,nMat,nInc,nDec,C,W,det0,FC00,nOsc,nYes,VibWind2,FCWind2)

const = Two*Pi/hbarcm
if (iPrint >= 4) then
  write(u6,*)
  write(u6,*) '  const =',const
  write(u6,*) '  dRho/cm =',dRho/auTocm
  write(u6,*) '  const*dRho=',const*dRho/auTocm
end if
const = const*dRho/auTocm

dSum = Zero
do ii=1,nYes
  dSum = dSum+FCWind2(ii)**2
end do
call mma_deallocate(FCWind2)

dSoc = Zero
do ii=1,3
  dSOC = dSOC+TranDip(ii)**2
end do

dRate = const*dSum*dSoc/dMinWind
dLT = One/dRate

if (iPrint >= 3) then
  write(u6,*) '  Sum of squares of FC factors =',dSum
  write(u6,*) '  Root-square of the sum =',sqrt(dSum)
  write(u6,*) '  dSOC =',dSOC
end if

if (iPrint >= 1) then
  write(u6,*)
  write(u6,*) ' InterSystem Crossing rate constant:'
  write(u6,*) ' ===================================='
  write(u6,'(a,e10.2,a)') '  ISC Rate Constant  ',dRate,' sec-1'
  write(u6,'(a,e10.2,a)') '  Lifetime           ',dLT,' sec'
  dLT = dLT*1.0e3_wp
  if ((dLT > One) .and. (dLT <= 1.0e3_wp)) write(u6,'(a19,f5.1,a)') ' ',dLT,' msec'
  dLT = dLT*1.0e3_wp
  if ((dLT > One) .and. (dLT <= 1.0e3_wp)) write(u6,'(a19,f5.1,a)') ' ',dLT,' microsec'
  dLT = dLT*1.0e3_wp
  if ((dLT > One) .and. (dLT <= 1.0e3_wp)) write(u6,'(a19,f5.1,a)') ' ',dLT,' nsec'
  dLT = dLT*1.0e3_wp
  if ((dLT > One) .and. (dLT <= 1.0e3_wp)) write(u6,'(a19,f5.1,a)') ' ',dLT,' psec'
  dLT = dLT*1.0e3_wp
  if ((dLT > One) .and. (dLT <= 1.0e3_wp)) write(u6,'(a19,f5.1,a)') ' ',dLT,' fsec'
  write(u6,*) ' ------------------------------------'
  write(u6,*)
  write(u6,*)
  call XFlush(u6)
end if

return

end subroutine ISCD_Rate

subroutine ISCD_Ene(iPrint,nOsc,max_nOrd,nYes,lNMAT,lnTabDim,GE1,GE2,harmfreq1,harmfreq2,x_anharm1,x_anharm2,dMinWind,dRho,nMat0, &
                    nTabDim,lVec)
! Calculate Energy of Levels  GG 30-Dec-08 - 08-Jan-09

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, auTocm
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, nOsc, max_nOrd, lNMAT, lnTabDim, nTabDim(0:lnTabDim)
integer(kind=iwp), intent(inout) :: nYes, lVec(0:lnTabDim)
real(kind=wp), intent(in) :: GE1, GE2, harmfreq1(nOsc), harmfreq2(nOsc), x_anharm1(nOsc,nOsc), x_anharm2(nOsc,nOsc), dRho
real(kind=wp), intent(inout) :: dMinWind
integer(kind=iwp), intent(out) :: nMat0(nOsc)
integer(kind=iwp) :: iIndex, iOrd, j, l_harm, loc_n_max, nYes_start
real(kind=wp) :: dEne, dWlow, dWup
logical(kind=iwp) :: lUpdate
integer(kind=iwp), allocatable :: level1(:), level2(:), lTVec(:)
real(kind=wp), allocatable :: EneMat(:)

if (dMinWind == Zero) then
  lUpDate = .true.
  dMinWind = One
else
  lUpDate = .false.
end if
call mma_allocate(lTVec,[0,max_nOrd],label='lTVec')
lTVec(:) = lVec(0:max_nOrd)

! Energy calculation

if (iPrint >= 4) then
  write(u6,*)
  write(u6,*) ' States in the preliminar window :'
  if (nOsc <= 24) then
    write(u6,'(a,a)') '  ',repeat('=',108)
    write(u6,*) '     jOrd    ene/au    ene/cm-1 Vibrational quantum numbers'
    write(u6,'(a,a)') '  ',repeat('-',108)
  else
    write(u6,'(a,a)') '  ',repeat('=',36)
    write(u6,*) '        #    jOrd   ene/au      ene/cm-1'
    write(u6,'(a,a)') '  ',repeat('-',36)
  end if
  call XFlush(u6)
end if

call mma_allocate(level1,nOsc,label='level1')
call mma_allocate(level2,nOsc,label='level2')
call mma_allocate(EneMat,[0,max_nOrd],label='EneMat')
level1(:) = 0
rewind(lNMAT)
do iOrd=0,max_nOrd
  if (lVec(iOrd) == 1) then
    iIndex = nTabDim(iOrd)
    call iDaFile(lNMAT,2,nMat0,nOsc,iIndex)
    level2(:) = nMat0
    l_harm = nOsc
    call TransEnergy(GE1,x_anharm1,harmfreq1,level1,GE2,x_anharm2,harmfreq2,level2,dEne,l_harm)
    EneMat(iOrd) = dEne
    if (iPrint >= 4) then
      if (nOsc <= 24) then
        loc_n_max = 0
        do j=1,nOsc
          loc_n_max = loc_n_max+nMat0(j)
        end do
        write(u6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,dEne*auTocm,loc_n_max,': ',(nMat0(j),j=1,nOsc)
      else
        write(u6,'(a2,i8,f11.6,f11.4,i4        )') ' ',iOrd,dEne,dEne*auTocm,loc_n_max
      end if
    end if
  end if
end do

! Energy selection

if (iPrint >= 3) then
  write(u6,*)
  write(u6,*) ' States in the window :'
  if (nOsc <= 24) then
    write(u6,'(a,a)') '  ',repeat('=',108)
    write(u6,*) '     jOrd    ene/au    ene/cm-1 Vibrational quantum numbers'
    write(u6,'(a,a)') '  ',repeat('-',108)
  else
    write(u6,'(a,a)') '  ',repeat('=',36)
    write(u6,*) '        #    jOrd   ene/au      ene/cm-1'
    write(u6,'(a,a)') '  ',repeat('-',36)
  end if
  call XFlush(u6)
end if

nYes_start = nYes
do
  dWlow = Half*dMinWind/dRho
  dWup = dWlow
  do iOrd=0,max_nOrd
    lVec(iOrd) = lTVec(iOrd)
    if (lVec(iOrd) == 1) then
      dEne = EneMat(iOrd)
      if ((dEne < -dWlow) .or. (dEne > dWup)) then
        lVec(iOrd) = 0
        nYes = nYes-1
      else
        lVec(iOrd) = 1
        if (iPrint >= 3) then
          if (nOsc <= 24) then
            iIndex = nTabDim(iOrd)
            call iDaFile(lNMAT,2,nMat0,nOsc,iIndex)
            loc_n_max = 0
            do j=1,nOsc
              loc_n_max = loc_n_max+nMat0(j)
            end do
            write(u6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,dEne*auTocm,loc_n_max,': ',(nMat0(j),j=1,nOsc)
          else
            write(u6,'(a2,i8,f11.6,f11.4,i4        )') ' ',iOrd,dEne,dEne*auTocm,loc_n_max
          end if
        end if
      end if
    end if
  end do
  if ((nYes > 1) .or. (.not. lUpDate)) exit
  dMinWind = dMinWind+One
  nYes = nYes_start
end do

call mma_deallocate(lTVec)
call mma_deallocate(level1)
call mma_deallocate(level2)
call mma_deallocate(EneMat)

if (iPrint >= 3) then
  if (nOsc <= 30) write(u6,'(a,a)') '  ',repeat('-',108)
  if (nOsc > 30) write(u6,'(a,a)') '  ',repeat('-',36)
  write(u6,'(a,f12.9,a,f12.9,a)') '  Window: ',-dWlow,' / ',dWup,' (au)'
  write(u6,'(a,f12.6,a,f12.6,a)') '  Window: ',-dWlow*auTocm,' / ',dWup*auTocm,' (cm-1)'
end if
if (iPrint >= 2) then
  write(u6,*) ' Final number of States=',nYes
end if
if ((dMinWind > One) .and. lUpDate .and. (iPrint >= 1)) then
  write(u6,*)
  write(u6,*) ' *** Warning: Expansion factor has been set to ',dMinWind
  write(u6,*)
end if
call XFlush(u6)

return

end subroutine ISCD_Ene

subroutine ISCD_LogEVec(iPrint,nOsc,max_nOrd,minQ,nYes,lNMAT,lnTabDim,nTabDim,nMaxQ,nMat0,lVec)
! Generate Logical Vector of useful States

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, nOsc, max_nOrd, minQ, lNMAT, lnTabDim, nTabDim(0:lnTabDim), nMaxQ(nOsc)
integer(kind=iwp), intent(out) :: nYes, nMat0(nOsc), lVec(0:lnTabDim)
integer(kind=iwp) :: iIndex, iOrd, iOsc, nSumQ

if (iPrint >= 3) then
  write(u6,*) ' Original number of States=',max_nOrd+1
end if
rewind(lNMAT)
iIndex = 0
do iOrd=0,max_nOrd
  iIndex = nTabDim(iOrd)
  call iDaFile(lNMAT,2,nMat0,nOsc,iIndex)
  nSumQ = 0
  lVec(iOrd) = 1
  do iOsc=1,nOsc
    if (nMat0(iOsc) > nMaxQ(iOsc)) lVec(iOrd) = 0
    nSumQ = nSumQ+nMat0(iOsc)
  end do
  if (nSumQ < minQ) lVec(iOrd) = 0
end do

nYes = 0
do iOrd=0,max_nOrd
  if (lVec(iOrd) == 1) nYes = nYes+1
end do

if (iPrint >= 3) then
  write(u6,*) ' Selected number of States=',nYes
end if

return

end subroutine ISCD_LogEVec

subroutine ISCD_MakeGraphs(m_max,maxOrd,Graph1,Graph2,nOsc)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: m_max, nOsc
integer(kind=iwp), intent(out) :: maxOrd, Graph1(m_max+1,nOsc+1), Graph2(m_max+1,m_max+1,nOsc)
integer(kind=iwp) :: i, iOsc, iQ1, iQ2, m, n, nQuanta, nTabDim
integer(kind=iwp), allocatable :: Num(:)

! Initialize.
if (m_max == 0) return
call TabDim(m_max,nOsc,nTabDim)
maxOrd = nTabDim-1

! Set up the vertex table
Graph1(:,:) = 0
Graph1(:,2) = 1
Graph1(1,:) = 1
if (nOsc > 1) then
  do iOsc=2,nOsc
    n = 0
    do nQuanta=0,m_max
      n = n+Graph1(nQuanta+1,iOsc)
      Graph1(nQuanta+1,iOsc+1) = n
    end do
  end do
end if

! set up the arc table
call mma_allocate(Num,[0,m_max],label='Number')
Num(0) = 0
N = 0
do m=1,m_max
  N = N+Graph1(m,nosc+1)
  Num(m) = N
end do
Graph2(:,:,:) = 0
do iOsc=1,nosc
  do iQ1=0,m_max      ! Where we are going
    do iQ2=0,iQ1-1    ! Where we came from
      do i=iQ2+1,iq1  ! Sum over preceding paths
        Graph2(iQ1+1,iQ2+1,iOsc) = Graph1(i+1,iOsc)+Graph2(iQ1+1,iQ2+1,iOsc)
      end do
    end do
  end do
end do

do iQ1=0,m_max  ! Where we are going
  do iQ2=0,iq1  ! Where we came from
    Graph2(iQ1+1,iQ2+1,nOsc) = Graph2(iQ1+1,iQ2+1,nOsc)+Num(iQ1)
  end do
end do

call mma_deallocate(Num)

return

end subroutine ISCD_MakeGraphs

subroutine ISCD_MakenMat(n_max,nOsc,lNMAT,lnTabDim,Graph2,nTabDim,nMat0)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n_max, nOsc, lNMAT, lnTabDim, Graph2(n_max+1,n_max+1,nOsc)
integer(kind=iwp), intent(out) :: nTabDim(0:lnTabDim), nMat0(nOsc)
integer(kind=iwp) :: i, iDet, iDNR, iIndex, iQ, iQuanta, nd, nD_0, nvTabDim
integer(kind=iwp), allocatable :: iVec(:)
integer(kind=iwp), external :: iDetnr

! Initialize.
nMat0(:) = 0
iIndex = 0
nTabDim(0) = iIndex
call iDaFile(lNMAT,1,nMat0,nOsc,iIndex)
call mma_allocate(iVec,nOsc,label='iVec')

! Macrocycle on iQuanta
nD_0 = 0
do iQuanta=1,n_max
  iVec(:) = 0
  iQ = -1
  iVec(1) = -1
  call TabDim(iQuanta,nOsc,nd)
  call TabDim(iQuanta-1,nOsc,nvTabDim)
  nd = nd-nvTabDim

  ! Microcycle on iDet
  do iDet=1,nD
    iVec(1) = iVec(1)+1
    iQ = iQ+1
    if (iQ > iQuanta) then
      do i=1,nOsc-1
        if (iQ <= iQuanta) exit
        iQ = iQ-iVec(i)+1
        iVec(i) = 0
        iVec(i+1) = iVec(i+1)+1
      end do
    end if
    iVec(nOsc) = iQuanta-iq
    iDNR = iDetnr(iVec,Graph2,nosc,n_max)
    iDNR = iDNR-nD_0
    nMat0(:) = iVec
    nTabDim(iDNR+nD_0) = iIndex
    call iDaFile(lNMAT,1,nMat0,nOsc,iIndex)
  end do
  nD_0 = nD_0+nD
end do
call mma_deallocate(iVec)

return

end subroutine ISCD_MakenMat

subroutine ISCD_MakenIncDec(n_max,nOrd,nOsc,lNMAT,lNINC,lNDEC,lBatch,nBatch,leftBatch,nIndex,Graph2,nMat,nInc,nDec)

use mula_global, only: maxMax_n
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: n_max, nOrd, nOsc, lNMAT, lNINC, lNDEC, lBatch, nBatch, leftBatch, Graph2(n_max+1,n_max+1,nOsc)
integer(kind=iwp), intent(inout) :: nIndex(3,0:maxMax_n)
integer(kind=iwp), intent(out) :: nMat(nOsc,lBatch), nInc(nOsc,lBatch), nDec(nOsc,lBatch)
integer(kind=iwp) :: iBatch, ii, iIndex, iOrd, iv, j, jIndex, kIndex
integer(kind=iwp), allocatable :: iVecD(:), iVecI(:)
integer(kind=iwp), external :: iDetnr

!GGt -------------------------------------------------------------------
!write(u6,*)
!write(u6,*) 'CGGt[ISCD_Mk_nIncDec] Infos:'
!write(u6,*) '     nMat(',nOsc,',',lBatch,')'
!write(u6,*) '     n_max,nOrd,nOsc==',n_max,nOrd,nOsc
!write(u6,*) '     lBatch,nBatch,leftBatch==',lBatch,nBatch,leftBatch
!write(u6,*) '----------------------------------------------'
!write(u6,*) '  The nIndex file:'
!do i=1,nBatch+1
!  write(u6,*) i,': ',nIndex(1,i)
!end do
!  write(u6,*) '----------------------------------------------'
!call XFlush(u6)
!GGt -------------------------------------------------------------------

! Initialize

call mma_allocate(iVecI,nOsc,label='iVecI')
call mma_allocate(iVecD,nOsc,label='iVecD')

! Macrocycle iBatch
! Reading nMat

iIndex = 0
jIndex = 0
do iBatch=1,nBatch
  nInc(:,:) = -1
  nDec(:,:) = -1
  kIndex = nIndex(1,iBatch)
  !write(u6,*)'          iBatch=',iBatch,'  kIndex=',kIndex
  call iDaFile(lNMAT,2,nMat,nOsc*lBatch,kIndex)
  !GGt -----------------------------------------------------------------
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,': ',(nMat(k,i+1),k=1,nOsc)
  !end do
  !write(u6,*) '----------------------------------------------'
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
  do ii=1,lBatch

    ! Create nInc.

    iVecI(:) = nMat(:,ii)
    do j=1,nOsc
      iVecI(j) = iVecI(j)+1
      nInc(j,ii) = iDetnr(iVecI,Graph2,nOsc,n_max)
      iVecI(j) = iVecI(j)-1
    end do

    ! Create nDec.

    do j=1,nOsc
      if (nMat(j,ii) /= 0) then
        iVecD(:) = nMat(:,ii)
        iVecD(j) = iVecD(j)-1
        nDec(j,ii) = iDetnr(iVecD,Graph2,nosc,n_max)
        do iv=1,nOsc
          iVecD(iv) = iVecD(j)+1
        end do
      else
        nDec(j,ii) = -1
      end if
    end do

  end do

  nIndex(2,iBatch) = iIndex
  call iDaFile(lNINC,1,nInc,nOsc*lBatch,iIndex)
  nIndex(3,iBatch) = jIndex
  call iDaFile(lNDEC,1,nDec,nOsc*lBatch,jIndex)
  !GGt -----------------------------------------------------------------
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,':I',(nInc(k,i+1),k=1,nOsc)
  !  write(u6,*) i+(iBatch-1)*lBatch,':D',(nDec(k,i+1),k=1,nOsc)
  !end do
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
  !  write(u6,*) '            nInc Written at ',nIndex(2,iBatch)
  !  write(u6,*) '            nDec Written at ',nIndex(3,iBatch)
  !  write(u6,*) '----------------------------------------------'
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------

end do

if (leftBatch > 0) then
  kIndex = nIndex(1,nBatch+1)
  !write(u6,*) '          nBatch+1',nBatch+1,'  kIndex=',kIndex
  call iDaFile(lNMAT,2,nMat,nOsc*lBatch,kIndex)
  !GGt -----------------------------------------------------------------
  !write(u6,*) '  --------- last Batch'
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,': ',(nMat(k,i+1),k=1,nOsc)
  !end do
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------

  !write(u6,*) 'CGGt nBatch*lBatch,nOrd==',nBatch*lBatch,nOrd
  !call XFlush(u6)
  nInc(:,:) = -1
  nDec(:,:) = -1
  do iOrd=nBatch*lBatch,nOrd
    ii = 1+iOrd-nBatch*lBatch

    ! Create nInc.

    iVecI(:) = nMat(:,ii)
    do j=1,nOsc
      iVecI(j) = iVecI(j)+1
      nInc(j,ii) = iDetnr(iVecI,Graph2,nOsc,n_max)
      iVecI(j) = iVecI(j)-1
    end do

    ! Create nDec.

    do j=1,nOsc
      if (nMat(j,ii) /= 0) then
        iVecD(:) = nMat(:,ii)
        iVecD(j) = iVecD(j)-1
        nDec(j,ii) = iDetnr(iVecD,Graph2,nosc,n_max)
        do iv=1,nOsc
          iVecD(iv) = iVecD(j)+1
        end do
      else
        nDec(j,ii) = -1
      end if
    end do

  end do

  nIndex(2,nBatch+1) = iIndex
  call iDaFile(lNINC,1,nInc,nOsc*lBatch,iIndex)
  nIndex(3,nBatch+1) = jIndex
  call iDaFile(lNDEC,1,nDec,nOsc*lBatch,jIndex)
  !GGt -----------------------------------------------------------------
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,':I',(nInc(k,i+1),k=1,nOsc)
  !  write(u6,*) i+(iBatch-1)*lBatch,':D',(nDec(k,i+1),k=1,nOsc)
  !end do
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
  !write(u6,*) '            nInc Written at ',nIndex(2,nBatch+1)
  !write(u6,*) '            nDec Written at ',nIndex(3,nBatch+1)
  !write(u6,*) '----------------------------------------------'
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
end if

!write(u6,*) '----------------------------------------------'
!call XFlush(u6)

call mma_deallocate(iVecI)
call mma_deallocate(iVecD)

return

end subroutine ISCD_MakenIncDec

subroutine ISCD_ReloadNMAT(lnTabDim,nOrd,nOsc,lNMAT0,lNMAT,lBatch,nBatch,leftBatch,nIndex,nTabDim,nMat0,nMat)

use mula_global, only: maxMax_n
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lnTabDim, nOrd, nOsc, lNMAT0, lNMAT, lBatch, nBatch, leftBatch, nTabDim(0:lnTabDim)
integer(kind=iwp), intent(out) :: nIndex(3,0:maxMax_n), nMat0(nOsc), nMat(nOsc,lBatch)
integer(kind=iwp) :: iBatch, ii, iIndex0, iOrd, iOsc, jIndex

! Initialize

!GGt -------------------------------------------------------------------
!write(u6,*)
!write(u6,*) 'CGGt[ISCD_ReloadNMAT] Infos:'
!write(u6,*) '     nMat(',nOsc,',',lBatch,')'
!write(u6,*) '     lnTabDim,nOrd,nOsc==',lnTabDim,nOrd,nOsc
!write(u6,*) '     lBatch,nBatch,leftBatch==',lBatch,nBatch,leftBatch
!write(u6,*) '     lnTabDim+1=',lnTabDim+1,':'
!do i=0,lnTabDim
!  write(u6,*) i,' read at ',nTabDim(i)
!  iIndex0 = nTabDim(i)
!  call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
!  write(u6,*) i,' read at',nTabDim(i),'  M:',(nMat0(j),j=1,nOsc)
!end do
!write(u6,*) '----------------------------------------------'
!call XFlush(u6)
!GGt -------------------------------------------------------------------
jIndex = 0
rewind(lNMAT0)
do iBatch=1,nBatch
  do ii=0,lBatch-1
    iOrd = ii+(iBatch-1)*lBatch
    iIndex0 = nTabDim(iOrd)
    call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
    do iOsc=1,nOsc
      nMat(iOsc,ii+1) = nMat0(iOsc)
    end do
  end do
  nIndex(1,iBatch) = jIndex
  call iDaFile(lNMAT,1,nMat,nOsc*lBatch,jIndex)
  !GGt -----------------------------------------------------------------
  !write(u6,*)'  --------- iBatch =',iBatch
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,':M',(nMat(k,i+1),k=1,nOsc)
  !end do
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
end do
if (leftBatch > 0) then
  do iOrd=nBatch*lBatch,nOrd
    ii = iOrd-nBatch*lBatch
    iIndex0 = nTabDim(iOrd)
    call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
    do iOsc=1,nOsc
      nMat(iOsc,ii+1) = nMat0(iOsc)
    end do
  end do
  nIndex(1,nBatch+1) = jIndex
  call iDaFile(lNMAT,1,nMat,nOsc*lBatch,jIndex)
  !GGt -----------------------------------------------------------------
  !write(u6,*) '  --------- last Batch'
  !do i=0,lBatch-1
  !  write(u6,*) i+(iBatch-1)*lBatch,':M',(nMat(k,i+1),k=1,nOsc)
  !end do
  !call XFlush(u6)
  !GGt -----------------------------------------------------------------
end if
!GGt -------------------------------------------------------------------
!write(u6,*) '----------------------------------------------'
!write(u6,*) '  The nIndex file:'
!do i=1,nBatch+1
! write(u6,*) i,': ',nIndex(1,i)
!end do
!write(u6,*) '----------------------------------------------'
!call XFlush(u6)
!GGt -------------------------------------------------------------------

return

end subroutine ISCD_ReloadNMAT

subroutine ISCD_FCval(iPrint,iMaxYes,lnTabDim,nnTabDim,lNMAT0,lNMAT,lNINC,lNDEC,lBatch,nIndex,C1,det1,r01,C2,W2,det2,r02,max_mOrd, &
                      max_nOrd,max_nOrd2,max_nInc,max_nInc2,nMat,nInc,nDec,C,W,det0,FC00,nOsc,nYes,VibWind2,FCWind2)

use mula_global, only: maxMax_n
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, iMaxYes, lnTabDim, nnTabDim(0:lnTabDim), lNMAT0, lNMAT, lNINC, lNDEC, lBatch, &
                                 nIndex(3,0:maxMax_n), max_mOrd, max_nOrd, max_nOrd2, max_nInc2, nOsc, nYes, VibWind2(nYes)
real(kind=wp), intent(in) :: C1(nOsc,nOsc), det1, r01(nOsc), C2(nOsc,nOsc), W2(nOsc,nOsc), det2, r02(nOsc), C(nOsc,nOsc), &
                             W(nOsc,nOsc), det0
integer(kind=iwp), intent(inout) :: max_nInc
integer(kind=iwp), intent(out) :: nMat(nOsc,lBatch), nInc(nOsc,lBatch), nDec(nOsc,lBatch)
real(kind=wp), intent(out) :: FC00, FCWind2(nYes)
integer(kind=iwp) :: i, iBatch, ii, iIndex, iIndex0, iiOrd, iOrd, j, jIndex, jjOrd, jOrd, kDelta, kIndex, kOsc, kOsc_start, &
                     loc_n_max, lOsc, n, nMaxMat, nTabDim
real(kind=wp) :: const, det, dFC, FC00_exp
integer(kind=iwp), allocatable :: nMat0(:)
real(kind=wp), allocatable :: A2(:,:), A2B2T(:,:), Alpha(:,:), Alpha1(:,:), Alpha2(:,:), B2(:,:), Beta(:,:), d2(:), L(:,:), &
                              r_temp1(:), r_temp2(:), sqr(:), temp(:,:), temp1(:,:), temp2(:,:), U(:,:)
real(kind=wp), external :: Ddot_

call mma_allocate(nMat0,nOsc,label='nMat0')

!GGt -------------------------------------------------------------------
!write(u6,*) 'CGGt[ISCD_FCval] Enter'
!write(u6,*) '     nYes = ',nYes
!write(u6,*) '     VibWind2 :',(VibWind2(i),i=1,nYes)
!write(u6,*) '     L matrix:',max_mOrd,max_nInc2
!write(u6,*) '     U matrix:',max_nOrd,max_nOrd2
!do i=0,iMaxYes
!  write(u6,*) (nInc(i,j),j=1,nOsc)
!end do
!write(u6,*) '     lnTabDim+1=',lnTabDim+1,':'
!do i=0,lnTabDim
!  iIndex0 = nnTabDim(i)
!  !call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
!  write(u6,*) i,' read at',nnTabDim(i),'  M:',(nMat0(j),j=1,nOsc)
!end do
!write(u6,*) '-----------------------------------------------'
!call XFlush(u6)
!GGt -------------------------------------------------------------------

! Initialize.
rewind(lnMAT0)
rewind(lnMAT)
rewind(lnINC)
rewind(lnDEC)
nMaxMat = max(max_mOrd+1,max_nOrd+1)
!write(u6,*) '            nMaxMat=',nMaxMat
nTabDim = max(nMaxMat,8)
!write(u6,*) '            nTabDim=',nTabDim
call mma_allocate(temp,nOsc,nOsc,label='temp')
call mma_allocate(temp1,nOsc,nOsc,label='temp1')
call mma_allocate(temp2,nOsc,nOsc,label='temp2')

! Setup sqr table.
n = nTabDim+1
call mma_allocate(sqr,[0,n],label='sqr')
do i=0,nTabDim+1
  sqr(i) = sqrt(real(i,kind=wp))
end do

! Calculate alpha1, alpha2 and alpha.
!write(u6,*) 'CGGt[FCVal] Calculate alpha(s)'
!call XFlush(u6)
call mma_allocate(Alpha,nOsc,nOsc,label='Alpha')
call mma_allocate(Alpha1,nOsc,nOsc,label='Alpha1')
call mma_allocate(Alpha2,nOsc,nOsc,label='Alpha2')
call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C1,nOsc,C1,nOsc,Zero,Alpha1,nOsc)
call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C2,nOsc,C2,nOsc,Zero,Alpha2,nOsc)
temp(:,:) = alpha1+alpha2
Alpha(:,:) = Half*temp

!call xxDgemul(C,nOsc,'T',C,nOsc,'N',alpha,nOsc,nOsc,nOsc,nOsc)
!alpha(:,:) = Half*alpha

! Calculate C using a Cholesky factorization of 2*alpha.
!call Cholesky(temp,C)

! Calculate W.
!call unitmat(W,nOsc)
!temp(:,:) = C
!call Dool_MULA(temp,W,det0)

! Calculate r00.
call mma_allocate(r_temp1,nOsc,label='r_temp1')
call mma_allocate(r_temp2,nOsc,label='r_temp2')

! Calculate beta.
!write(u6,*) 'CGGt[FCVal] Calculate beta.'
!call XFlush(u6)
call mma_allocate(Beta,nOsc,nOsc,label='Beta')
do i=1,nOsc
  do j=1,nOsc
    temp1(j,i) = C1(i,j)
  end do
  !write(u6,*) 'CGGt C1(',i,',j)=',(C1(i,jj),jj=1,nOsc)
end do
!call XFlush(u6)
temp1(:,:) = Alpha1
temp(:,:) = Two*Alpha

call Dool_MULA(temp,nOsc,nOsc,temp1,nOsc,nOsc,det)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,Alpha2,nOsc,temp1,nOsc,Zero,Beta,nOsc)

call mma_deallocate(Alpha1)
call mma_deallocate(Alpha2)

! Calculate FC00.
!r_temp1(:) = r02-r01
r_temp1(:) = r01-r02

call DGEMM_('N','N',nOsc,1,nOsc,One,Beta,nOsc,r_temp1,nOsc,Zero,r_temp2,nOsc)
FC00_exp = Ddot_(nOsc,r_temp1,1,r_temp2,1)
FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)
!write(u6,*) 'CGGt[FCVal] FC00_exp,FC00=',FC00_exp,FC00
!call XFlush(u6)

call mma_deallocate(Beta)

! Calculate A, B and d matrices.
call mma_allocate(A2,nOsc,nOsc,label='A2')
call mma_allocate(B2,nOsc,nOsc,label='B2')
call mma_allocate(d2,nOsc,label='d2')

call DGEMM_('N','N',nOsc,nOsc,nOsc,One,C2,nOsc,W,nOsc,Zero,A2,nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,One,W2,nOsc,C,nOsc,Zero,temp,nOsc)
B2(:,:) = A2-temp

const = -sqr(8)
call DGEMM_('T','N',nOsc,1,nOsc,const,W2,nOsc,r_temp2,nOsc,Zero,d2,nOsc)

! Calculate A2B2T.
call mma_allocate(A2B2T,nOsc,nOsc,label='A2B2T')

call DGEMM_('N','T',nOsc,nOsc,nOsc,One,A2,nOsc,B2,nOsc,Zero,A2B2T,nOsc)

! Initialize L matrix.
call mma_allocate(L,[0,max_mOrd],[0,max_nInc2],label='L')
L(:,:) = Zero
L(0,0) = One

! If max_mOrd > 0 then set up L(m,0).
if (max_mOrd > 0) then
  write(u6,*) '*****************************************'
  write(u6,*) ' Hot initial state not implemented yet !'
  write(u6,*) '*****************************************'
  write(u6,*)
  call Quit_OnUserError()
end if

! Initialize U matrix.
!write(u6,*) 'CGGt[FCVal] Initialize U matrix.'
!call XFlush(u6)
call mma_allocate(U,[0,max_nOrd],[0,max_nOrd2],label='U')
U(:,:) = Zero
U(0,0) = One
!do kOsc=1,nOsc
!  write(u6,*) 'CGGt[FCVal] d2(..)=',d2(kOsc)
!  call XFlush(u6)
!end do

! Reading for first batch nMat, nInc, nDec

iIndex = nIndex(1,1)
jIndex = nIndex(2,1)
kIndex = nIndex(3,1)
!write(u6,*) 'CGGt[FCVal] Reading nMat at ',iIndex
!call XFlush(u6)
call iDaFile(lNMAT,2,nMat,nOsc*lBatch,iIndex)
!GGt -------------------------------------------------------------------
!write(u6,*) 'CGGt-nMat:'
!do i=0,lBatch-1
!  write(u6,*) i,'-M:',(nMat(j,i+1),j=1,nOsc)
!end do
!write(u6,*) '-----------------------------------------------'
!GGt -------------------------------------------------------------------

!write(u6,*) 'CGGt[FCVal] Reading nInc at ',jIndex
!call XFlush(u6)
call iDaFile(lNINC,2,nInc,nOsc*lBatch,jIndex)
!GGt -------------------------------------------------------------------
!write(u6,*) 'CGGt-nInc:'
!do i=0,lBatch-1
!  write(u6,*) i,'-I:',(nInc(j,i+1),j=1,nOsc)
!end do
!write(u6,*) '-----------------------------------------------'
!GGt -------------------------------------------------------------------

!write(u6,*) 'CGGt[FCVal] Reading nDec at ',kIndex
!call XFlush(u6)
call iDaFile(lNDEC,2,nDec,nOsc*lBatch,kIndex)
!GGt -------------------------------------------------------------------
!write(u6,*) 'CGGt-nDec:'
!do i=0,lBatch-1
!  write(u6,*) i,'-D:',(nDec(j,i+1),j=1,nOsc)
!end do
!write(u6,*) '-----------------------------------------------'
!GGt -------------------------------------------------------------------

! If max_nOrd > 0 then set up U(n,0).

!GGt -------------------------------------------------------------------
!write(u6,*) 'CGGt[FCVal] max_nOrd > 0 then set up U(n,0).'
!call XFlush(u6)
!GGt -------------------------------------------------------------------
do kOsc=1,nOsc
  U(nInc(kOsc,1),0) = d2(kOsc)
end do
max_nInc = min(max_nInc,iMaxYes)
!GGt -------------------------------------------------------------------
!do kOsc=1,nOsc
!  write(u6,*) 'CGGt[FCVal] U(',nInc(kOsc,1),',0) =',U(nInc(kOsc,1),0)
!end do
!write(u6,*) 'CGGt[FCVal] max_nInc=',max_nInc
!write(u6,*) 'CGGt[FCVal]  iMaxYes=',iMaxYes
!call XFlush(u6)
!GGt -------------------------------------------------------------------

iBatch = 1
if (max_nInc > 0) then
  do iOrd=1,max_nInc
    !write(u6,*) '              iOrd=',iOrd,'    iBatch=',iBatch
    !call XFlush(u6)
    iiOrd = iOrd-(iBatch-1)*lBatch+1
    if (iiOrd == (lBatch+1)) then
      iBatch = iBatch+1
      iIndex = nIndex(1,iBatch)
      jIndex = nIndex(2,iBatch)
      kIndex = nIndex(3,iBatch)
      !write(u6,*) 'CGGt[] iBatch=',iBatch,': Reading nMat at ',iIndex
      !call XFlush(u6)
      call iDaFile(lNMAT,2,nMat,nOsc*lBatch,iIndex)
      !write(u6,*) 'CGGt[] iBatch=',iBatch,': Reading nInc at ',iIndex
      !call XFlush(u6)
      call iDaFile(lNINC,2,nInc,nOsc*lBatch,jIndex)
      !write(u6,*) 'CGGt[] iBatch=',iBatch,': Reading nDec at ',iIndex
      !call XFlush(u6)
      call iDaFile(lNDEC,2,nDec,nOsc*lBatch,kIndex)
      iiOrd = iOrd-(iBatch-1)*lBatch+1
    end if
    !write(u6,*) '              iOrd=',iOrd,'  iiOrd=',iiOrd
    !write(u6,*) iOrd,'-MID:',(nMat(j,iiOrd),j=1,nOsc),' /',(nInc(j,iiOrd),j=1,nOsc),' /',(nDec(j,iiOrd),j=1,nOsc)
    !call XFlush(u6)
    kOsc_start = nOsc
    do while ((nMat(kOsc_start,iiOrd) == 0) .and. (kOsc_start > 1))
      kOsc_start = kOsc_start-1
    end do
    do kOsc=kOsc_start,nOsc
      do lOsc=1,nOsc
        !write(u6,*) 'iOrd,kOsc_start,kOsc,lOsc==',iOrd,kOsc_start,kOsc,lOsc
        !call XFlush(u6)
        if (nMat(lOsc,iiOrd) > 0) then
          U(nInc(kOsc,iiOrd),0) = U(nInc(kOsc,iiOrd),0)+sqr(nMat(lOsc,iiOrd))*A2B2T(kOsc,lOsc)*U(nDec(lOsc,iiOrd),0)
          !write(u6,*) iOrd,' U(',nInc(kOsc,iiOrd),')=',U(nInc(kOsc,iiOrd),0)
          !call XFlush(u6)
        end if
      end do
      !write(u6,*) iOrd,' jOrd=',nInc(kOsc,iiOrd)
      jOrd = nInc(kOsc,iiOrd)
      jjOrd = jOrd-(iBatch-1)*lBatch+1
      if ((jjOrd < 1) .or. (jjOrd > lBatch)) then
        iIndex0 = nnTabDim(jOrd)
        call iDaFile(lNMAT0,2,nMat0,nOsc,iIndex0)
        !write(u6,*) jOrd,' read at',nnTabDim(jOrd),'  M:',(nMat0(j),j=1,nOsc),'  jjOrd=',jjOrd
        kDelta = nMat0(kOsc)
      else
        !GGn kDelta = nMat(kOsc,nInc(kOsc,iiOrd)+1)
        kDelta = nMat(kOsc,jjOrd)
      end if
      !Write(u6,*) '         ',iOrd,' >',kDelta
      !Write(u6,*) iOrd,' jOrd=',jOrd  ,'  kDelta=',kDelta
      U(nInc(kOsc,iiOrd),0) = (U(nInc(kOsc,iiOrd),0)+d2(kOsc)*U(iOrd,0))/sqr(kDelta) ! nMat(kOsc,nInc(kOsc,iiOrd)))
    end do
  end do
end if

! Use recursion formula to obtain the rest of U.
!write(u6,*) '            Use recursion ... rest of U.'
!write(u6,*) '            max_nOrd2=',max_nOrd2
!call XFlush(u6)
!!do jOrd=1,max_nOrd2
!!  lOsc = nOsc
!!  do while ((nMat(jOrd,lOsc) == 0) .and. (lOsc > 1))
!!    lOsc = lOsc-1
!!  end do
!!  do iOrd=0,max_nOrd
!!    do kOsc=1,nOsc
!!      if (nMat(iOrd,kOsc) > 0) then
!!        !write(u6,*) '              ',iOrd,kOsc,nMat(iOrd,kOsc)
!!        U(iOrd,jOrd) = U(iOrd,jOrd)+sqr(nMat(iOrd,kOsc))/sqr(nMat(jOrd,lOsc))*A2(kOsc,lOsc)*U(nDec(iOrd,kOsc),nDec(jOrd,lOsc))
!!      end if
!!    end do
!!  end do
!!end do

call mma_deallocate(A2)
call mma_deallocate(B2)
call mma_deallocate(d2)
call mma_deallocate(A2B2T)

! Calculate Franck-Condon factors.
if (iPrint >= 3) then
  write(u6,*) ' Franck-Condon factors for States in the Window:'
  write(u6,'(a,a)') '  ',repeat('=',36)
  write(u6,*) '     #     jOrd   FC factor     jSum'
  write(u6,'(a,a)') '  ',repeat('-',36)
end if
do ii=1,nYes
  jOrd = VibWind2(ii)
  dFC = FC00*L(0,0)*U(jOrd,0)
  FCWind2(ii) = dFC
  if (iPrint >= 3) then
    loc_n_max = 0
    kIndex = nnTabDim(jOrd)
    call iDaFile(lNMAT0,2,nMat0,nOsc,kIndex)
    do j=1,nOsc
      loc_n_max = loc_n_max+nMat0(j)
    end do
    write(u6,'(a2,i5,i9,e15.6,a2,i4)') ' ',ii,jOrd,FCWind2(ii),' ',loc_n_max
  end if
end do
if (iPrint >= 3) then
  write(u6,'(a,a)') '  ',repeat('-',36)
  write(u6,*) ' FC_00 =',FC00
  write(u6,*)
end if

if (iPrint >= 4) then
  write(u6,*)
  write(u6,*) ' Full Franck-Condon factors (FC_00=',FC00,'):'
  write(u6,*) ' =================================================='
  write(u6,*) '    jOrd   FC            level'
  write(u6,*) ' --------------------------------------------------'
  do jOrd=0,max_nInc ! max_nOrd
    loc_n_max = 0
    kIndex = nnTabDim(jOrd)
    call iDaFile(lNMAT0,2,nMat0,nOsc,kIndex)
    do j=1,nOsc
      loc_n_max = loc_n_max+nMat0(j)
    end do
    dFC = FC00*L(0,0)*U(jOrd,0)
    write(u6,'(a,i8,e15.6,a2,i4,a2,24i3)') ' ',jOrd,dFC,' ',loc_n_max,' ',(nMat0(j),j=1,nOsc)
  end do
  write(u6,*) ' --------------------------------------------------'
end if

call mma_deallocate(r_temp1)
call mma_deallocate(r_temp2)
call mma_deallocate(Alpha)
call mma_deallocate(sqr)
call mma_deallocate(temp)
call mma_deallocate(temp1)
call mma_deallocate(temp2)
call mma_deallocate(L)
call mma_deallocate(U)
call mma_deallocate(nMat0)

!write(u6,*) 'CGGt[FCVal] Exit'
!call XFlush(u6)

end subroutine ISCD_FCval
