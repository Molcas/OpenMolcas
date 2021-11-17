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
