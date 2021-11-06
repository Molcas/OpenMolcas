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

subroutine ISCD_Rate(iPrint,nOsc,max_nOrd,iMx_nOrd,iMaxYes,nYes,dMinWind,lBatch,nBatch,leftBatch,nIndex,VibWind2,lNMAT0,lNMAT, &
                     lNINC,lNDEC,lnTabDim,nnTabDim,C1,C2,W1,W2,det0,det1,det2,C,W,r01,r02,r00,m_max,n_max,max_dip,nnsiz,FC00, &
                     FCWind2,dRho,mTabDim,mMat,mInc,mDec,nMat,nInc,nDec)

implicit real*8(a-h,o-z)
implicit integer(i-n)
#include "Constants_mula.fh"
#include "inout.fh"
#include "WrkSpc.fh"
#include "io_mula.fh"
integer nIndex(3,0:maxMax_n)
real*8 C1(nOsc,nOsc), C2(nOsc,nosc), W1(nOsc,nOsc), W2(nOsc,nOsc), C(nOsc,nOsc), W(nOsc,nOsc)
real*8 r01(nOsc), r02(nOsc), r00(nOsc), det0, det1, det2, FC00
real*8 FCWind2(nYes)
integer VibWind2(nYes), nnTabDim(0:lnTabDim)
integer mMat(0:mTabDim,nOsc), nMat(nOsc,lBatch)
integer mInc(0:mTabDim,nOsc), nInc(nOsc,lBatch)
integer mDec(0:mTabDim,nOsc), nDec(nOsc,lBatch)

call TabDim2_drv(m_max,nosc,nvTabDim)
call TabDim2_drv(n_max,nosc,nvTabDim)
max_nOrd = nvTabDim-1
call TabDim2_drv(m_max,nosc,nvTabDim)
m_max_ord = nvTabDim-1
call TabDim2_drv(min(n_max,m_max+1),nosc,nvTabDim)
mx_max_ord = nvTabDim-1
call TabDim2_drv(min(m_max,n_max+1),nosc,nvTabDim)
nx_max_ord = nvTabDim-1
call TabDim2_drv(m_max-1,nosc,nvTabDim)
max_mInc = nvTabDim-1
call TabDim2_drv(n_max-1,nosc,nvTabDim)
max_nInc = nvTabDim-1
call TabDim2_drv(n_max,nosc,nvTabDim)
n_max_ord = nvTabDim-1

mx_max_ord = 0 ! CGGn
if (iPrint >= 3) write(6,*) ' Memory allocated for U matrix:',(n_max_ord+1)*(mx_max_ord+1),' words,  ', &
                            8*(n_max_ord+1)*(mx_max_ord+1)/1024/1024,' MB.     '
call XFlush(6)
call GetMem('L','Allo','Real',ipL,(m_max_ord+1)*(nx_max_ord+1))
call GetMem('U','Allo','Real',ipU,(n_max_ord+1)*(mx_max_ord+1))
call GetMem('alpha1','Allo','Real',ipAlpha1,nOsc*nOsc)
call GetMem('alpha2','Allo','Real',ipAlpha2,nOsc*nOsc)
call GetMem('beta','Allo','Real',ipBeta,nOsc*nOsc)
call GetMem('MAT0','Allo','Inte',ipnMat0,nOsc)

!GGt -------------------------------------------------------------------
!write(6,*) '     lnTabDim+1=',lnTabDim+1,':'
!do i=0,lnTabDim
!  iIndex0 = nnTabDim(i)
!  call iDaFile(lNMAT0,2,iWork(ipnMat0),nOsc,iIndex0)
!  write(6,*) i,' read at',nnTabDim(i),'  M:',(iWork(ipnMat0+j),j=0,nOsc-1)
!end do
!write(6,*) '-----------------------------------------------'
!GGt -------------------------------------------------------------------
!call GetMem('Test_2','LIST','INTE',iDum,iDum)
!call XFlush(6)
call ISCD_FCval(iPrint,iMaxYes,lnTabDim,nnTabDim,lNMAT0,lNMAT,lNINC,lNDEC,lBatch,nBatch,leftBatch,nIndex,C1,W1,det1,r01,C2,W2, &
                det2,r02,m_max_ord,n_max_ord,mx_max_ord,max_mInc,max_nInc,nx_max_ord,mMat,nMat,mInc,nInc,mDec,nDec,C,W,det0,r00, &
                Work(ipL),Work(ipU),FC00,Work(ipAlpha1),Work(ipAlpha2),Work(ipBeta),nOsc,nnsiz,iMx_nOrd,nYes,VibWind2,FCWind2, &
                iWork(ipnMat0))

const = 2.0d0*rpi/5.309d-12
if (iPrint >= 4) then
  write(6,*)
  write(6,*) '  const =',const
  write(6,*) '  dRho/cm =',dRho/HarToRcm
  write(6,*) '  const*dRho=',const*dRho/HarToRcm
end if
const = const*dRho/HarToRcm

dSum = 0.0d0
do ii=1,nYes
  dSum = dSum+FCWind2(ii)**2
end do

dSoc = 0.0d0
do ii=1,3
  dSOC = dSOC+TranDip(ii)**2
end do

dRate = const*dSum*dSoc/dMinWind
dLT = 1.0d0/dRate

if (iPrint >= 3) then
  write(6,*) '  Sum of squares of FC factors =',dSum
  write(6,*) '  Root-square of the sum =',sqrt(dSum)
  write(6,*) '  dSOC =',dSOC
end if

if (iPrint >= 1) then
  write(6,*)
  write(6,*) ' InterSystem Crossing rate constant: '
  write(6,*) ' ===================================='
  write(6,'(a,e10.2,a)') '  ISC Rate Constant  ',dRate,' sec-1'
  write(6,'(a,e10.2,a)') '  Lifetime           ',dLT,' sec'
  dLT = dLT*1.0d3
  if ((dLT > 1.0d0) .and. (dLT <= 1.0d3)) write(6,'(a19,f5.1,a)') ' ',dLT,' msec'
  dLT = dLT*1.0d3
  if ((dLT > 1.0d0) .and. (dLT <= 1.0d3)) write(6,'(a19,f5.1,a)') ' ',dLT,' microsec'
  dLT = dLT*1.0d3
  if ((dLT > 1.0d0) .and. (dLT <= 1.0d3)) write(6,'(a19,f5.1,a)') ' ',dLT,' nsec'
  dLT = dLT*1.0d3
  if ((dLT > 1.0d0) .and. (dLT <= 1.0d3)) write(6,'(a19,f5.1,a)') ' ',dLT,' psec'
  dLT = dLT*1.0d3
  if ((dLT > 1.0d0) .and. (dLT <= 1.0d3)) write(6,'(a19,f5.1,a)') ' ',dLT,' fsec'
  write(6,*) ' ------------------------------------'
  write(6,*)
  write(6,*)
  call XFlush(6)
end if

call GetMem('MAT0','Free','Inte',ipnMat0,nOsc)
call GetMem('beta','Free','Real',ipBeta,nOsc*nOsc)
call GetMem('alpha2','Free','Real',ipAlpha2,nOsc*nOsc)
call GetMem('alpha1','Free','Real',ipAlpha1,nOsc*nOsc)
call GetMem('U','Free','Free',ipU,(n_max_ord+1)*(mx_max_ord+1))
call GetMem('L','Free','Real',ipL,(m_max_ord+1)*(nx_max_ord+1))

return

! Avoid unused argument warnings
if (.false.) call Unused_integer(max_dip)

end subroutine ISCD_Rate
