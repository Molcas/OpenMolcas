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
! Copyright (C) 2008,2009, Giovanni Ghigo                              *
!***********************************************************************

!  InterSystem Crossing rate evaluation: "Engine" routines
!  Author: Giovanni Ghigo
!          Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
!          28 Dec-08 - 06 Jan-09

subroutine ISC_Rho(iPrint,nOsc,new_n_max,dRho,energy1,energy2,minQ,dMinWind,nMaxQ,harmfreq1,harmfreq2)
! Calculate State Density  dRho  GG 30-Dec-08
! Formula (86) taken from  M. Bixon, J. Jortner  JCP,48,715 (1969)

use Constants, only: Zero, One, Two, Six, Twelve, Half
use Definitions, only: wp, u6

implicit real*8(a-h,o-z)
implicit integer(i-n)
#include "Constants_mula.fh"
real*8 energy1, energy2, harmfreq1(nOsc), harmfreq2(nOsc)
real*8 GE1, GE2, dMinWind, dMinWind0
integer nMaxQ(nOsc)

if (iPrint >= 2) then
  write(u6,*)
  write(u6,*) ' State Density data:                         '
  write(u6,*) ' ============================================'
end if

dMinWind0 = dMinWind
if (dMinWind == Zero) dMinWind = One
minQ = 0
dRho = Two/rpi
dRho = sqrt(dRho*nOsc)
dRho = dRho*(One-One/(Twelve*nOsc))
avFreq = Zero
avFreqSq = Zero
dLambda = One
dZPE1 = Zero
dZPE2 = Zero
dMinFreq2 = 1.0e99_wp
dMaxFreq2 = Zero
do iOsc=1,nOsc
  avFreq = avFreq+harmfreq2(iOsc)
  avFreqSq = avFreqSq+(harmfreq2(iOsc))**2
  dZPE1 = dZPE1+Half*harmfreq1(iOsc)
  dZPE2 = dZPE2+Half*harmfreq2(iOsc)
  dMinFreq2 = min(dMinFreq2,harmfreq2(iOsc))
  dMaxFreq2 = max(dMaxFreq2,harmfreq2(iOsc))
end do
GE1 = energy1+dZPE1
GE2 = energy2+dZPE2
T0 = abs(GE2-GE1)
avFreq = avFreq/nOsc
avFreqSq = avFreqSq/nOsc

do jOsc=1,nOsc
  dLambda = dLambda*harmfreq2(jOsc)/avFreq
end do
new_n_max = int(Half+(GE1-GE2)/dMinFreq2)
dAlpha = avFreqSq/(avFreq**2)
dBeta = ((nOsc-1)*(nOsc-2)*dAlpha-nOsc**2)/(Six*nOsc)
dEtha = abs(GE1-GE2)/dZPE2
dRho = dRho/avFreq
dRho = dRho*dAlpha
dRho = dRho/(1+dEtha)
dDn = (One+(Two/dEtha))
dDn = dDn**(dEtha/Two)
dDn = dDn*(One+(dEtha/Two))
dDn = dDn**nOsc
dRho = dRho*dDn
dFE = (One+dEtha)**2
dFE = (One-One/dFE)**dBeta
dRho = dRho*dFE

do jOsc=1,nOsc
  nMaxQ(jOsc) = int(Half+(T0+dMinWind/dRho)/harmfreq2(jOsc))
end do
minQ = int(Half+(T0-dMinWind/dRho)/dMaxFreq2)
if (minQ < 0) then
  write(u6,*)
  write(u6,*) ' ***** ERROR ******'
  write(u6,*) ' Window too large !'
  write(u6,*) ' ******************'
  call Quit_OnUserError()
end if

if (iPrint >= 2) then
  write(u6,'(a,f11.6,a)') '  T_0  = ',T0,' (au)'
  write(u6,'(a,f11.3,a)') '  T_0  = ',T0*HarToRcm,' (cm-1)'
  write(u6,'(a,f11.3,a)') '  T_0  = ',T0*auToeV,' (eV)'
  write(u6,'(a,d14.3,a)') '  State Density (dRho) = ',dRho,' (au-1)'
  write(u6,'(a,g14.3,a)') '  State Density (dRho) = ',dRho/HarToRcm,' (cm)'
  write(u6,'(a,g17.9,a)') '  1/dRho = ',HarToRcm/dRho,' (cm-1)'
  write(u6,'(a,f7.3,a)') '  Expansion factor =',dMinWind
  write(u6,'(a,g17.9,a)') '  Window = (+/-)',Half*dMinWind*HarToRcm/dRho,' (cm-1)'
end if
if (iPrint >= 3) then
  write(u6,*) ' Maximum quantum numbers:',(nMaxQ(i),i=1,nOsc)
  write(u6,*) ' Minimum quantum number: ',minQ
  write(u6,*) ' Suggested n_max (new_n_max)=',new_n_max
  write(u6,*)
end if
call XFlush(u6)
dMinWind = dMinWind0

return

end subroutine ISC_Rho
!####
subroutine ISC_Ene(iPrint,nOsc,max_nOrd,nYes,nMat,nTabDim,GE1,GE2,harmfreq1,harmfreq2,x_anharm1,x_anharm2,dMinWind,dRho,lVec)
! Calculate Energy of Levels  GG 30-Dec-08 - 08-Jan-09

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: u6

implicit real*8(a-h,o-z)
implicit integer(i-n)
#include "Constants_mula.fh"
real*8 GE1, GE2, harmfreq1(nOsc), harmfreq2(nOsc)
real*8 x_anharm1(nOsc,nOsc), x_anharm2(nOsc,nOsc)
real*8 dMinWind, dRho, dWlow, dWup
real*8 dEne
integer nMat(0:nTabDim,nOsc), lVec(0:nTabDim)
logical lUpdate
integer, allocatable :: level1(:), level2(:), lTVec(:)
real*8, allocatable :: EneMat(:)

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
    write(u6,'(a,108a)') '  ',('=',i=1,108)
    write(u6,*) '     jOrd    ene/au    ene/cm-1 Vibrational quantum numbers'
    write(u6,'(a,108a)') '  ',('-',i=1,108)
  else
    write(u6,'(a,36a)') '  ',('=',i=1,36)
    write(u6,*) '        #    jOrd   ene/au      ene/cm-1 '
    write(u6,'(a,36a)') '  ',('-',i=1,36)
  end if
  call XFlush(u6)
end if

call mma_allocate(level1,nOsc,label='level1')
call mma_allocate(level2,nOsc,label='level2')
call mma_allocate(EneMat,[0,max_nOrd],label='EneMat')
level1(:) = 0
do iOrd=0,max_nOrd
  if (lVec(iOrd) == 1) then
    level2(:) = nMat(iOrd,:)
    l_harm = nOsc
    call TransEnergy(GE1,x_anharm1,harmfreq1,level1,GE2,x_anharm2,harmfreq2,level2,dEne,l_harm)
    EneMat(iOrd) = dEne
    if (iPrint >= 4) then
      if (nOsc <= 24) then
        loc_n_max = 0
        do j=1,nOsc
          loc_n_max = loc_n_max+nMat(iOrd,j)
        end do
        write(u6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,dEne*HarToRcm,loc_n_max,': ',(nMat(iOrd,j),j=1,nOsc)
      else
        write(u6,'(a2,i8,f11.6,f11.4,i4        )') ' ',iOrd,dEne,dEne*HarToRcm,loc_n_max
      end if
    end if
  end if
end do

! Energy selection

if (iPrint >= 3) then
  write(u6,*)
  write(u6,*) ' States in the window :'
  if (nOsc <= 24) then
    write(u6,'(a,108a)') '  ',('=',i=1,108)
    write(u6,*) '     jOrd    ene/au    ene/cm-1 Vibrational quantum numbers'
    write(u6,'(a,108a)') '  ',('-',i=1,108)
  else
    write(u6,'(a,36a)') '  ',('=',i=1,36)
    write(u6,*) '        #    jOrd   ene/au      ene/cm-1 '
    write(u6,'(a,36a)') '  ',('-',i=1,36)
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
            loc_n_max = 0
            do j=1,nOsc
              loc_n_max = loc_n_max+nMat(iOrd,j)
            end do
            write(u6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,dEne*HarToRcm,loc_n_max,': ',(nMat(iOrd,j),j=1,nOsc)
          else
            write(u6,'(a2,i8,f11.6,f11.4,i4        )') ' ',iOrd,dEne,dEne*HarToRcm,loc_n_max
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
  if (nOsc <= 30) write(u6,'(a,108a)') '  ',('-',i=1,108)
  if (nOsc > 30) write(u6,'(a,36a)') '  ',('-',i=1,36)
  write(u6,'(a,f12.9,a,f12.9,a)') '  Window: ',-dWlow,' / ',dWup,' (au)'
  write(u6,'(a,f12.6,a,f12.6,a)') '  Window: ',-dWlow*HarToRcm,' / ',dWup*HarToRcm,' (cm-1)'
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

end subroutine ISC_Ene
!####
subroutine ISC_Rate(iPrint,nOsc,max_nOrd,iMx_nOrd,iMaxYes,nYes,dMinWind,VibWind2,C1,C2,W1,W2,det0,det1,det2,C,W,r01,r02,r00, &
                    mTabDim,mMat,nTabDim,nMat,mInc,mDec,nInc,nDec,m_max,n_max,max_dip,nnsiz,FC00,dRho)
! Estimate ISC rate  GG 30-Dec-08

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, u6

implicit real*8(a-h,o-z)
implicit integer(i-n)
#include "Constants_mula.fh"
#include "inout.fh"
real*8 C1(nOsc,nOsc), C2(nOsc,nosc), W1(nOsc,nOsc), W2(nOsc,nOsc), C(nOsc,nOsc), W(nOsc,nOsc)
real*8 r01(nOsc), r02(nOsc), r00(nOsc), det0, det1, det2, FC00
integer VibWind2(nYes)
integer mMat(0:mTabDim,nOsc), nMat(0:nTabDim,nOsc)
integer mInc(0:mTabDim,nOsc), nInc(0:iMaxYes,nOsc)
integer mDec(0:mTabDim,nOsc), nDec(0:iMaxYes,nOsc)
real*8, allocatable :: FCWind2(:)

call TabDim(m_max,nosc,nvTabDim)
call TabDim(n_max,nosc,nvTabDim)
max_nOrd = nvTabDim-1
call TabDim(m_max,nosc,nvTabDim)
m_max_ord = nvTabDim-1
call TabDim(min(n_max,m_max+1),nosc,nvTabDim)
mx_max_ord = nvTabDim-1
call TabDim(min(m_max,n_max+1),nosc,nvTabDim)
nx_max_ord = nvTabDim-1
call TabDim(m_max-1,nosc,nvTabDim)
max_mInc = nvTabDim-1
call TabDim(n_max-1,nosc,nvTabDim)
max_nInc = nvTabDim-1
call TabDim(n_max,nosc,nvTabDim)
n_max_ord = nvTabDim-1

mx_max_ord = 0 ! CGGn
if (iPrint >= 3) write(u6,*) ' Memory allocated for U matrix:',(n_max_ord+1)*(mx_max_ord+1),' words,  ', &
                             8*(n_max_ord+1)*(mx_max_ord+1)/1048576,' MB.     '
call XFlush(u6)

call mma_allocate(FCWind2,nYes,label='FCWind2')
call ISC_FCval(iPrint,iMaxYes,nTabDim,C1,W1,det1,r01,C2,W2,det2,r02,m_max_ord,n_max_ord,mx_max_ord,max_mInc,max_nInc,nx_max_ord, &
               mMat,nMat,mInc,nInc,mDec,nDec,C,W,det0,r00,FC00,nOsc,nnsiz,iMx_nOrd,nYes,VibWind2,FCWind2)

! where does this number come from?
const = Two*rpi/5.309e-12_wp
if (iPrint >= 4) then
  write(u6,*)
  write(u6,*) '  const =',const
  write(u6,*) '  dRho/cm =',dRho/HarToRcm
  write(u6,*) '  const*dRho=',const*dRho/HarToRcm
end if
const = const*dRho/HarToRcm

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
  write(u6,*) ' InterSystem Crossing rate constant: '
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
! Avoid unused argument warnings
if (.false.) call Unused_integer(max_dip)

end subroutine ISC_Rate
