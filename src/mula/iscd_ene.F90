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

subroutine ISCD_LogEVec(iPrint,nOsc,max_nOrd,minQ,nYes,lNMAT,lnTabDim,nTabDim,nMaxQ,nMat0,lVec)
! Generate Logical Vector of useful States

use Definitions, only: u6

implicit real*8(a-h,o-z)
implicit integer(i-n)
dimension nMat0(nOsc)
integer nTabDim(0:lnTabDim), lVec(0:lnTabDim)
integer nMaxQ(nOsc)

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
!####
subroutine ISCD_Ene(iPrint,nOsc,max_nOrd,nYes,lNMAT,lnTabDim,GE1,GE2,harmfreq1,harmfreq2,x_anharm1,x_anharm2,dMinWind,dRho,nMat0, &
                    nTabDim,lVec)
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
integer nMat0(nOsc), nTabDim(0:lnTabDim)
integer lVec(0:lnTabDim)
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
        write(u6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,dEne*HarToRcm,loc_n_max,': ',(nMat0(j),j=1,nOsc)
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
            iIndex = nTabDim(iOrd)
            call iDaFile(lNMAT,2,nMat0,nOsc,iIndex)
            loc_n_max = 0
            do j=1,nOsc
              loc_n_max = loc_n_max+nMat0(j)
            end do
            write(u6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,dEne*HarToRcm,loc_n_max,': ',(nMat0(j),j=1,nOsc)
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

end subroutine ISCD_Ene
