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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

!***********************************************************************
!*                                                                     *
!*  CNFINI   := Set NVBR, NDETVB, NDETVB2, MNION, MXION,               *
!*              NCONFION, and IFSC.                                    *
!*                                                                     *
!***********************************************************************
subroutine cnfini_cvb(iconfs,nconf1,nel1,nS,i2s,nMs,nalf1,nvbr1,ndetvb1,ndetvb21,mnion1,mxion1,nconfion,ifsc1)

use casvb_global, only: noe, norb
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nconf1, iconfs(noe,nconf1), nel1, nS, i2s(nS), nMs, nalf1(nMs)
integer(kind=iwp), intent(out) :: nvbr1, ndetvb1, ndetvb21, mnion1, mxion1, ifsc1
integer(kind=iwp), intent(_OUT_) :: nconfion(0:*)
integer(kind=iwp) :: i, iconf, iMs, ion, iorb, iretval, iretval1, iretval2, iS

! Main loop over configurations:
mnion1 = nel1/2
mxion1 = 0
nconfion(0:nel1/2) = 0
ndetvb1 = 0
ndetvb21 = 0
nvbr1 = 0
do iconf=1,nconf1
  ion = 0
  do iorb=1,norb
    if (iconfs(iorb,iconf) == 2) ion = ion+1
  end do
  if (ion < mnion1) mnion1 = ion
  if (ion > mxion1) mxion1 = ion
  nconfion(ion) = nconfion(ion)+1
  do iS=1,nS
    call icomb_cvb(nel1-2*ion,(nel1-i2s(iS))/2-ion,iretval1)
    call icomb_cvb(nel1-2*ion,(nel1-i2s(iS))/2-ion-1,iretval2)
    nvbr1 = nvbr1+iretval1-iretval2
  end do
  do iMs=1,nMs
    call icomb_cvb(nel1-2*ion,nalf1(iMs)-ion,iretval)
    ndetvb1 = ndetvb1+iretval
    ndetvb21 = ndetvb21+(iretval+1)/2
  end do
end do
if ((norb == nel1) .and. (nconf1 == 1)) then
  ifsc1 = 1
  do i=1,nel1
    if (iconfs(i,nconf1) /= 1) then
      ifsc1 = 0
      exit
    end if
  end do
else
  ifsc1 = 0
end if

return

end subroutine cnfini_cvb
