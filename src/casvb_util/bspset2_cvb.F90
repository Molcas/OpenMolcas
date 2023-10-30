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

subroutine bspset2_cvb(ikcoff,nel,kbasis,need)

use casvb_global, only: i2s_fr, mnion_fr, mxion_fr, nalf_fr, nel_fr, nfrag, nMs_fr, nS_fr
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nel, kbasis
integer(kind=iwp), intent(inout) :: ikcoff(0:nel,0:nel,0:nel)
integer(kind=iwp), intent(out) :: need
integer(kind=iwp) :: i2s1, ifrag, iMs, ion, iS, nalf1, nalf1_spin, nalfsing, ndet, nel1, nelsing
integer(kind=iwp), external :: ifns_cvb

do ifrag=1,nfrag
  do ion=mnion_fr(ifrag),mxion_fr(ifrag)
    nelsing = nel_fr(ifrag)-2*ion
    if (nelsing < 0) cycle
    do iMs=1,nMs_fr(ifrag)
      nalfsing = nalf_fr(iMs,ifrag)-ion
      if (nalfsing < 0) cycle
      do iS=1,nS_fr(ifrag)
        if ((i2s_fr(iS,ifrag) <= nelsing) .and. (i2s_fr(iS,ifrag) >= 2*nalfsing-nelsing)) &
          ikcoff(nelsing,nalfsing,i2s_fr(iS,ifrag)) = 1
      end do
    end do
  end do
end do
need = 0
do nel1=0,nel
  do nalf1=0,nel
    do i2s1=0,nel
      if (ikcoff(nel1,nalf1,i2s1) == 1) then
        ikcoff(nel1,nalf1,i2s1) = need
        nalf1_spin = (nel1+i2s1)/2
        call icomb_cvb(nel1,nalf1,ndet)
        need = need+ifns_cvb(nel1,nalf1_spin,kbasis)*ndet
      end if
    end do
  end do
end do

return

end subroutine bspset2_cvb
