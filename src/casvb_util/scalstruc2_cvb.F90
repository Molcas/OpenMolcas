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

subroutine scalstruc2_cvb(orbs,cvb,iconfs,ifnss)

use casvb_global, only: i2s_fr, nconf, nconfion_fr, nel, nfrag, noe, norb, nS_fr, nvb, sc
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: orbs(norb,norb)
real(kind=wp), intent(inout) :: cvb(nvb)
integer(kind=iwp), intent(in) :: iconfs(noe,nconf), ifnss(0:nel,0:nel)
integer(kind=iwp) :: i, iconf_off, ifrag, ion, iorb, iS, istr, nelsing, nss
real(kind=wp) :: fac, fac1, fac2
real(kind=wp), external :: ddot_

if (sc) then
  fac = One
  do iorb=1,norb
    fac2 = ddot_(norb,orbs(:,iorb),1,orbs(:,iorb),1)
    fac = fac*sqrt(fac2)
  end do
  cvb(:) = fac*cvb(:)
else
  do iorb=1,norb
    fac2 = ddot_(norb,orbs(:,iorb),1,orbs(:,iorb),1)
    fac1 = sqrt(fac2)
    istr = 0
    iconf_off = 0
    do ifrag=1,nfrag
      do iS=1,nS_fr(ifrag)
        do ion=0,nel/2
          nelsing = nel-2*ion
          nss = ifnss(nelsing,i2s_fr(iS,ifrag))
          do i=iconf_off+1,iconf_off+nconfion_fr(ion,ifrag)
            if (iconfs(iorb,i) == 1) then
              cvb(istr+1:istr+nss) = fac1*cvb(istr+1:istr+nss)
            else if (iconfs(iorb,i) == 2) then
              cvb(istr+1:istr+nss) = fac2*cvb(istr+1:istr+nss)
            end if
            istr = istr+nss
          end do
          iconf_off = iconf_off+nconfion_fr(ion,ifrag)
        end do
      end do
    end do
    if (istr /= nvb) then
      write(u6,*) ' ISTR not equal to NVB in SCALSTRUC! ',istr,nvb
      call abend_cvb()
    end if
  end do
end if

return

end subroutine scalstruc2_cvb
