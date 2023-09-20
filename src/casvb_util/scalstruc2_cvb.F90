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

use casvb_global, only: i2s_fr, nconfion_fr, nfrag, nS_fr

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
dimension orbs(norb,norb), cvb(nvb)
dimension iconfs(noe,nconf), ifnss(0:nel,0:nel)

if (sc) then
  fac = one
  do iorb=1,norb
    fac2 = ddot_(norb,orbs(1,iorb),1,orbs(1,iorb),1)
    fac = fac*sqrt(fac2)
  end do
  call dscal_(nvb,fac,cvb,1)
else
  do iorb=1,norb
    fac2 = ddot_(norb,orbs(1,iorb),1,orbs(1,iorb),1)
    fac1 = sqrt(fac2)
    istr = 0
    iconf_off = 0
    do ifrag=1,nfrag
      do iS=1,nS_fr(ifrag)
        do ion=0,nel/2
          nelsing = nel-2*ion
          do i=iconf_off+1,iconf_off+nconfion_fr(ion,ifrag)
            if (iconfs(iorb,i) == 1) then
              call dscal_(ifnss(nelsing,i2s_fr(iS,ifrag)),fac1,cvb(istr+1),1)
            else if (iconfs(iorb,i) == 2) then
              call dscal_(ifnss(nelsing,i2s_fr(iS,ifrag)),fac2,cvb(istr+1),1)
            end if
            istr = istr+ifnss(nelsing,i2s_fr(iS,ifrag))
          end do
          iconf_off = iconf_off+nconfion_fr(ion,ifrag)
        end do
      end do
    end do
    if (istr /= nvb) then
      write(6,*) ' ISTR not equal to NVB in SCALSTRUC! ',istr,nvb
      call abend_cvb()
    end if
  end do
end if

return

end subroutine scalstruc2_cvb
