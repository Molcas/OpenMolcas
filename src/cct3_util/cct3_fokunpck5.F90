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

subroutine cct3_fokunpck5(symp,foka,fokb,dpa,dpb,dimfok,rc)
! this routine produces dpa,dpb from foka,fokb
! for some cases
! shifto,shiftv will be also added
!
! symp   - symmtry of this block
! foka   - Fok aa matrix (I)
! fokb   - Fok bb matrix (I)
! dpa    - Diagonal part alfa vector (O)
! dpa    - Diagonal part beta vector (O)
! dimfok - dimension for Fok matrix - norb (I)
! rc     - return (error) code

use CCT3_global, only: eps, keysa, noa, nob, norb, shifto, shiftv, typden
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: symp, dimfok, rc
real(kind=wp) :: foka(dimfok,dimfok), fokb(dimfok,dimfok), dpa(dimfok), dpb(dimfok)
integer(kind=iwp) :: nhelp1, nhelp2, p

rc = 0

if (typden == 0) then
  !1 diagonal elements are required

  do p=1,dimfok
    dpa(p) = foka(p,p)
    dpb(p) = fokb(p,p)
  end do

else if (typden == 1) then
  !2 (faa+fbb)/2 are required

  do p=1,dimfok
    dpa(p) = (foka(p,p)+fokb(p,p))/2
    dpb(p) = dpa(p)
  end do

else if (typden == 2) then
  !3 orbital energies are required

  !3.1 def shift
  if (symp == 1) then
    nhelp1 = 0
  else
    nhelp1 = 0
    do nhelp2=1,symp-1
      nhelp1 = nhelp1+norb(nhelp2)
    end do
  end if

  !3.2 map oe to dp
  do p=1,dimfok
    dpa(p) = eps(nhelp1+p)
    dpb(p) = eps(nhelp1+p)
  end do

else
  ! RC=1 : invalid key (NCI/Stup)
  rc = 1
end if

if ((keysa == 3) .or. (keysa == 4)) then
  ! for full adaptation scheme only D and V orbitals are shifted

  do p=1,nob(symp)
    dpa(p) = dpa(p)-shifto
    dpb(p) = dpb(p)-shifto
  end do

  do p=1+noa(symp),norb(symp)
    dpa(p) = dpa(p)+shiftv
    dpb(p) = dpb(p)+shiftv
  end do

else
  ! for other schemes all orbitals are shifted

  do p=1,noa(symp)
    dpa(p) = dpa(p)-shifto
  end do

  do p=1,nob(symp)
    dpb(p) = dpb(p)-shifto
  end do

  do p=1+noa(symp),norb(symp)
    dpa(p) = dpa(p)+shiftv
  end do

  do p=1+nob(symp),norb(symp)
    dpb(p) = dpb(p)+shiftv
  end do

end if

return

end subroutine cct3_fokunpck5
