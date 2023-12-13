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
! symp   - symmetry of this block
! foka   - Fok aa matrix (I)
! fokb   - Fok bb matrix (I)
! dpa    - Diagonal part alfa vector (O)
! dpa    - Diagonal part beta vector (O)
! dimfok - dimension for Fok matrix - norb (I)
! rc     - return (error) code

use CCT3_global, only: eps, keysa, noa, nob, norb, shifto, shiftv, typden
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: symp, dimfok
real(kind=wp), intent(in) :: foka(dimfok,dimfok), fokb(dimfok,dimfok)
real(kind=wp), intent(out) :: dpa(dimfok), dpb(dimfok)
integer(kind=iwp), intent(out) :: rc
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
    dpa(p) = (foka(p,p)+fokb(p,p))*Half
  end do
  dpb(:) = dpa

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
  dpa(:) = eps(nhelp1+1:nhelp1+dimfok)
  dpb(:) = eps(nhelp1+1:nhelp1+dimfok)

else
  ! RC=1 : invalid key (NCI/Stup)
  rc = 1
end if

if ((keysa == 3) .or. (keysa == 4)) then
  ! for full adaptation scheme only D and V orbitals are shifted

  dpa(1:nob(symp)) = dpa(1:nob(symp))-shifto
  dpa(noa(symp)+1:norb(symp)) = dpa(noa(symp)+1:norb(symp))+shiftv

  dpb(1:nob(symp)) = dpb(1:nob(symp))-shifto
  dpb(noa(symp)+1:norb(symp)) = dpb(noa(symp)+1:norb(symp))+shiftv

else
  ! for other schemes all orbitals are shifted

  dpa(1:noa(symp)) = dpa(1:noa(symp))-shifto
  dpa(noa(symp)+1:norb(symp)) = dpa(noa(symp)+1:norb(symp))+shiftv

  dpb(1:nob(symp)) = dpb(1:nob(symp))-shifto
  dpb(nob(symp)+1:norb(symp)) = dpb(nob(symp)+1:norb(symp))+shiftv

end if

return

end subroutine cct3_fokunpck5
