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

subroutine fokunpck5(symp,foka,fokb,dpa,dpb,dimfok,rc)
! this routine produces dpa,dpb from foka,fokb
! for some cases
! shifto,shiftv will be also added
!
! symp   - symmetry of this block
! foka   - Fok aa matrix (I)
! fokb   - Fok bb matrix (I)
! dpa    - Diagonal part alpha vector (O)
! dpb    - Diagonal part beta vector (O)
! dimfok - dimension for Fok matrix - norb (I)
! rc     - return (error) code

use ccsd_global, only: eps, fullprint, keysa, noa, nob, norb, shifto, shiftv, typden
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: symp, dimfok
real(kind=wp), intent(in) :: foka(dimfok,dimfok), fokb(dimfok,dimfok)
real(kind=wp), intent(inout) :: dpa(dimfok), dpb(dimfok)
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp) :: nhelp1, p

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
  end do
  dpb(1:dimfok) = dpa(1:dimfok)

else if (typden == 2) then
  !3 orbital energies are required

  !3.1 def shift
  if (symp == 1) then
    nhelp1 = 0
  else
    nhelp1 = sum(norb(1:symp-1))
  end if

  !3.2 map oe to dp
  dpa(1:dimfok) = eps(nhelp1+1:nhelp1+dimfok)
  dpb(1:dimfok) = eps(nhelp1+1:nhelp1+dimfok)

else
  ! RC=1 : invalid key (NCI/Stup)
  rc = 1
end if

if ((keysa == 3) .or. (keysa == 4)) then
  ! for full adaptation scheme only D and V orbitals are shifted

  dpa(1:nob(symp)) = dpa(1:nob(symp))-shifto
  dpb(1:nob(symp)) = dpb(1:nob(symp))-shifto

  dpa(1:noa(symp):norb(symp)) = dpa(1+noa(symp):norb(symp))+shiftv
  dpb(1:noa(symp):norb(symp)) = dpb(1+noa(symp):norb(symp))+shiftv

else
  ! for other schemes all orbitals are shifted

  dpa(1:noa(symp)) = dpa(1:noa(symp))-shifto

  dpb(1:nob(symp)) = dpb(1:nob(symp))-shifto

  dpa(1+noa(symp):norb(symp)) = dpa(1+noa(symp):norb(symp))+shiftv

  dpb(1+nob(symp):norb(symp)) = dpb(1+nob(symp):norb(symp))+shiftv

end if

if (fullprint >= 2) then
  write(u6,*) ' Diagonal part Dp aa, bb for irrep: ',symp
  do p=1,norb(symp)
    write(u6,99) p,dpa(p),dpb(p)
  end do
end if

return

99 format(2x,i4,2(f20.14,2x))

end subroutine fokunpck5
