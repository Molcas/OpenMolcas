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

subroutine ccsd_exc(key)
! check, if there is at least one determinant in CCSD expansion
! key=0 - no determinant in expansion
!     1 - only monoexcitations in expansion
!     2 - both mono and biexcitations in expansion

use ccsd_global, only: mmul, noa, nob, nsym, nva, nvb
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: key
integer(kind=iwp) :: asym, bsym, ijsym, isym, jsym, naa, naaaa, nab, nabab, nbb, nbbbb, nij

!1.1 calc # of monoexcitations
!    taking into account also symmetry

naa = 0
nbb = 0
do isym=1,nsym
  asym = isym
  naa = naa+noa(isym)*nva(asym)
  nbb = nbb+nob(isym)*nvb(asym)
end do

!1.2 calc # of biexcitation
!    taking into account also symmetry

naaaa = 0
do isym=1,nsym
  do jsym=1,isym
    ijsym = mmul(isym,jsym)
    if (isym == jsym) then
      nij = noa(isym)*(noa(isym)-1)/2
    else
      nij = noa(isym)*noa(jsym)
    end if
    do asym=1,nsym
      bsym = mmul(ijsym,asym)
      if (bsym < asym) then
        nab = nva(asym)*nva(bsym)
      else if (bsym == asym) then
        nab = nva(asym)*(nva(asym)-1)/2
      else
        nab = 0
      end if
      naaaa = naaaa+nij*nab
    end do
  end do
end do

nbbbb = 0
do isym=1,nsym
  do jsym=1,isym
    ijsym = mmul(isym,jsym)
    if (isym == jsym) then
      nij = nob(isym)*(nob(isym)-1)/2
    else
      nij = nob(isym)*nob(jsym)
    end if
    do asym=1,nsym
      bsym = mmul(ijsym,asym)
      if (bsym < asym) then
        nab = nvb(asym)*nvb(bsym)
      else if (bsym == asym) then
        nab = nvb(asym)*(nvb(asym)-1)/2
      else
        nab = 0
      end if
      nbbbb = nbbbb+nij*nab
    end do
  end do
end do

nabab = 0
do isym=1,nsym
  do jsym=1,isym
    ijsym = mmul(isym,jsym)
    nij = noa(isym)*nob(jsym)
    do asym=1,nsym
      bsym = mmul(ijsym,asym)
      nab = nva(asym)*nvb(bsym)
      nabab = nabab+nij*nab
    end do
  end do
end do

!2 set key

if ((naaaa+nbbbb+nabab) == 0) then
  if ((naa+nbb) == 0) then
    key = 0
  else
    key = 1
  end if
else
  key = 2
end if

return

end subroutine ccsd_exc
