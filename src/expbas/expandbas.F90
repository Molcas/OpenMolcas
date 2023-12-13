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

subroutine expandbas(Bas1,nBas1,Bas2,nBas2,Orb1,Orb2,occ1,eorb1,indt1,occ2,eorb2,indt2)
! This subroutine expands the MOs for a given symmetry
! Orb1 are the input orbitals of dimension nBas1*bas1 (file INPORB)
! Orb2 are the output orbitals of dimension nBas2*bas2 (file EXPORB)
! Bas1 and Bas2 are the basis set specifications for the old and
! new basis, respectively. They have dimensions nBas1 and nBas2.

#include "intent.fh"

use info_expbas_mod, only: LenIn
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
character(len=LenIn+8), intent(in) :: Bas1(*), Bas2(*)
integer(kind=iwp), intent(in) :: nBas1, nBas2, indt1(*)
integer(kind=iwp), intent(_OUT_) :: indt2(*)
real(kind=wp), intent(in) :: Orb1(*), occ1(*), eorb1(*)
real(kind=wp), intent(_OUT_) :: Orb2(*), occ2(*), eorb2(*)
integer(kind=iwp) :: i, Ibas1, Ibas2, imo, lmo1, lmo2, Nzero
integer(kind=iwp), allocatable :: Izero(:)

! Loop through the new basis labels and compare with the old.
! If they are equal copy orbital coefficients
! If not, add zeros until they fit again

call mma_allocate(Izero,nBas2,label='Izero')

Nzero = 0
Ibas2 = 1
Ibas1 = 1
do while (Ibas1 <= nBas1)
  if (Bas2(Ibas2) == Bas1(Ibas1)) then
    lmo1 = 0
    lmo2 = 0
    do imo=1,nBas1
      Orb2(lmo2+Ibas2) = Orb1(lmo1+Ibas1)
      lmo1 = lmo1+nBas1
      lmo2 = lmo2+nBas2
    end do
    Ibas1 = Ibas1+1
    Ibas2 = Ibas2+1
  else
    Nzero = Nzero+1
    Izero(Nzero) = Ibas2
    lmo2 = 0
    do imo=1,nBas1
      Orb2(lmo2+Ibas2) = Zero
      lmo2 = lmo2+nBas2
    end do
    Ibas2 = Ibas2+1
  end if
end do

! Add zeros at the end of each basis function

if (Ibas2 <= nBas2) then
  do i=Ibas2,nBas2
    Nzero = Nzero+1
    Izero(Nzero) = i
    lmo2 = 0
    do imo=1,nBas1
      Orb2(lmo2+i) = Zero
      lmo2 = lmo2+nBas2
    end do
  end do
end if

! Add new basis functions at the end and expand occ, eorb, and indt

if (nBas1 /= 0) then
  do imo=1,nBas1
    occ2(imo) = occ1(imo)
    eorb2(imo) = eorb1(imo)
    indt2(imo) = indt1(imo)
  end do
end if

if (nBas1 < nBas2) then
  Nzero = 0
  do imo=nBas1+1,nBas2
    Nzero = Nzero+1
    do Ibas2=1,nBas2
      Orb2(nBas2*(imo-1)+Ibas2) = Zero
    end do
    Orb2(nBas2*(imo-1)+Izero(Nzero)) = One
    occ2(imo) = Zero
    eorb2(imo) = Zero
    indt2(imo) = 6
  end do
end if

call mma_deallocate(Izero)

return

end subroutine expandbas
