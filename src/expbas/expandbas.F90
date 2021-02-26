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

subroutine expandbas(Bas1,Nbas1,Bas2,Nbas2,Orb1,Orb2,occ1,eorb1,indt1,occ2,eorb2,indt2)

!     This subroutine expands the MOs for a given symmetry
!     Orb1 are the input orbitals of dimension Nbas1*bas1 (file INPORB)
!     Orb2 are the output orbitals of dimension Nbas2*bas2 (file EXPORB)
!     Bas1 and Bas2 are the basis set specifications for the old and
!     new basis, respectively. They have dimensions Nbas1 and Nbas2.

implicit real*8(a-h,o-z)
#include "Molcas.fh"
character*(LENIN8) Bas1(*), Bas2(*)
integer indt1(*), indt2(*)
dimension Orb1(*), Orb2(*), Izero(Nbas2), occ1(*), eorb1(*),occ2(*), eorb2(*)

!     Loop through the new basis labels and compare with the old.
!     If they are equal copy orbital coefficients
!     If not, add zeros until they fit again

Nzero = 0
Ibas2 = 1
if (Nbas1 == 0) go to 200
Ibas1 = 1
100 continue
if (Bas2(Ibas2) == Bas1(Ibas1)) then
  lmo1 = 0
  lmo2 = 0
  do imo=1,Nbas1
    Orb2(lmo2+Ibas2) = Orb1(lmo1+Ibas1)
    lmo1 = lmo1+nBas1
    lmo2 = lmo2+nBas2
  end do
  Ibas1 = Ibas1+1
  Ibas2 = Ibas2+1
else if (Bas2(Ibas2) /= Bas1(Ibas1)) then
  Nzero = Nzero+1
  Izero(Nzero) = Ibas2
  lmo2 = 0
  do imo=1,Nbas1
    Orb2(lmo2+Ibas2) = 0.d0
    lmo2 = lmo2+nBas2
  end do
  Ibas2 = Ibas2+1
end if
if (Ibas1 <= Nbas1) go to 100
200 continue

! Add zeros at the end of each basis function

if (Ibas2 <= Nbas2) then
  do i=Ibas2,Nbas2
    Nzero = Nzero+1
    Izero(Nzero) = i
    lmo2 = 0
    do imo=1,Nbas1
      Orb2(lmo2+i) = 0.d0
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

if (Nbas1 < Nbas2) then
  Nzero = 0
  do imo=Nbas1+1,Nbas2
    Nzero = Nzero+1
    do Ibas2=1,Nbas2
      Orb2(Nbas2*(imo-1)+Ibas2) = 0.d0
    end do
    Orb2(Nbas2*(imo-1)+Izero(Nzero)) = 1.d0
    occ2(imo) = 0.d0
    eorb2(imo) = 0.d0
    indt2(imo) = 6
  end do
end if

return

end subroutine expandbas
