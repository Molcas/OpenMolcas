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
!***********************************************************************

module rasdef

! NRSPRT = Nr of RAS partitions
! NRAS(ISYM,IP)=Nr of orbitals with symmetry ISYM in each part.
! NRASEL(IP)=Min nr of accumulated electrons

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: MXRSPRT = 3
integer(kind=iwp) :: NRAS(8,MXRSPRT), NRASEL(MXRSPRT), NRS1(8), NRS1T, NRS2(8), NRS2T, NRS3(8), NRS3T, NRSPRT

public :: NRAS, NRASEL, NRS1, NRS1T, NRS2, NRS2T, NRS3, NRS3T, NRSPRT

end module rasdef
