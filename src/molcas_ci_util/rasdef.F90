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
integer, parameter :: MXRSPRT = 3
integer :: NRS1(8), NRS2(8), NRS3(8), NRS1T, NRS2T, NRS3T, NRSPRT, NRAS(8,MXRSPRT), NRASEL(MXRSPRT)

public :: NRS1, NRS2, NRS3, NRS1T, NRS2T, NRS3T, NRSPRT, NRAS, NRASEL

end module rasdef
