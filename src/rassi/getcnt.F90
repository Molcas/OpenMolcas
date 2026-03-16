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

subroutine GetCnt(NGROUP,IGROUP,NATOMS,ATLBL)
! Purpose: Read data from ONEINT

use Molcas, only: LenIn

implicit none
integer NGROUP, IGROUP(8), NATOMS
character(len=LenIn) ATLBL(*)

! Read NGROUP
call Get_iScalar('nSym',NGROUP)

! Read SYMMETRY GROUP ELEMENTS (Symmetry operations)
call Get_iArray('Symmetry operations',IGROUP,NGROUP)

! Read NATOMS, nr of symmetry-unique atoms
call Get_iScalar('Unique atoms',NATOMS)

! Read ATLBL, an array of atom labels
call Get_cArray('Unique Atom Names',ATLBL,LenIn*NATOMS)

! Read COOR, cartesian coordinates of each atom
!call Get_dArray('Unique Coordinates',COOR,3*nAtoms)

end subroutine GetCnt
