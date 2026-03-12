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
      Subroutine GetCnt(NGROUP,IGROUP,NATOMS,ATLBL)
      use Molcas, only: LenIn
      Implicit None
      Integer NGROUP,IGROUP(8),NATOMS
      Character(LEN=LENIN) ATLBL(*)
! Purpose: Read data from ONEINT

! Read NGROUP
      Call Get_iScalar('nSym',NGROUP)

! Read SYMMETRY GROUP ELEMENTS (Symmetry operations)
      Call Get_iArray('Symmetry operations',IGROUP,NGROUP)

! Read NATOMS, nr of symmetry-unique atoms
      Call Get_iScalar('Unique atoms',NATOMS)

! Read ATLBL, an array of atom labels
      Call Get_cArray('Unique Atom Names',ATLBL,LENIN*NATOMS)

! Read COOR, cartesian coordinates of each atom
!     Call Get_dArray('Unique Coordinates',COOR,3*nAtoms)

      End Subroutine GetCnt
