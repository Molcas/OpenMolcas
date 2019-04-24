************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine GetCnt(NGROUP,IGROUP,NATOMS,ATLBL)
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "WrkSpc.fh"
      Dimension IGROUP(8)
      Character*(LENIN) ATLBL(*)
C Purpose: Read data from ONEINT

C Read NGROUP
      Call Get_iScalar('nSym',NGROUP)

C Read SYMMETRY GROUP ELEMENTS (Symmetry operations)
      Call Get_iArray('Symmetry operations',IGROUP,NGROUP)

C Read NATOMS, nr of symmetry-unique atoms
      Call Get_iScalar('Unique atoms',NATOMS)

C Read ATLBL, an array of atom labels
      Call Get_cArray('Unique Atom Names',ATLBL,LENIN*NATOMS)

C Read COOR, cartesian coordinates of each atom
!     Call Get_dArray('Unique Coordinates',COOR,3*nAtoms)

      Return
      End
