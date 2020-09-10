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
      Subroutine Get_LblCnt_All(xLblCnt)
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
#include "Molcas.fh"
      Real*8, Allocatable:: Coord(:,:)
      Character*(LENIN) xLblCnt(*), xLblCnt_Unique(MxAtom)
*
      Call Get_iScalar('Unique atoms',nAtoms)
      Call mma_allocate(Coord,3,nAtoms,Label='Coord')
      Call Get_dArray('Unique Coordinates',Coord,3*nAtoms)
      Call Get_Name(xLblCnt_Unique)
      Call Get_cArray('Unique Atom Names',xLblCnt_Unique,LENIN*nAtoms)
      Call Get_Name_All_(Coord,nAtoms,nAtoms_all,xLblCnt_Unique,xLblCnt)
      Call mma_deallocate(Coord)
*
      Return
      End
