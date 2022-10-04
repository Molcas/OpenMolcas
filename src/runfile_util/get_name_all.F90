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

subroutine Get_Name_All(Element)

use Definitions, only: iwp

implicit none
character(len=2) :: Element(*)
#include "WrkSpc.fh"
#include "Molcas.fh"
integer(kind=iwp) :: ipCoord, nAtoms, nAtoms_all
character(len=2) :: Element_Unique(MxAtom)

call Get_iScalar('Unique atoms',nAtoms)
call Allocate_Work(ipCoord,3*nAtoms)
call Get_dArray('Unique Coordinates',Work(ipCoord),3*nAtoms)
call Get_Name(Element_Unique)
call Get_Name_All_(Work(ipCoord),nAtoms,nAtoms_all,Element_Unique,Element)
call Free_Work(ipCoord)

return

end subroutine Get_Name_All
