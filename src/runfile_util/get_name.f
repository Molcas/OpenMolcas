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
      Subroutine Get_Name(Element)
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
#include "periodic_table.fh"
      Character*2 Element(*)
      Real*8, Allocatable :: Chrg(:)
*
      Call Get_iScalar('Unique atoms',nAtoms)
      Call mma_Allocate(Chrg,nAtoms)
*
      Call Get_dArray('Nuclear charge',Chrg,nAtoms)
*
      Do i = 1, nAtoms
         iElement_Nr=INT(Chrg(i))
         If (iElement_Nr.ge.0 .and. iElement_Nr.le.Num_Elem) Then
            Element(i)=PTab(iElement_Nr)
         Else
            Element(i)=' X'
         End If
      End Do
      call mma_deallocate(Chrg)
*
      Return
      End
