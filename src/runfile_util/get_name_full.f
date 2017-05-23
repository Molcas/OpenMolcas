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
      Subroutine Get_Name_Full(Element)
      Implicit None
#include "Molcas.fh"
#include "WrkSpc.fh"
      Character*2 Element(*)
      Integer nAtom, nAtMM, ipLabMM, i
      Logical Found
*
      Call Get_nAtoms_All(nAtom)
      Call Get_Name_All(Element)
*
      Call Qpg_cArray('MMO Labels',Found,nAtMM)
      If (Found) Then
        nAtMM=nAtMM/LENIN
        Call GetMem('MMO Labels','ALLO','CHAR',ipLabMM,LENIN*nAtMM)
        Call Get_cArray('MMO Labels',cWork(ipLabMM),LENIN*nAtMM)
        Do i=1,nAtMM
          Element(nAtom+i)=cWork(ipLabMM+(i-1)*LENIN)//
     &                     cWork(ipLabMM+(i-1)*LENIN+1)
          If (cWork(ipLabMM+(i-1)*LENIN+1).eq.'_')
     &      Element(nAtom+i)(2:2)=' '
        End Do
        Call GetMem('MMO Labels','FREE','CHAR',ipLabMM,LENIN*nAtMM)
      End If
*
      Return
      End
