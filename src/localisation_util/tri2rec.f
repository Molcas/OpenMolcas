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
      Subroutine Tri2Rec(OvlTri,OvlRec,nBas,Debug)
c
      Implicit Real*8 (a-h,o-z)
      Real*8 OvlTri(*),OvlRec(nBas,nBas)
      Logical Debug
c
      ioffset=1
      Do nlig=1,nBas
        nElem=nlig
        ioffset=ioffset+(nlig-1)
        call dcopy_(nElem,OvlTri(ioffset),1,OvlRec(1,nlig),1)
      End Do
c
      Do iBas=1,nBas
        Do jBas=nBas,iBas,-1
          OvlRec(jBas,iBas)=OvlRec(iBas,jBas)
        End Do
      End Do
c
      If (Debug) Call RecPrt('OvlRec ',' ',OvlRec,nBas,nBas)
c
      Return
      End
