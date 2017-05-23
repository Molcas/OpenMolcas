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
      Subroutine MxLbls(GrdMax,StpMax,GrdLbl,StpLbl,nInter,Grad,Shift,
     &                  Lbl)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Shift(nInter), Grad(nInter)
      Character Lbl(nInter)*8,  GrdLbl*8, StpLbl*8
*
*     Call RecPrt('MxLbls:Shift',' ',Shift,nInter,1)
*     Call RecPrt('MxLbls:Grad',' ',Grad,nInter,1)
*
      GrdMax=Zero
      StpMax=Zero
      Do i = 1, nInter
         If(Abs(Grad(i)).gt.Abs(GrdMax)) Then
           GrdMax=Grad(i)
           GrdLbl=Lbl(i)
         End If
         If(Abs(Shift(i)).gt.Abs(StpMax)) Then
           StpMax=Shift(i)
           StpLbl=Lbl(i)
         End If
      End Do
*     Write (*,*) ' Tmp output in MxLbls'
*     Write (*,*) GrdLbl,' ',GrdMax
*     Write (*,*) StpLbl,' ',StpMax
*
      Return
      End
