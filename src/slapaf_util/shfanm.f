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
      Subroutine ShfANM(nInter,nIter,rInt,Shift,iPrint)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 rInt(nInter,nIter), Shift(nInter,nIter)
*
      If (nIter.eq.1) Return
      If (iPrint.ge.19) Call RecPrt(' ShfANM: rInt',' ',
     &                              rInt,nInter,nIter)
*
*-----Shifts: dq = q   -q
*               n   n+1  n
*
      Do Iter=1,nIter-1
         call dcopy_(nInter,rInt(1,Iter+1),1,Shift(1,Iter),1)
         Call DaXpY_(nInter,-One,rInt(1,Iter),1,Shift(1,Iter),1)
      End Do
      If (iPrint.ge.19) Call RecPrt(' In ShfANM: New Shifts',' ',
     &                              Shift,nInter,nIter-1)
*
      Return
      End
