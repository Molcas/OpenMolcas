************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine ModInd(IndxZt,IndZt,nZeta,CutZt,IncZet,SkipZt,A,B)
************************************************************************
* Object : to modify the index array.                                  *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN.                              *
*             June '91                                                 *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 A(3), B(3)
      Integer IndZt(nZeta), IndxZt(nZeta)
      Logical CutZt, SkipZt
*
      iRout = 246
      iPrint = nPrint(iRout)
*
      AB2 = (A(1)-B(1))**2 + (A(2)-B(2))**2 + (A(3)-B(3))**2
      CutZt = .False.
      SkipZt = .False.
*
*     If (AB2.gt.Zero) Then
         If (IncZet.eq.nZeta) Then
            Call ICopy(nZeta,IndxZt,1,IndZt,1)
            Ind = Abs(IndZt(nZeta))
         Else
            Ind = 0
            Jnd = 0
            Do 10 iZeta = 1, nZeta
               If (Mod(iZeta-1,IncZet).eq.0) Jnd = 0
               IndZt(iZeta) = -Jnd
*              Check if primitive function can be omitted
               If (IndxZt(iZeta).gt.0) Then
                  Ind = Ind + 1
                  Jnd = Jnd + 1
                  IndZt(iZeta) = Jnd
               End If
 10         Continue
         End If
*        Set flag to indicate if arguments were ignored
         CutZt = nZeta .ne. Ind
*        Set flag to indicate that all arguments were ignored.
         SkipZt= Ind.eq.0
*     End If
*
      Return
      End
