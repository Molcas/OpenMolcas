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
*               1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine PrMx2(Label,iComp,lOper,Result,Mem)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              TriPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, Sweden, January '91                  *
*                                                                      *
*     Modified by AB 950620                                            *
************************************************************************
      use Symmetry_Info, only: iOper
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
*     Local arrays
      Real*8 Result(Mem)
      Character Label*(*), Line*80
      Integer  lOper
      Logical Type
*
*     Call qEnter('PrMtrx')
*
         ip1=1
         Type = .True.
         Do 30 iIrrep = 0, nIrrep - 1
            If (nBas(iIrrep).le.0) Go To 30
            Do 40 jIrrep=0,iIrrep
             lop=iEor(iOper(iIrrep),iOper(jIrrep))
               if (lop.ne.loper) Go To 40
               If (nBas(jIrrep).le.0) Go To 40
               If (Type) Then
                  Type = .False.
                  Write (6,*)
                  Write (6,*)
                  Write (6,'(A,A,A,I2)')
     &      ' SO Integral gradients of the ', Label,' Component ',
     &                     iComp
               End If
               Line=''
               If (iIrrep.eq.jIrrep) Then
                  Write (Line,'(1X,A,I1)')
     &            ' Diagonal Symmetry Block ', iIrrep+1
*                  Call TriPrt(Line,' ',Result(ip1),nBas(iIrrep))
                  ip1 = ip1 + nBas(iIrrep)*(nBas(iIrrep)+1)/2
               Else
                  Write (Line,'(1X,A,I1,A,I1)')
     &            ' Off-diagonal Symmetry Block ',
     &            iIrrep+1, ',' , jIrrep+1
                  Call RecPrt(Line,' ',
     &                        Result(ip1),nBas(iIrrep),nBas(jIrrep))
                  ip1 = ip1 + nBas(iIrrep)*nBas(jIrrep)
               End If
 40         Continue
 30      Continue
*
*     Call qExit('PrMtrx')
      Return
      End
