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
      SubRoutine PrMtrx(Label,lOper,nComp,ip,Matrix)
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
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
*     Local arrays
      Real*8 Matrix(*)
      Character Label*(*), Line*80
      Integer ip(nComp), lOper(nComp)
      Logical Type
*
      Call qEnter('PrMtrx')
*
      Do 10 iComp = 1, nComp
         ip1 = ip(iComp)
         iSmLbl = lOper(iComp)
         If (Prprt) iSmLbl = iAnd(1,iSmLbl)
         Type = .True.
         Do 30 iIrrep = 0, nIrrep - 1
            If (nBas(iIrrep).le.0) Go To 30
            Do 40 jIrrep = 0, iIrrep
               If (nBas(jIrrep).le.0) Go To 40
               If (iAnd(iSmLbl,2**iEor(iIrrep,jIrrep)).eq.0)
     &            Go To 40
               If (Type) Then
                  Type = .False.
                  Write (6,*)
                  Write (6,*)
                  Write (6,'(A,A,A,I2)')
     &                  ' SO Integrals of type ', Label,' Component ',
     &                     iComp
               End If
               Line=Bline
               If (iIrrep.eq.jIrrep) Then
                  Write (Line,'(1X,A,I1)')
     &            ' Diagonal Symmetry Block ', iIrrep+1
                  Call TriPrt(Line,' ',
     &                 Matrix(ip1),nBas(iIrrep))
                  ip1 = ip1 + nBas(iIrrep)*(nBas(iIrrep)+1)/2
               Else
                  Write (Line,'(1X,A,I1,A,I1)')
     &            ' Off-diagonal Symmetry Block ',
     &            iIrrep+1, ',' , jIrrep+1
                  Call RecPrt(Line,' ',
     &                        Matrix(ip1),nBas(iIrrep),nBas(jIrrep))
                  ip1 = ip1 + nBas(iIrrep)*nBas(jIrrep)
               End If
 40         Continue
 30      Continue
 10   Continue
*
      Call qExit('PrMtrx')
      Return
      End
