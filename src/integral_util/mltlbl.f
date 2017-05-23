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
      Integer Function MltLbl(Lbl1,Lbl2,nIrrep)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             February '91                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
*
      iRout =213
      iPrint = nPrint(iRout)
*     Call qEnter('MltLbl')
*
      MltLbl = 0
      Do 10 iIrrep = 0, nIrrep - 1
         If (iAnd(Lbl1,2**iIrrep).eq.0) Go To 10
         Do 20 jIrrep = 0, nIrrep - 1
            If (iAnd(Lbl2,2**jIrrep).eq.0) Go To 20
            ijSym = iEor(iIrrep,jIrrep)
            If (iAnd(MltLbl,2**ijSym).eq.0) MltLbl = MltLbl + 2**ijSym
 20      Continue
 10   Continue
*
*     Call qExit('MltLbl')
      Return
      End
