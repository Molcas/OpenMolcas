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
      Integer Function MltLbl(Lbl1,Lbl2)
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
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      Integer Lbl1, Lbl2
*
      Integer iIrrep, jIrrep, ijSym
*
      MltLbl = 0
      Do iIrrep = 0, nIrrep - 1
         If (iAnd(Lbl1,2**iIrrep).eq.0) Cycle
         Do jIrrep = 0, nIrrep - 1
            If (iAnd(Lbl2,2**jIrrep).eq.0) Cycle
            ijSym = iEor(iIrrep,jIrrep)
            If (iAnd(MltLbl,2**ijSym).eq.0) MltLbl = MltLbl + 2**ijSym
         End Do
      End Do
*
      Return
      End
