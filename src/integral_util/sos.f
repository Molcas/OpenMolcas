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
      SubRoutine SOS(iStabO,nStabO,lOper)
************************************************************************
*                                                                      *
* Object: to generate the stabilizer S for the operator O.             *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             February '91                                             *
************************************************************************
      use Symmetry_Info, only: nIrrep, iChTbl, iOper
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Integer iStabO(8)
*
#ifdef _DEBUGPRINT_
      Write (6,*) ' In SOS'
      Write (6,*) ' lOper=',lOper
      Do 1 iS = 0, nIrrep-1
         Write(6,'(8I5)') (iChTbl(iIrrep,iS),iIrrep=0,nIrrep-1)
 1    Continue
#endif
      If (lOper.lt.0.or.lOper.gt.255) Then
         Call WarningMessage(2,'SOS: Symmetry label is corrupted.')
         Write (6,*) 'lOper=',lOper
         Call Abend()
      End If
      nStabO = 0
      Do 10 iS = 0, nIrrep-1
         Do 20 iIrrep = 0, nIrrep-1
            If (iAnd(lOper,2**iIrrep).eq.0) Go To 20
            If (iChTbl(iIrrep,iS).ne.1) Go To 10
 20      Continue
         nStabO = nStabO + 1
         iStabO(nStabO) = iOper(iS)
 10   Continue
*
      Return
      End
