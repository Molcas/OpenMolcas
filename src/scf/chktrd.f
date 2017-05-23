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
* Copyright (C) 1996, Martin Schuetz                                   *
************************************************************************
      SubRoutine ChkTrD(nSym,nBas,nOrb,Occ,Dlt,Ovl)
************************************************************************
*                                                                      *
*     purpose: Compute trace of density matrix and compare with sum    *
*              over occupation numbers...                              *
*                                                                      *
*     input:                                                           *
*       nSym    : number of symmetries                                 *
*       nBas(i) : number of basis functions (i = 1, nSym)              *
*       Occ     : occupation numbers                                   *
*       Dlt     : density matrix in triangular storrage                *
*       Ovl     : overlap matrix                                       *
*                                                                      *
*     called from: ??????                                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*
      Implicit Real*8 (a-h,o-z)
*
*     declaration of subroutine parameters...
      Real*8 Occ(*),Dlt(*),Ovl(*)
      Integer nBas(nSym),nOrb(nSym)
*
*     declaration of some local variables...
      Real*8 ThrDif
      Data ThrDif /1.0d-7/
*
#include "real.fh"
*
#ifdef _DEBUG_
      Call QEnter('ChkTrD')
#endif
      ipDlt = 1
      ipOvl = 1
      ipOcc = 0
      Scale= One
      Do iSym = 1, nSym
        nBs = nBas(iSym)
        nOr = nOrb(iSym)
        lth = nBs*(nBs + 1)/2
*       count occupation number...
        SumOcc=Zero
        Do iOr = 1, nOr
          SumOcc=SumOcc+Occ(ipOcc+iOr)*Scale
        End Do
*       do trace of PS for symmetry block...
        TrDns=DDOT_(lth,Dlt(ipDlt),1,Ovl(ipOvl),1)
        ipDlt = ipDlt+lth
        ipOvl = ipOvl+lth
        ipOcc = ipOcc+nOr
        If (Abs(SumOcc-TrDns).gt.ThrDif) Then
*         print Warning...
          Call WarningMessage(1,
     &   'WARNING: trace of density is inconsistent with occupation !')
          Write(6,'(A,I1,A,3F12.7)') 'SymBlock: ',iSym,' deviation: ',
     &     SumOcc-TrDns,SumOcc,TrDns
        End If
      End Do
#ifdef _DEBUG_
      Call QExit('ChkTrD')
#endif
      Return
      End
