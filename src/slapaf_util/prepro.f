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
      Subroutine PrePro(nLines,iInt,nFix,nAtom,nInter,Coor)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "info_slapaf.fh"
#include "print.fh"
      Real*8 Coor(3,nAtom)
      Logical CofM
*
      Call QEnter('PrePro')
      iRout=134
      iPrint=nPrint(iRout)
*
      CofM = Iter.eq.1 .and. lNmHss
      Call Allocate_Work(ipTR,18*nAtom)
      Call FZero(Work(ipTR),18*nAtom)
      Call TRPGen(nDimBC,nAtom,Coor,Degen,nSym,iOper,
     &            Smmtrc,mTR,Work(ipCM),CofM,Work(ipTR))
      Call Free_Work(ipTR)
      If (lNmHss) Then
         If (Iter.eq.1) mTROld=mTR
         If (iter.le.2*(nDimBC-mTROld)+1.and.iter.ne.1)
     &    mTR=mTROld
      Else
         mTROld=mTR
      End If
*
*-----Operate according to two modes
*     nLines.gt.0 : user supplied internal coordinates
*     nLines.le.0 : Cartesian or Internal Coordinates
*
      nRowH = 0
      nInter = nDimBC - mTR
      If (nLines.gt.0) Then
*
*--------Find the number of active and frozen internal coordinates.
*
         Call Rd_UDIC(nLines,iInt,nFix,nRowH)
         nQQ=iInt+nFix
         If (nRowH.GT.0) then
            lRowH=.True.
            Call Rd_UDIC_RowH(nQQ,nRowH,mRowH)
         EndIf
         If (nQQ.gt.nInter) Redundant=.True.
*
      Else
*
         nFix=0
      End If
*
*-----Initiate the force constant matrix in internal coordinate
*     basis, excluding rotation and translation.
*
      Call IntFcm(ipH,nQQ,lOld,lOld_Implicit,nAtom,iter)
*
*     Write to runfile only on the first iteration and that there
*     was not an already defined Hessian.
*
      If (iter.eq.1.and.lOld) Then
         Call Put_dArray('Hss_Q',Work(ipH),nQQ**2)
         Call Put_dArray('Hss_upd',Work(ip_Dummy),0)
         Call Free_Work(ipH)
      End If
*
*-----Symmetrize forces
*
      If (LSup) Then
         Call SupSym(Work(ipGrd),nAtom,cMass,Coor,nSupSy,
     &               iWork(ipNSup),iWork(ipAtom),iOper,nSym)
         Call GetMem('iAtom ','Free','Inte',ipAtom,nAtom)
         Call GetMem(' NSUP ','Free','Inte',ipNSup,nSupSy)
      End If
*
      Call QExit('PrePro')
      Return
      End
