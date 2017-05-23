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
      Subroutine BMtrx_Internal(
     &                 nLines,ipBMx,nAtom,nInter,ip_rInt,Coor,nDim,
     &                 rMass,Name,nSym,iOper,Smmtrc,Degen,BSet,HSet,
     &                 nIter,ip_drInt,Gx,Cx,mTtAtm,iAnr,nStab,jStab,
     &                 Numerical,HWRS,Analytic_Hessian,iOptC,PrQ,mxdc,
     &                 iCoSet,lOld,nFix,iIter,mTR,TRVec,ip_TabAI,
     &                 ip_TabA,ip_TabB,nBonds,nMax,
     &                 iRef,ip_KtB,nQQ,Redundant,nqInt,MaxItr,nWndw)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Coor(3,nAtom), rMass(nAtom), Degen(3*nAtom),
     &       Gx(3*nAtom,nIter), Cx(3*nAtom,nIter), TRVec(nDim,mTR)
      Character Name(nAtom)*(LENIN)
      Integer   iOper(0:nSym-1), iAnr(nAtom),
     &          nStab(nAtom), jStab(0:7,nAtom), iCoSet(0:7,nAtom)
      Logical Smmtrc(3*nAtom), BSet, HSet, Redundant,
     &        Numerical, HWRS, Analytic_Hessian, PrQ, lOld
*                                                                      *
************************************************************************
*                                                                      *
      iRout=133
      iPrint=nPrint(iRout)
*                                                                      *
************************************************************************
*                                                                      *
*-----Recompute the B matrix once each macroIteration, this is
*     not done if a numerical Hessian is computed.
*
      Call GetMem('Proj','Allo','Real',ipProj,nDim)
*
      Call CurviL(nAtom,nDim,Cx,Gx,nIter,iIter,iRef,nStab,iOper,
     &            nSym,jStab,Degen,Smmtrc,mTR,TRVec,
     &            ip_rInt,ip_drInt,HSet,BSet,ipBMx,
     &            Numerical,iANr,HWRS,Analytic_Hessian,iOptC,
     &            Name,PrQ,Work(ipProj),rMass,iCoSet,
     &            iWork(ip_TabB),iWork(ip_TabA),nBonds,nMax,
     &            iWork(ip_TabAI),mTtAtm,lOld,ip_KtB,nQQ,nqInt,MaxItr,
     &            nWndw)
*
      Call GetMem('Proj','Free','Real',ipProj,nDim)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nLines)
         Call Unused_integer(nInter)
         Call Unused_real_array(Coor)
         Call Unused_integer(mxdc)
         Call Unused_integer(nFix)
         Call Unused_real_array(TRVec)
         Call Unused_logical(Redundant)
      End If
      End
