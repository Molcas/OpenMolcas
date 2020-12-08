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
     &                 ipBMx,nAtom,nDim,
     &                 Name,Smmtrc,Degen,BSet,HSet,
     &                 nIter,Gx,mTtAtm,iAnr,nStab,jStab,
     &                 Numerical,HWRS,Analytic_Hessian,iOptC,PrQ,
     &                 iCoSet,lOld,iIter,mTR,TRVec,iTabAI,
     &                 iTabAtoms,iTabBonds,nBonds,nMax,
     &                 iRef,nQQ,MaxItr,nWndw)
      use Slapaf_Info, only: Cx, dMass
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      Real*8 Degen(3*nAtom), Gx(3,nAtom,nIter), TRVec(nDim,mTR)
      Character Name(nAtom)*(LENIN)
      Integer   iAnr(nAtom), nStab(nAtom), jStab(0:7,nAtom),
     &          iCoSet(0:7,nAtom)
      Logical Smmtrc(3*nAtom), BSet, HSet, Numerical, HWRS,
     &        Analytic_Hessian, PrQ, lOld
      Integer iTabBonds(3,nBonds), iTabAtoms(0:nMax,nAtom),
     &        iTabAI(2,mTtAtm)

*                                                                      *
************************************************************************
*                                                                      *
*-----Recompute the B matrix once each macroIteration, this is
*     not done if a numerical Hessian is computed.
*
      Call CurviL(nAtom,nDim,Cx,Gx,nIter,iIter,iRef,nStab,
     &            jStab,Degen,Smmtrc,mTR,TRVec,HSet,BSet,ipBMx,
     &            Numerical,iANr,HWRS,Analytic_Hessian,iOptC,
     &            Name,PrQ,dMass,iCoSet,
     &            iTabBonds,iTabAtoms,nBonds,nMax,
     &            iTabAI,mTtAtm,lOld,nQQ,MaxItr,
     &            nWndw)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
