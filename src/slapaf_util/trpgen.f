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
      Subroutine TRPGen(nDim,nAtom,Coor,Degen,nSym,iOper,
     &                  Smmtrc,mTR,dMass,CofM,TRVec)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
      Real*8 Coor(3,nAtom), Degen(3*nAtom),
     &       dMass(nAtom), TRVec(3*nAtom*6)
      Integer   iOper(0:nSym-1)
      Logical Smmtrc(3*nAtom), CofM
      Logical g12K
      Data g12K/.True./
      Save g12K
*
      iRout=135
      iPrint = nPrint(iRout)
*
      Call GetMem('TRVec','Allo','Real',ipTR,18*nAtom)
*
*-----Compute the symmetric translations and rotations
*
*     B    (nTR x nDim)
*      tr
*
      Call TRMake(Work(ipTR),Coor,nAtom,nTR,Degen,iOper,nSym,
     &            Smmtrc,nDim,dMass,CofM)
*
      call dcopy_(nTR*nDim,Work(ipTR),1,TRVec,1)
*
      Call GetMem('Scrt','Allo','Real',ipScrt,(3*nAtom)*nTR)
      Call GetMem('GMtrx','Allo','Real',ipG,nTR**2)
      Call GetMem('EVal','Allo','Real',ipEVal,nTR*(nTR+1)/2)
      Call GetMem('EVec','Allo','Real',ipEVec,nTR**2)
*
*-----Eliminate redundancy and produce an orthogonal representation.
*
*        -1/2
*     K g        (nTR x mTR)
*
      Call GetMem('uMtrx','Allo','Real',ipu,nDim)
      call dcopy_(nDim,One,0,Work(ipu),1)
      i=0
      Do iX = 1, 3*nAtom
         If (Smmtrc(iX)) Then
            i = i + 1
            ii = (i-1)*nDim + i
            Call DScal_(nTR,Sqrt(Degen(iX)),
     &                 TRVec((i-1)*nTR+1),1)
         End If
      End Do
*
      Thr_ElRed=1.0D-12
      Call ElRed(TRVec,nTR,nDim,Work(ipG),Work(ipEVal),
     &           Work(ipEVec),mTR,Work(ipU),Work(ipScrt),g12K,
     &           Thr_ElRed)
*
      If (mTR.gt.0) Then
         Call FZero(TRVec,3*nAtom*nTR)
         Call DGEMM_('T','N',
     &               nDim,mTR,nTR,
     &               1.0d0,Work(ipTR),nTR,
     &                     Work(ipEVec),nTR,
     &               0.0d0,TRVec,nDim)
      End If
*
      Call GetMem('uMtrx','Free','Real',ipu,nDim**2)
      Call GetMem('EVec','Free','Real',ipEVec,nTR**2)
      Call GetMem('EVal','Free','Real',ipEVal,nTR*(nTR+1)/2)
      Call GetMem('GMtrx','Free','Real',ipG,nTR**2)
      Call GetMem('TRVec','Free','Real',ipTR,18*nAtom)
      Call GetMem('Scrt','Free','Real',ipScrt,(3*nAtom)*nTR)
*
      Return
      End
