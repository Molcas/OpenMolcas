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
      Subroutine Reset_ThrGrd(nAtom,nDim,dMass,nSym,iOper,Smmtrc,
     &                 Degen,nIter,Cx,mTtAtm,iAnr,DDV_Schlegel,iOptC,
     &                 rHidden,ThrGrd)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 dMass(nAtom), Degen(3*nAtom), Cx(3*nAtom,nIter)
      Integer iOper(0:nSym-1), iANr(nAtom)
      Logical Smmtrc(3*nAtom),DDV_Schlegel,Found
*                                                                      *
************************************************************************
*                                                                      *
      Call qpg_dArray('Saddle',Found,nSaddle)
      If (.NOT.Found) Return
      iIter = nIter           ! Normal Computation
*                                                                      *
************************************************************************
*                                                                      *
*---- Find the translational and rotational eigenvectors for the
*     current structure.
*
      Call Allocate_Work(ipTR,18*nAtom)
      Call FZero(Work(ipTR),18*nAtom)
*
      Call TRPGen(nDim,nAtom,Cx(1,iIter),Degen,Smmtrc,mTR,dMass,.False.,
     &            Work(ipTR))
*
*     Call RecPrt('Work(ipTR)',' ',Work(ipTR),nDim,mTR)
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('TabAI','Allo','Inte',ip_TabAI,2*mTtAtm)
      Call GetMem('Vect','Allo','Real',ipVec,3*mTtAtm*nDim)
      Call GetMem('AN','Allo','Inte',ipAN,mTtAtm)
      Call GetMem('Coor','Allo','Real',ipCoor,3*mTtAtm)
*
*-----Generate Grand atoms list
*
      Call GenCoo(Cx(1,iIter),nAtom,Work(ipCoor),iOper,nSym,
     &            mTtAtm,Work(ipVec),Smmtrc,nDim,iAnr,iWork(ipAN),
     &            iWork(ip_TabAI),Degen)
*
*---- Are there some hidden frozen atoms ?
*
      nHidden = 0
      nMDstep = 0
      If (rHidden.ge.Two) Call Hidden(mTtAtm,ipCoor,ipAN,nHidden,
     &                                rHidden,nMDstep)
*
*-----Generate bond list
*
      ThrB=0.0D0  ! dummy
      mTtAtm = mTtAtm+nHidden
      Call Box(Work(ipCoor),mTtAtm,iWork(ipAN),iOptC,
     &         ddV_Schlegel,ip_TabB,ip_TabA,nBonds,nMax,ThrB)
      mTtAtm = mTtAtm-nHidden
*                                                                      *
************************************************************************
*                                                                      *
*     If there are some bond types 2, and we are in saddle
*     far from the TS, let us get a reduced threshold to avoid
*     wasting our time
*
      Call Allocate_Work(ipTmp,nSaddle)
      Call Get_dArray('Saddle',Work(ipTmp),nSaddle)
      Found=.false.
      If (Work(ipTmp+nSaddle-2).gt.0.50d0) Then
         Do i=0,nBonds-1
            If (iWork(ip_TabB+3*i+2).eq.2) Then
               Found=.true.
               Go To 20
            EndIf
         End Do
 20      Continue
         If (Found) Then
*           ThrGrd=0.03D0
            ThrGrd=Ten*ThrGrd
            Call WarningMessage(1,
     &             'Molecule composed of many fragments '//
     &             'Convergence threshold reduced')
         EndIf
      EndIf
      Call Free_Work(ipTmp)
*
      Call Free_iWork(ip_TabA)
      Call Free_iWork(ip_TabB)
      Call GetMem('Coor','Free','Real',ipCoor,3*mTtAtm)
      Call GetMem('AN','Free','Inte',ipAN,mTtAtm)
      Call GetMem('Vect','Free','Real',ipVec,3*mTtAtm*nDim)
      Call GetMem('TabAI','Free','Inte',ip_TabAI,2*mTtAtm)
      Call Free_Work(ipTR)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
