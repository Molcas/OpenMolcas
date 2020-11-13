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
      Subroutine Reset_ThrGrd(nAtom,nDim,dMass,Smmtrc,
     &                 Degen,nIter,Cx,mTtAtm,iAnr,DDV_Schlegel,iOptC,
     &                 rHidden,ThrGrd)
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Real*8 dMass(nAtom), Degen(3*nAtom), Cx(3*nAtom,nIter)
      Integer iANr(nAtom)
      Logical Smmtrc(3*nAtom),DDV_Schlegel,Found
      Integer, Allocatable:: TabAI(:), AN(:)
      Real*8, Allocatable:: TR(:), Vec(:), Coor(:,:), Tmp(:)
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
      Call mma_allocate(TR,18*nAtom,Label='TR')
      TR(:)=Zero
*
      Call TRPGen(nDim,nAtom,Cx(1,iIter),Degen,Smmtrc,mTR,dMass,.False.,
     &            TR)
*
*     Call RecPrt('TR',' ',TR,nDim,mTR)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_Allocate(TabAI,2*mTtAtm,Label='TabAI')
      Call mma_Allocate(Vec,3*mTtAtm*nDim,Label='Vec')
      Call mma_Allocate(AN,mTtAtm,Label='AN')
      Call mma_Allocate(Coor,3,mTtAtm,Label='Coor')
*
*-----Generate Grand atoms list
*
      Call GenCoo(Cx(1,iIter),nAtom,Coor,mTtAtm,Vec,Smmtrc,nDim,iAnr,
     &            AN,TabAI,Degen)
*
*---- Are there some hidden frozen atoms ?
*
      nHidden = 0
      nMDstep = 0
      ipCoor = ip_of_Work(Coor(1,1))
      ipAN   = ip_of_iWork(AN(1))
      If (rHidden.ge.Two) Call Hidden(mTtAtm,ipCoor,ipAN,nHidden,
     &                                rHidden,nMDstep)
*
*-----Generate bond list
*
      ThrB=0.0D0  ! dummy
      mTtAtm = mTtAtm+nHidden
      Call Box(Coor,mTtAtm,AN,iOptC,ddV_Schlegel,ip_TabB,ip_TabA,nBonds,
     &         nMax,ThrB)
      mTtAtm = mTtAtm-nHidden
*                                                                      *
************************************************************************
*                                                                      *
*     If there are some bond types 2, and we are in saddle
*     far from the TS, let us get a reduced threshold to avoid
*     wasting our time
*
      Call mma_allocate(Tmp,nSaddle,Label='Tmp')
      Call Get_dArray('Saddle',Tmp,nSaddle)
      Found=.false.
      If (Tmp(nSaddle-1).gt.0.50d0) Then
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
      Call mma_deallocate(Tmp)
*
      Call Free_iWork(ip_TabA)
      Call Free_iWork(ip_TabB)
      Call mma_deallocate(Coor)
      Call mma_deallocate(AN)
      Call mma_deallocate(Vec)
      Call mma_deallocate(TabAI)
      Call mma_deallocate(TR)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
