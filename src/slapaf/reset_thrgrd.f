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
      Subroutine Reset_ThrGrd(nIter,mTtAtm,ThrGrd)
      use Slapaf_Info, only: Cx
      use Slapaf_Parameters, only: nDimBC, rHidden
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
      Logical Found
      Integer, Allocatable:: TabAI(:), AN(:)
      Real*8, Allocatable:: TR(:), Vec(:), Coor(:,:), Tmp(:)
      Integer, Allocatable:: TabB(:,:), TabA(:,:,:)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
        Subroutine Box(Coor,nsAtom,iANr,TabB,TabA,nBonds,nMax)
        Integer nsAtom
        Real*8 Coor(3,nsAtom)
        Integer iANr(nsAtom)
        Integer, Allocatable:: TabB(:,:), TabA(:,:,:)
        Integer nBonds, nMax
        End Subroutine Box
        Subroutine Hidden(mTtAtm,Coor,AN,nHidden)
        Integer mTtAtm
        Real*8, Allocatable:: Coor(:,:)
        Integer, Allocatable:: AN(:)
        Integer nHidden
        End Subroutine Hidden
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      nsAtom=SIZE(Cx,2)
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
      Call mma_allocate(TR,18*nsAtom,Label='TR')
      TR(:)=Zero
*
      Call TRPGen(nDimBC,nsAtom,Cx(1,1,iIter),mTR,.False.,TR)
*
*     Call RecPrt('TR',' ',TR,nDimBC,mTR)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_Allocate(TabAI,2*mTtAtm,Label='TabAI')
      Call mma_Allocate(Vec,3*mTtAtm*nDimBC,Label='Vec')
      Call mma_Allocate(AN,mTtAtm,Label='AN')
      Call mma_Allocate(Coor,3,mTtAtm,Label='Coor')
*
*-----Generate Grand atoms list
*
      Call GenCoo(Cx(1,1,iIter),nsAtom,Coor,mTtAtm,Vec,nDimBC,AN,TabAI)
*
*---- Are there some hidden frozen atoms ?
*
      nHidden = 0
      If (rHidden.ge.Two) Call Hidden(mTtAtm,Coor,AN,nHidden)
*
*-----Generate bond list
*
      mTtAtm = mTtAtm+nHidden
      Call Box(Coor,mTtAtm,AN,TabB,TabA,nBonds,nMax)
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
         Do i=1,nBonds
            If (TabB(3,i).eq.2) Then
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
      Call mma_deallocate(TabA)
      Call mma_deallocate(TabB)
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
