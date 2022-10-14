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
*#define _DEBUGPRINT_
      Subroutine OptClc_QNR(CInter,nCI,nD,Grd1,Xnp1,mOV,Ind,MxOptm,
     &                      kOptim,kOV)
      use LnkLst, only: LLGrad,LLx
      Implicit None
#include "file.fh"
#include "stdalloc.fh"
      Integer nCI,nD,mOV,MxOptm,kOptim,kOV(2)
      Integer iSt, iEnd
      Real*8 CInter(nCI,nD), Grd1(mOV), Xnp1(mOV)
      Integer Ind(MxOptm)
*
      Real*8, Allocatable:: Aux(:)
      Integer inode,ivec,iD,i
*
*-----QNR/DIIS case: compute extrapolated Gradient grd'(n),
*     extrapolated Orb Rot Param x'(n), and from this, the
*     new, extrapolated displacement direction delta(n)
      Call mma_allocate(Aux,mOV,Label='Aux')
      Aux(:)=0.0D0
*
*     get last gradient grad(n) from LList
*
      Call GetVec(Ind(kOptim),LLGrad,inode,Grd1,mOV)
      Call GetVec(Ind(kOptim),LLx,   inode,Xnp1,mOV)
*
      iEnd=0
      Do iD = 1, nD
         iSt=iEnd + 1
         iEnd = iEnd + kOV(iD)
         Call DSCAL_(kOV(iD),CInter(kOptim,iD),Grd1(iSt:iEnd),1)
         Call DSCAL_(kOV(iD),CInter(kOptim,iD),Xnp1(iSt:iEnd),1)
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'Initial scalled entities.'
      Call NrmClc(Grd1(:),mOV,'OptClc_qNR','Grd1')
      Call NrmClc(Xnp1(:),mOV,'OptClc_qNR','Xnp1')
      Write (6,*)
#endif
*
      Do i=1,kOptim-1
         ivec=Ind(i)
*
*        get proper gradient from LList.
         Call GetNod(ivec,LLGrad,inode)
         If (inode.eq.0) GoTo 555
         Call iVPtr(Aux,mOV,inode)
         iEnd = 0
         Do iD = 1, nD
            iSt=iEnd + 1
            iEnd = iEnd + kOV(iD)
            Call Daxpy_(kOV(iD),CInter(i,iD),Aux(iSt:iEnd),1,
     &                                      Grd1(iSt:iEnd),1)
         End Do
*
*        get proper X-vector from LList.
         Call GetNod(ivec,LLx,inode)
         If (inode.eq.0) GoTo 555
         Call iVPtr(Aux,mOV,inode)
         iEnd = 0
         Do iD = 1, nD
            iSt=iEnd + 1
            iEnd = iEnd + kOV(iD)
            Call Daxpy_(kOV(iD),CInter(i,iD),Aux(iSt:iEnd),1,
     &                                      Xnp1(iSt:iEnd),1)
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*)
      Call NrmClc(Grd1(:),mOV,'OptClc_qNR','Grd1')
      Call NrmClc(Xnp1(:),mOV,'OptClc_qNR','Xnp1')
      Write (6,*)
#endif
*
      Call mma_deallocate(Aux)
*
      Return
*
*-----Error handling
*
*     Hmmm, no entry found in LList, that's strange
 555  Write (6,*) 'DIIS: no entry found in LList!'
      Call Abend()
      End
