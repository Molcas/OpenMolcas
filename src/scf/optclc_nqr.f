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
      Subroutine OptClc_QNR(CInter,nCI,nD,Grd1,Xnp1,nOV,Ind,MxOptm,
     &                      kOptim)
      Implicit None
#include "llists.fh"
#include "file.fh"
#include "stdalloc.fh"
      Integer nCI,nD,nOV,MxOptm,kOptim
      Real*8 CInter(nCI,nD), Grd1(nOV,nD), Xnp1(nOV,nD)
      Integer Ind(MxOptm)
*
      Real*8, Dimension(:,:), Allocatable:: Aux
      Integer inode,ivec,iD,i
*
*-----QNR/DIIS case: compute extrapolated Gradient grd'(n),
*     extrapolated Orb Rot Param x'(n), and from this, the
*     new, extrapolated displacement direction delta(n)
      Call mma_allocate(Aux,nOV,nD,Label='Aux')
      Call FZero(Aux,nOV*nD)
*
*     get last gradient grad(n) from LList
*
      Call GetVec(LuGrd,Ind(kOptim),LLGrad,inode,Grd1,nOV*nD)
      Call GetVec(Lux,  Ind(kOptim),LLx,   inode,Xnp1,nOV*nD)
*
      Do iD = 1, nD
         Call DSCAL_(nOV,CInter(kOptim,iD),Grd1(1,iD),1)
         Call DSCAL_(nOV,CInter(kOptim,iD),Xnp1(1,iD),1)
      End Do
*
      Do i=1,kOptim-1
         ivec=Ind(i)
*
*        get proper gradient from LList.
         Call GetNod(ivec,LLGrad,inode)
         If (inode.eq.0) GoTo 555
         Call iVPtr(LuGrd,Aux,nOV*nD,inode)
         Do iD = 1, nD
            Call Daxpy_(nOV,CInter(i,iD),Aux(1,iD),1,Grd1(1,iD),1)
         End Do
*
*        get proper X-vector from LList.
         Call GetNod(ivec,LLx,inode)
         If (inode.eq.0) GoTo 555
         Call iVPtr(Lux,Aux,nOV*nD,inode)
         Do iD = 1, nD
            Call Daxpy_(nOV,CInter(i,iD),Aux(1,iD),1,Xnp1(1,iD),1)
         End Do
      End Do
*
      Call mma_deallocate(Aux)
*
      Return
*
*-----Error handling
*
*     Hmmm, no entry found in LList, that's strange
 555  Write (6,*) 'DIIS: no entry found in LList!'
      Call QTrace
      Call Abend()
      End
