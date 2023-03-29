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
      Subroutine OptClc_X(CInter,nCI,nD,Array,mOV,Ind,MxOptm,
     &                    kOptim,kOV,LL)
      use Files
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nCI,nD,mOV,MxOptm,kOptim,kOV(2), LL
      Integer iSt, iEnd
      Real*8 CInter(nCI,nD), Array(mOV)
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
*     get last array(n) from LL
*
      Call GetVec(Ind(kOptim),LL,inode,Array,mOV)
*
      iEnd=0
      Do iD = 1, nD
         iSt=iEnd + 1
         iEnd = iEnd + kOV(iD)
         Call DSCAL_(kOV(iD),CInter(kOptim,iD),Array(iSt:iEnd),1)
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'Initial scaled entities.'
      Call NrmClc(Array(:),mOV,'OptClc_X','Array')
      Write (6,*)
#endif
*
      Do i=1,kOptim-1
         ivec=Ind(i)
*
*        get proper gradient from LList.
         Call GetNod(ivec,LL,inode)
         If (inode.eq.0) Then
            Write (6,*) 'DIIS: no entry found in LList!'
            Call Abend()
         End If
         Call iVPtr(Aux,mOV,inode)
         iEnd = 0
         Do iD = 1, nD
            iSt=iEnd + 1
            iEnd = iEnd + kOV(iD)
            Call Daxpy_(kOV(iD),CInter(i,iD),Aux(iSt:iEnd),1,
     &                                      Array(iSt:iEnd),1)
         End Do
*
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*)
      Call NrmClc(Array(:),mOV,'OptClc_X','Array')
      Write (6,*)
#endif
*
      Call mma_deallocate(Aux)
*
      Return
      End Subroutine OptClc_X
