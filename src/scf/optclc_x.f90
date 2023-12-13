!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!#define _DEBUGPRINT_
      Subroutine OptClc_X(CInter,nCI,nD,Array,mOV,Ind,MxOptm,kOptim,kOV,LL,DD)
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only:Zero
      Implicit None
      Integer nCI,nD,mOV,MxOptm,kOptim,kOV(2), LL
      Real*8 CInter(nCI,nD), Array(mOV)
      Integer Ind(MxOptm)
      Real*8, Optional:: DD
!
      Integer iSt, iEnd
      Real*8, Allocatable:: Aux(:)
      Integer inode,ivec,iD,i
!
!-----QNR/DIIS case: compute extrapolated Gradient grd'(n),
!     extrapolated Orb Rot Param x'(n), and from this, the
!     new, extrapolated displacement direction delta(n)
      Call mma_allocate(Aux,mOV,Label='Aux')
      Aux(:)=Zero
!
!     get last array(n) from LL
!
      Call GetVec(Ind(kOptim),LL,inode,Array,mOV)
!
      iEnd=0
      Do iD = 1, nD
         iSt=iEnd + 1
         iEnd = iEnd + kOV(iD)
         Array(iSt:iEnd)=CInter(kOptim,iD)*Array(iSt:iEnd)
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'Initial scaled entities.'
      Call NrmClc(Array(:),mOV,'OptClc_X','Array')
      Write (6,*)
#endif
!
      Do i=1,kOptim-1
         ivec=Ind(i)
!
!        get proper gradient from LList.
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
            Call Daxpy_(kOV(iD),CInter(i,iD),Aux(iSt:iEnd),1,Array(iSt:iEnd),1)
         End Do
!
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*)
      Call NrmClc(Array(:),mOV,'OptClc_X','Array')
      Write (6,*)
#endif
      If (Present(DD)) Then
         DD = Zero
         iEnd = 0
         Do iD = 1, nD
            iSt=iEnd + 1
            iEnd = iEnd + kOV(iD)
            Do i = iSt, iEnd
               DD = DD + Array(i)**2
            End Do
         End Do
         DD = Sqrt(DD)
      End If
!
      Call mma_deallocate(Aux)
!
      Return
      End Subroutine OptClc_X
