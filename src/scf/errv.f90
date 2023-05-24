!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Martin Schuetz                                         *
!               2017, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
Subroutine ErrV(lvec,ivec,QNRstp,ErrVec)
!***********************************************************************
!                                                                      *
!     computes error vector for DIIS, in 1st order case, this          *
!     is just grad(ivec): in 2nd order case (QNRstp.eq..TRUE.)         *
!     this is -H(iterso)*grad(ivec)                                    *
!     the pointer to the proper error vector is returned as function   *
!     val.                                                             *
!                                                                      *
!***********************************************************************
use LnkLst, only: LLGrad
use stdalloc, only: mma_allocate, mma_deallocate
Implicit None
!
Integer lVec, iVec
Real*8 ErrVec(lVec)
Logical QNRstp
!
!     local vars
Integer inode
Real*8, Allocatable:: Grad(:)
!
Call GetNod(ivec,LLGrad,inode)
If (inode.eq.0) Then
!        Hmmm, no entry found in LList, that's strange
   Write (6,*) 'ErrV: no entry found in LList!'
   Call Abend()
End If
!
If (QNRstp) Then
!
!       for qNR step compute delta = - H^{-1}g
!
  Call mma_allocate(Grad,lvec,Label='Grad')
  Call iVPtr(Grad,lvec,inode)
#ifdef _DEBUGPRINT_
  Write (6,*) 'ErrV, iVec=',iVec
  Call NrmClc(Grad,lVec,'ErrV','Grad')
#endif
  Call SOrUpV(Grad,lvec,ErrVec,'DISP','BFGS')
#ifdef _DEBUGPRINT_
  Call NrmClc(ErrVec,lVec,'ErrV','ErrVec')
#endif
  Call mma_deallocate(Grad)
!
Else
!
!       Pick up the gradient
!
  Call iVPtr(ErrVec,lvec,inode)
!
End If
!
Return
!
End Subroutine ErrV
