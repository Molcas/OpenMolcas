************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Martin Schuetz                                         *
*               2017, Roland Lindh                                     *
************************************************************************
      Subroutine ErrV(lvec,ivec,QNRstp,ErrVec,HDiag)
************************************************************************
*                                                                      *
*     computes error vector for DIIS, in 1st order case, this          *
*     is just grad(ivec): in 2nd order case (QNRstp.eq..TRUE.)         *
*     this is -H(iterso)*grad(ivec)                                    *
*     the pointer to the proper error vector is returned as function   *
*     val.                                                             *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
*
      Real*8 HDiag(lVec), ErrVec(lVec)
      Integer lvec
      Logical QNRstp
*
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "infso.fh"
#include "stdalloc.fh"
#include "file.fh"
#include "llists.fh"
*
*     local vars
      Integer inode
      Real*8, Dimension(:), Allocatable:: Grad
*
      Call GetNod(ivec,LLGrad,inode)
      If (inode.eq.0) GoTo 555
*
      If (QNRstp) Then
*
*       we eventually need one more vector
*
        Call mma_allocate(Grad,lvec,Label='Grad')
        Call iVPtr(LuGrd,Grad,lvec,inode)
        Call SOrUpV(MemRsv,Grad,HDiag,lvec,ErrVec,'DISP','BFGS')
        Call mma_deallocate(Grad)
*
      Else
*
        Call iVPtr(LuGrd,ErrVec,lvec,inode)
*
      End If
*
      Return
*
*-----Error handling
*
*     Hmmm, no entry found in LList, that's strange
 555  Write (6,*) 'ErrV: no entry found in LList!'
      Call QTrace
      Call Abend()
      End
