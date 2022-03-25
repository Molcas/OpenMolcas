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
*
*     Compute the difference with the previous gradient
*
      Subroutine dGrd()
      use LnkLst, only: SCF_V
      Implicit None
#include "mxdm.fh"
#include "real.fh"
#include "infscf.fh"
#include "stdalloc.fh"
#include "file.fh"
#include "llists.fh"
      Integer nD,jpgrd,inode
      Real*8, Dimension(:,:), Allocatable:: Scr
      Integer, External :: LstPtr
      If (iter.eq.1) Return
      nD=iUHF+1
      Call mma_allocate(Scr,nOV,nD,Label='Scr')
      jpgrd=LstPtr(iter,LLGrad)
      Call GetNod(iter-1,LLGrad,inode)
      If (inode.eq.0) Then
         Write (6,*) 'inode.eq.0'
         Call Abend()
      End If
      Call iVPtr(Scr,nOV*nD,inode)
      Call DaXpY_(nOV*nD,-One,SCF_V(jpgrd)%A,1,Scr,1)
      Call DScal_(nOV*nD,-One,Scr,1)
      Call PutVec(Scr,nOV*nD,iter-1,'NOOP',LLdGrd)
      Call mma_deallocate(Scr)
      Return
      End Subroutine dGrd
