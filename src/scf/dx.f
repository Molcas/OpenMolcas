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
*     Compute the dx parameters value
*
      Subroutine dX()
      use LnkLst, only: SCF_V, LLx, LLDelt
      use InfSCF, only: Iter, Iter_Start, mOV
      Implicit None
#include "real.fh"
#include "stdalloc.fh"
#include "file.fh"
      Integer jpgrd,inode, i
      Real*8, Dimension(:), Allocatable:: Scr
      Integer, External :: LstPtr

      Call mma_allocate(Scr,mOV,Label='Scr')

!     Loop over all iterations starting at Iter_Start

      Do i = Iter_Start, Iter-1

!        dX(i)=X(i+1)-X(i)

         jpgrd=LstPtr(i+1,LLx)   ! Pointer to X(i+1)

         Call GetNod(i,LLx,inode) ! X(i)
         If (inode.eq.0) Then
            Write (6,*) 'inode.eq.0'
            Call Abend()
         End If
         Call iVPtr(Scr,mOV,inode)

         Scr(:)=SCF_V(jpgrd)%A(:)-Scr(:)

         Call PutVec(Scr,mOV,i,'OVWR',LLDelt)
      End Do

      Call mma_deallocate(Scr)

      Return
      End Subroutine dX
