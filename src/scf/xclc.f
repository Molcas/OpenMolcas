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
*     Compute the x parameters value as a function of Iter_ref
*
      Subroutine XClc()
      use LnkLst, only: SCF_V, LLx
      use InfSCF, only: Iter, Iter_Start, mOV, Iter_Ref
      Implicit None
#include "real.fh"
#include "stdalloc.fh"
#include "file.fh"
      Integer jpgrd,inode, i
      Real*8, Dimension(:), Allocatable:: Scr, Scr_ref
      Integer, External :: LstPtr

      Call mma_allocate(Scr,mOV,Label='Scr')

      jpgrd=LstPtr(Iter_Ref,LLx)   ! Pointer to X_old(i_ref)
*     Write (*,*) 'iter=',iter
*     Write (*,*) 'iter_Start=',iter_Start
*     Write (*,*) 'iter_ref=',iter_ref
*     Call RecPrt('x_old(i_Ref)',' ',SCF_V(jpgrd)%A,1,mOV)

!     Loop over all iterations starting at Iter_Start+1

      Do i = Iter_Start, Iter

!        X_new(i)=X_old(i)-X_old(i_ref)

         Call GetNod(i,LLx,inode)
         If (inode.eq.0) Then
            Write (6,*) 'inode.eq.0'
            Call Abend()
         End If
         Call iVPtr(Scr,mOV,inode)

         Call DaXpY_(mOV,-One,SCF_V(jpgrd)%A,1,Scr,1)
         Call PutVec(Scr,mOV,i,'OVWR',LLx)
      End Do

      Call mma_deallocate(Scr)

      Return
      End Subroutine XClc
