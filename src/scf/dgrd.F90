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
!
!     Compute the difference between consecutive gradients
!
      Subroutine dGrd()
      use LnkLst, only: SCF_V, LLGrad, LLdGrd
      use InfSCF, only: iter, Iter_Start, mOV
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer jpgrd,inode, i
      Real*8, Dimension(:), Allocatable:: Scr
      Integer, External :: LstPtr

      Call mma_allocate(Scr,mOV,Label='Scr')

!     Loop over all iterations starting at Iter_Start+1

      Do i = Iter_Start+1, Iter

!        dg(i-1)=g(i)-g(i-1)

         jpgrd=LstPtr(i,LLGrad)
         Call GetNod(i-1,LLGrad,inode)
         If (inode.eq.0) Then
            Write (6,*) 'inode.eq.0'
            Call Abend()
         End If
         Call iVPtr(Scr,mOV,inode)

         Scr(:) =  SCF_V(jpgrd)%A(:) - Scr (:)

         Call PutVec(Scr,mOV,i-1,'OVWR',LLdGrd)
      End Do

      Call mma_deallocate(Scr)

      Return
      End Subroutine dGrd
