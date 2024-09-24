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
!
!     Compute the x parameters value as a function of Iter_ref
!
      Subroutine XClc()
      use LnkLst, only: SCF_V, LLx
      use InfSCF, only: Iter, Iter_Start, mOV, Iter_Ref
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer jpgrd,inode, i
      Real*8, Dimension(:), Allocatable:: Scr
      Integer, External :: LstPtr

      Call mma_allocate(Scr,mOV,Label='Scr')

      jpgrd=LstPtr(Iter_Ref,LLx)   ! Pointer to X_old(i_ref)
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'iter=',iter
      Write (6,*) 'iter_Start=',iter_Start
      Write (6,*) 'iter_ref=',iter_ref
      Call NrmClc(SCF_V(jpgrd)%A(:),mOV,'XClc','X(i_ref)(:)')
      Write (6,*)
#endif

!     Loop over all iterations starting at Iter_Start+1

      Do i = Iter_Start, Iter

!        X_new(i)=X_old(i)-X_old(i_ref)

         Call GetNod(i,LLx,inode)
         If (inode.eq.0) Then
            Write (6,*) 'inode.eq.0'
            Call Abend()
         End If
         Call iVPtr(Scr,mOV,inode)
#ifdef _DEBUGPRINT_
         Write (6,*)
         Write (6,*) 'X(i) before  i=',i
         Call NrmClc(Scr(:),mOV,'XClc','Scr(:)')
#endif
         Scr(:)=Scr(:)-SCF_V(jpgrd)%A(:)
#ifdef _DEBUGPRINT_
         Write (6,*) 'X(i) after  i=',i
         Call NrmClc(Scr(:),mOV,'XClc','Scr(:)')
#endif

         Call PutVec(Scr,mOV,i,'OVWR',LLx)
      End Do

      Call mma_deallocate(Scr)

      Return
      End Subroutine XClc
