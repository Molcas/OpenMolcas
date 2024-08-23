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
      Subroutine DiagMtrx_x(H,nH,iNeg)
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nH, iNeg
      Real*8 H(nH,nH)

      Real*8, Allocatable :: EVal(:), EVec(:,:), Diag(:,:), HU(:,:)
      Real*8 SumHii, Temp
      Integer i, j, ij, ii
!
!     Lu=6
!
      Call mma_allocate(EVal,nH*(nH+1)/2,label='EVal')
      Call mma_allocate(EVec,nH,nH,label='EVec')
!
!---- Copy elements for H
!
      SumHii=Zero
      Do i = 1, nH
         Do j = 1, i
            ij = i*(i-1)/2 + j
            EVal(ij)=H(i,j)
         End Do
         SumHii=SumHii+H(i,i)
      End Do
!     Write (Lu,*) ' SumHii=',SumHii
!
!---- Set up a unit matrix
!
      call dcopy_(nH*nH,[Zero],0,EVec,1)
      call dcopy_(nH,[One],0,EVec,nH+1)
!
!---- Compute eigenvalues and eigenvectors
!
      Call Jacob (EVal,EVec,nH,nH)
      Call Jacord(EVal,EVec,nH,nH)
!
!---- Print out the result
!
      iNeg=0
      Do i = 1, nH
         ii = i*(i+1)/2
         If (EVal(ii).lt.Zero) iNeg=iNeg+1
      End Do
!
      Call mma_allocate(Diag,nH,nH,label='Diag')
      Call mma_allocate(HU,nH,nH,label='HU')
!
      call dcopy_(nH*nH,[Zero],0,Diag,1)
      Do i = 1, nH
         ii=i*(i+1)/2
         temp = EVal(ii)
!        Write (Lu,'(A,G10.4)') 'Hii=',temp
         Diag(i,i)=Max(Abs(temp),1.0D-15)
      End Do
!
      Call DGEMM_('N','N',                                              &
     &            nH,nH,nH,                                             &
     &            1.0d0,EVec,nH,                                        &
     &                  Diag,nH,                                        &
     &            0.0d0,HU,nH)
      Call DGEMM_('N','T',                                              &
     &            nH,nH,nH,                                             &
     &            1.0d0,HU,nH,                                          &
     &                  EVec,nH,                                        &
     &            0.0d0,H,nH)
!
      Call mma_deallocate(HU)
      Call mma_deallocate(Diag)
      Call mma_deallocate(EVec)
      Call mma_deallocate(EVal)
!
      Return
      End Subroutine DiagMtrx_x
