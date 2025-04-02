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

subroutine FckMat()
!***********************************************************
!                                                          *
!   Driver for calculation of optimized fock matrix.       *
!                                                          *
!***********************************************************

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: FAMO, FIMO, F0SQMO, INT2
use MCLR_Data, only: nDens2
use MCLR_Data, only: nrec
use input_mclr, only: nSym, nAsh, nIsh, nBas, iMethod
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero

implicit none
real*8, allocatable :: Q(:), Tmp2(:,:), T3(:)
integer nm, nmm, nmmm, iS, nAtri

!                                                                      *
!***********************************************************************
!                                                                      *
! Read density matrix

nm = 0
nmm = 0
nmmm = 0
do iS=1,nSym
  nm = nAsh(is)+nm
  nMM = max(nMM,nAsh(is)+nIsh(iS))
  nMMM = max(nmmM,nBas(is))
end do
nAtri = nTri_Elem(nm)
nAtri = nTri_Elem(nAtri)
nmmm = ((nmmm-1)/nRec+1)*nRec
nmm = nmm*nMMM
nmm = nmm**2
call mma_allocate(F0SQMO,ndens2,Label='F0SQMO')
call mma_allocate(FIMO,ndens2,Label='FIMO')
if (iMethod == 2) then
  call mma_allocate(Int2,nAtri,Label='Int2')
else
  call mma_allocate(Int2,1,Label='Int2')
end if
Int2(:) = Zero
call mma_allocate(FAMO,nDens2,Label='FAMO')
call mma_allocate(Q,nDens2,Label='Q')
call mma_allocate(Tmp2,ndens2,2,Label='Tmp2')
call mma_allocate(T3,ndens2,Label='T3')
!                                                                      *
!***********************************************************************
!                                                                      *
! Calculate two-electron contribution

call Read22_2(Int2,F0SQMO,Q,FIMO,FAMO,Tmp2(:,1),Tmp2(:,2),T3)
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(T3)
call mma_deallocate(Tmp2)
call mma_deallocate(Q)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine FckMat
