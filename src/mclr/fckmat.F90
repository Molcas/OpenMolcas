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
use MCLR_Data, only: F0SQMO, FAMO, FIMO, INT2, nDens, nrec
use input_mclr, only: iMethod, nAsh, nBas, nIsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nAtri, nm, nmm, nmmm
real(kind=wp), allocatable :: Q(:), T3(:), Tmp2(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Read density matrix

nm = sum(nAsh(1:nSym))
nmm = max(0,maxval(nAsh(1:nSym)+nIsh(1:nSym)))
nmmm = max(0,maxval(nBas(1:nSym)))
nAtri = nTri_Elem(nm)
nAtri = nTri_Elem(nAtri)
nmmm = ((nmmm-1)/nRec+1)*nRec
nmm = nmm*nMMM
nmm = nmm**2
call mma_allocate(F0SQMO,nDens,Label='F0SQMO')
call mma_allocate(FIMO,nDens,Label='FIMO')
if (iMethod == 2) then
  call mma_allocate(Int2,nAtri,Label='Int2')
else
  call mma_allocate(Int2,1,Label='Int2')
end if
Int2(:) = Zero
call mma_allocate(FAMO,nDens,Label='FAMO')
call mma_allocate(Q,nDens,Label='Q')
call mma_allocate(Tmp2,nDens,2,Label='Tmp2')
call mma_allocate(T3,nDens,Label='T3')
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
