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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

subroutine CMSRHSGDMat(GDMat)

use Index_Functions, only: iTri, nTri_Elem
use ipPage, only: W
use MCLR_Data, only: ipCI, n1Dens, n2Dens, nNA, XISPSM
use input_mclr, only: nCSF, nRoots, State_Sym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: GDMat(nTri_Elem(nRoots),nnA,nnA)
integer(kind=iwp) :: I, iL, iR, J, nConfL, nConfR, NIJ
real(kind=wp), allocatable :: CIL(:), CIR(:), GDArray(:,:), rdum(:)

! I:state index, "excited state" of state J when I /= J

! (D^IJ)_pq = <I|E_pq|J>, setting I>=J
!  <I|E_pq|J>=<J|E_qp|I>
call mma_allocate(GDArray,nnA,nnA)
call mma_allocate(rdum,n2Dens)
iL = state_sym
iR = state_sym
nConfR = max(ncsf(iR),nint(xispsm(iR,1)))
nConfL = max(ncsf(iL),nint(xispsm(iL,1)))
call mma_allocate(CIR,nConfR)
call mma_allocate(CIL,nConfL)
do I=1,nRoots
  call CSF2SD(W(ipCI)%A(1+(I-1)*ncsf(iL)),CIL,iL)
  do J=1,I
    call CSF2SD(W(ipCI)%A(1+(J-1)*ncsf(iR)),CIR,iR)
    call Densi2_mclr(1,GDArray,rdum,CIL,CIR,0,0,0,n1Dens,n2Dens)
    NIJ = iTri(I,J)
    GDMat(NIJ,:,:) = GDArray(:,:)
  end do
end do
call mma_deallocate(GDArray)
call mma_deallocate(rdum)
call mma_deallocate(CIL)
call mma_deallocate(CIR)

end subroutine CMSRHSGDMat
