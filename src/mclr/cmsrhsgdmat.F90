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

use CI_Info, only: CI_St
use Index_Functions, only: iTri, nTri_Elem
use ipPage, only: W
use MCLR_Data, only: ipCI, n1Dens, n2Dens, nNA
use input_mclr, only: nCSF, nRoots
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

use CI_Operations, only: CI_Trans_Dens
implicit none
integer(kind=iwp), parameter :: iState = 1
real(kind=wp), intent(out) :: GDMat(nTri_Elem(nRoots),nnA,nnA)
integer(kind=iwp) :: I, iL, iR, J, NIJ
real(kind=wp), allocatable :: GDArray(:,:), rdum(:)

! I:state index, "excited state" of state J when I /= J

! (D^IJ)_pq = <I|E_pq|J>, setting I>=J
!  <I|E_pq|J>=<J|E_qp|I>
call mma_allocate(GDArray,nnA,nnA)
call mma_allocate(rdum,n2Dens)
iL = CI_St(iState)%STSym
iR = CI_St(iState)%STSym
do I=1,nRoots
  do J=1,I
    call CI_Trans_Dens(iState,1,0,iL,iR,ncsf(iL),W(ipCI)%A(1+(I-1)*ncsf(iL):I*ncsf(iL)),ncsf(iR), &
                       W(ipCI)%A(1+(J-1)*ncsf(iR):J*ncsf(iR)),n1Dens,GDArray,n2Dens,rdum)
    NIJ = iTri(I,J)
    GDMat(NIJ,:,:) = GDArray(:,:)
  end do
end do
call mma_deallocate(GDArray)
call mma_deallocate(rdum)

end subroutine CMSRHSGDMat
