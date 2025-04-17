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
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine QaaP2MO(G2q,ng2,GDMat,IKL,IKK,ILL)

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: nNA
use input_mclr, only: nRoots
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Two, Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nG2, IKL, IKK, ILL
real(kind=wp), intent(out) :: G2q(nG2)
real(kind=wp), intent(in) :: GDMat(nTri_Elem(nRoots),nnA,nnA)
integer(kind=iwp) :: i, ij, ijkl, j, k, kl, l, lMax, nD
real(kind=wp) :: Fact
real(kind=wp), allocatable :: Ddif(:), Dsum(:)

! Calculating dQ_aa/dX_KL, original purpose of this subroutine
nD = nTri_Elem(nnA)
call mma_allocate(Dsum,nD)
call mma_allocate(Ddif,nD)
do i=1,nnA
  do j=1,i
    ij = iTri(i,j)
    Dsum(ij) = GDMat(IKL,i,j)+GDMat(IKL,j,i)
    Ddif(ij) = GDMat(IKK,i,j)-GDMat(ILL,i,j)
  end do
end do
ijkl = 0
do i=1,nna
  do j=1,i
    ij = iTri(i,j)
    do k=1,i
      if (i == k) then
        lmax = j
      else
        lmax = k
      end if
      do l=1,lmax
        kl = iTri(k,l)
        ijkl = ijkl+1
        fact = Half
        if (k == l) fact = Quart
        G2q(ijkl) = fact*(Dsum(ij)*Ddif(kl)+Dsum(kl)*Ddif(ij))*Two
      end do
    end do
  end do
end do
call mma_deallocate(Dsum)
call mma_deallocate(Ddif)

end subroutine QaaP2MO
