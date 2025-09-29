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

subroutine QaaVerif(G2q,ng2,PUVX,NPUVX,IndTUVX)

use MCLR_Data, only: nNA
use input_mclr, only: ntAsh
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nG2, nPUVX, IndTUVX(ntAsh,ntAsh,ntAsh,ntAsh)
real(kind=wp), intent(in) :: G2q(nG2), PUVX(NPUVX)
integer(kind=iwp) :: I, IJKL, J, K, L, lMax
real(kind=wp) :: dQdX

ijkl = 0
dQdX = Zero
do i=1,nna
  do j=1,i
    do k=1,i
      if (i == k) then
        lmax = j
      else
        lmax = k
      end if
      do l=1,lmax
        ijkl = ijkl+1
        dQdX = dQdX+G2q(ijkl)*PUVX(IndTUVX(I,J,K,L))
      end do
    end do
  end do
end do

write(u6,*) 'dQdX in QaaVerif=',dQdX

end subroutine QaaVerif
