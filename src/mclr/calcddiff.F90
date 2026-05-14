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

subroutine CalcDdiff(Ddiff,GDMat,M,K,nnA,nRoots)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: M, K, nnA, nRoots
real(kind=wp), intent(out) :: Ddiff(nnA,nnA)
real(kind=wp), intent(in) :: GDMat(nTri_Elem(nRoots),nnA,nnA)
integer(kind=iwp) :: iKK, iMM, it

iMM = nTri_Elem(M)
iKK = nTri_Elem(K)

do it=1,nnA
  Ddiff(:,it) = GDMat(iMM,it,:)-GDMat(iKK,it,:)
end do

end subroutine CalcDdiff
