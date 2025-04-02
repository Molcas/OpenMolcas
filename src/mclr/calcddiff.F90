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

implicit none
integer nnA, nRoots, M, K
real*8, dimension(nTri_Elem(nRoots),nnA,nnA) :: GDMat
real*8, dimension(nnA**2) :: Ddiff
integer it, iu, iMM, iKK

iMM = nTri_Elem(M)
iKK = nTri_Elem(K)

do it=1,nnA
  do iu=1,nnA
    Ddiff((it-1)*nnA+iu) = GDMat(iMM,it,iu)-GDMat(iKK,it,iu)
  end do
end do

end subroutine CalcDdiff
