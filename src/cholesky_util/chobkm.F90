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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************

! info for Cholesky bookmarks
module ChoBkm

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: nCol_BkmThr = 0, nCol_BkmVec = 0, nRow_BkmThr = 0, nRow_BkmVec = 0
integer(kind=iwp), allocatable :: BkmVec(:,:)
real(kind=wp), allocatable :: BkmThr(:,:)

public :: BkmThr, BkmVec, nCol_BkmThr, nCol_BkmVec, nRow_BkmThr, nRow_BkmVec

end module ChoBkm
