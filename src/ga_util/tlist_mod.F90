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

module TList_Mod

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: iEnd_TList, igaTsk, iStrt_TList, iTCnSt, iTskCan, mTasks, nTasks
real(kind=wp) :: P, PQ, QLast(2)
logical(kind=iwp) :: GT_Status = .false., PP_Status = .false.
integer(kind=iwp), allocatable, target :: TskL(:)
real(kind=wp), allocatable :: TskM(:,:), TskQ(:,:)
real(kind=wp), parameter :: Not_Used = -One

public :: GT_Status, iEnd_TList, igaTsk, iStrt_TList, iTCnSt, iTskCan, mTasks, Not_Used, nTasks, P, PP_Status, PQ, QLast, TskL, &
          TskM, TskQ

end module TList_Mod
