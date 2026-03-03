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

real*8, parameter :: Not_Used = -1.0d0
real*8 QLast(2), P, PQ
integer nTasks, igaTsk, iTCnSt, mTasks, iStrt_TList, iEnd_TList, iTskCan
real*8, allocatable :: TskQ(:,:)
real*8, allocatable :: TskM(:,:)
integer, allocatable, target :: TskL(:)
logical :: PP_Status = .false.
logical :: GT_Status = .false.

end module TList_Mod
