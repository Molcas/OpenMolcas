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
Module TList_Mod
Real*8, Parameter:: Not_Used=-1.0D0
Real*8 QLast(2),P,PQ
Integer nTasks, igaTsk, iTCnSt, mTasks, iStrt_TList, iEnd_TList, iTskCan
Real*8, Allocatable:: TskQ(:,:)
Real*8, Allocatable:: TskM(:,:)
Integer, Allocatable, Target:: TskL(:)
Logical:: PP_Status=.False.
End Module TList_Mod
