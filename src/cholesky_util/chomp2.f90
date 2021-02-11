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
Module ChoMP2
Implicit none
Private
Public:: ChoMP2_allocated, ChoMP2g_allocated, OldVec, EFrozT, EOccuT, EVirtT
Public:: AdrR1, AdrR2

Logical:: ChoMP2_allocated=.FALSE.
Logical:: ChoMP2g_allocated=.FALSE.
Real*8, Allocatable:: OldVec(:)

Type Pointers
     Real*8, Pointer:: A(:,:)=>Null()
End Type Pointers

Real*8, Allocatable:: EFrozT(:)
Real*8, Allocatable:: EOccuT(:)
Real*8, Allocatable:: EVirtT(:)

Integer, Allocatable:: AdrR1(:,:,:), AdrR2(:,:,:)
End Module ChoMP2
