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
Public:: Pointer_2D
Public:: ChoMP2_allocated, ChoMP2g_allocated, OldVec, EFrozT, EOccuT, EVirtT
Public:: AdrR1, AdrR2
Public:: MP2W_full, MP2W
Public:: MP2D_full, MP2D
Public:: MP2W_e_full, MP2W_e
Public:: MP2D_e_full, MP2D_e

Public:: iFirst, iFirstS, NumOcc, LnOcc, LnT1am, LiT1am, LnMatij, LiMatij
Public:: lUnit, NumBatOrb, LnBatOrb
Public:: LnPQprod, LiPQprod

Logical:: ChoMP2_allocated=.FALSE.
Logical:: ChoMP2g_allocated=.FALSE.
Real*8, Allocatable:: OldVec(:)

Type Pointer_2D
     Real*8, Pointer:: A(:,:)=>Null()
End Type Pointer_2D

Real*8, Allocatable, Target:: MP2D_full(:)
Type (Pointer_2D):: MP2D(8)
Real*8, Allocatable, Target:: MP2W_full(:)
Type (Pointer_2D):: MP2W(8)
Real*8, Allocatable, Target:: MP2D_e_full(:)
Type (Pointer_2D):: MP2D_e(8)
Real*8, Allocatable, Target:: MP2W_e_full(:)
Type (Pointer_2D):: MP2W_e(8)

Real*8, Allocatable:: EFrozT(:)
Real*8, Allocatable:: EOccuT(:)
Real*8, Allocatable:: EVirtT(:)

Integer, Allocatable:: AdrR1(:,:,:), AdrR2(:,:,:)

Integer, Allocatable:: iFirst(:), iFirstS(:,:), NumOcc(:), LnOcc(:,:)
Integer, Allocatable:: LnT1am(:,:), LiT1am(:,:,:)
Integer, Allocatable:: LnMatij(:,:), LiMatij(:,:,:)
Integer, Allocatable:: lUnit(:,:)
Integer, Allocatable:: NumBatOrb(:), LnBatOrb(:,:)
Integer, Allocatable:: LnPQprod(:,:), LiPQprod(:,:,:)
End Module ChoMP2
