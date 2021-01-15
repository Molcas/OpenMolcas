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
Module ChoArr
Implicit none
Private
Public:: iSOShl, iShlSO, iBasSh, nBasSh, nBstSh, iSP2F, iAtomShl, iRS2F, IntMap, iScr
Public:: nDimRS
Integer, Allocatable:: iSOShl(:)
Integer, Allocatable:: iShlSO(:)
Integer, Allocatable:: iBasSh(:,:)
Integer, Allocatable:: nBasSh(:,:)
Integer, Allocatable:: nBstSh(:)
Integer, Allocatable:: iSP2F(:)
Integer, Allocatable:: iAtomShl(:)
Integer, Allocatable:: iRS2F(:,:)
Integer, Allocatable:: IntMap(:)
Integer, Allocatable:: iScr(:)
Integer, Allocatable:: nDimRS(:,:)
End Module ChoArr
