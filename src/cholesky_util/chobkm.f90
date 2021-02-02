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
Module ChoBkm
Implicit none
Private
Public:: BkmVec, BkmThr, nRow_BkmVec, nCol_BkmVec, nRow_BkmThr, nCol_BkmThr
Integer, Allocatable:: BkmVec(:,:)
Real*8, Allocatable:: BkmThr(:,:)
Integer:: nRow_BkmVec=0, nCol_BkmVec=0
Integer:: nRow_BkmThr=0, nCol_BkmThr=0
End Module ChoBkm
