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
! Copyright (C) 2019, Roland Lindh                                     *
!***********************************************************************

module RASSI_GLOBAL_ARRAYS

implicit none

real*8, allocatable :: HAM(:,:), SFDYS(:,:,:), SODYSAMPS(:,:), SODYSAMPSR(:,:), SODYSAMPSI(:,:), EIGVEC(:,:), PROP(:,:,:), &
                       ESHFT(:), HDIAG(:)
integer, allocatable :: JBNUM(:), LROOT(:)
integer, allocatable :: PART(:)
integer, allocatable :: OrbTab(:)
integer, allocatable :: SSTAB(:)
integer, allocatable, target :: REST1(:), REST2(:)
integer, pointer :: REST(:)
integer, allocatable, target :: CNFTAB1(:), CNFTAB2(:)
integer, pointer :: CNFTAB(:)
integer, allocatable, target :: FSBTAB1(:), FSBTAB2(:)
integer, pointer :: FSBTAB(:)
integer, allocatable :: FSBARR(:)
integer, allocatable, target :: SPNTAB1(:), SPNTAB2(:)
integer, pointer :: SPNTAB(:)
real*8, allocatable, target :: TRANS1(:), TRANS2(:)
real*8, pointer :: TRANS(:)
integer, allocatable, target :: FSBANN1(:), FSBANN2(:)
integer, allocatable, target :: FSBANN3(:), FSBANN4(:)
integer, pointer :: FSBANN(:)

end module RASSI_GLOBAL_ARRAYS
