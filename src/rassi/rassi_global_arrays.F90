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

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), allocatable :: FSBARR(:), JBNUM(:), LROOT(:), OrbTab(:), PART(:), SSTAB(:)
integer(kind=iwp), allocatable, target :: CNFTAB1(:), CNFTAB2(:), FSBANN1(:), FSBANN2(:), FSBANN3(:), FSBANN4(:), FSBTAB1(:), &
                                          FSBTAB2(:), REST1(:), REST2(:), SPNTAB1(:), SPNTAB2(:)
real(kind=wp), allocatable :: EIGVEC(:,:), ESHFT(:), HAM(:,:), HDIAG(:), SFDYS(:,:,:), SODYSAMPS(:,:), SODYSAMPSI(:,:), &
                              SODYSAMPSR(:,:)
real(kind=wp), allocatable, target :: TRANS1(:), TRANS2(:)

public :: CNFTAB1, CNFTAB2, EIGVEC, ESHFT, FSBANN1, FSBANN2, FSBANN3, FSBANN4, FSBARR, FSBTAB1, FSBTAB2, HAM, HDIAG, JBNUM, LROOT, &
          OrbTab, PART, REST1, REST2, SFDYS, SODYSAMPS, SODYSAMPSI, SODYSAMPSR, SPNTAB1, SPNTAB2, SSTAB, TRANS1, TRANS2

end module RASSI_GLOBAL_ARRAYS
