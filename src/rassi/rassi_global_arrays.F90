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
      Module RASSI_GLOBAL_ARRAYS
      Real*8, Allocatable:: HAM(:,:), SFDYS(:,:,:), &
                            SODYSAMPS(:,:), SODYSAMPSR(:,:), &
                            SODYSAMPSI(:,:), EIGVEC(:,:), &
                            PROP(:,:,:), ESHFT(:), HDIAG(:)
      Integer, Allocatable:: JBNUM(:), LROOT(:)
      Integer, Allocatable:: PART(:)
      Integer, Allocatable:: OrbTab(:)
      Integer, Allocatable:: SSTAB(:)
      Integer, Allocatable, Target:: REST1(:), REST2(:)
      Integer, Pointer:: REST(:)
      Integer, Allocatable, Target:: CNFTAB1(:), CNFTAB2(:)
      Integer, Pointer:: CNFTAB(:)
      Integer, Allocatable, Target:: FSBTAB1(:), FSBTAB2(:)
      Integer, Pointer:: FSBTAB(:)
      Integer, Allocatable:: FSBARR(:)
      Integer, Allocatable, Target:: SPNTAB1(:), SPNTAB2(:)
      Integer, Pointer:: SPNTAB(:)
      Real*8, Allocatable, Target:: TRANS1(:), TRANS2(:)
      Real*8, Pointer:: TRANS(:)
      Integer, Allocatable, Target:: FSBANN1(:), FSBANN2(:)
      Integer, Allocatable, Target:: FSBANN3(:), FSBANN4(:)
      Integer, Pointer:: FSBANN(:)
      End Module RASSI_GLOBAL_ARRAYS
