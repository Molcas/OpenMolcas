!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module RASSCF_LUCIA

private

integer, public :: kvec3_length = 0, ini_h0, Memory_Needed_Lucia = 0
logical, public :: Sigma_on_disk = .false.
real*8, allocatable, public :: CIVec(:)
real*8, allocatable, public :: PAtmp(:)
real*8, allocatable, public :: Pscr(:)
real*8, allocatable, public :: Ptmp(:)
real*8, allocatable, public :: DStmp(:)
real*8, allocatable, public :: Dtmp(:)
real*8, allocatable, public :: RF1(:)
real*8, allocatable, public :: RF2(:)

end module RASSCF_LUCIA
