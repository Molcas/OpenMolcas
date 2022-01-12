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

! This module contains seward parameters which are not supposed to be carried between programs.
module Temporary_Parameters

use Definitions, only: iwp

implicit none
private

integer(kind=iwp) :: IsChi = 0
logical(kind=iwp) :: DirInt = .false., &
                     Expert = .true., &
                     Fake_ERIs = .false., &
                     force_out_of_core = .false., &
                     force_part_c = .false., &
                     force_part_p = .false., &
                     IfAllOrb = .false., &
                     Onenly = .false., &
                     Primitive_Pass = .true., &
                     PrPrt = .false., &
                     Short = .true., &
                     Test = .false.

public :: DirInt, Expert, Fake_ERIs, force_out_of_core, force_part_c, force_part_p, IfAllOrb, IsChi, Onenly, Primitive_Pass, &
          PrPrt, Short, Test

end module Temporary_Parameters
