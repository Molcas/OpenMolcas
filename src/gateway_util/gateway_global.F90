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

module Gateway_global

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), parameter :: G_Mode = 1, S_Mode = 2, GS_Mode = 3
integer(kind=iwp) :: iPack = 0, IsChi = 0, Run_Mode
logical(kind=iwp) :: asymptotic_Rys = .false., &
                     DirInt = .false., &
                     Expert = .true., &
                     Fake_ERIs = .false., &
                     FMM_shortrange = .false., &
                     force_out_of_core = .false., &
                     force_part_c = .false., &
                     force_part_p = .false., &
                     IfAllOrb = .false., &
                     NoTab = .false., &
                     Onenly = .false., &
                     Primitive_Pass = .true., &
                     PrPrt = .false., &
                     Short = .true., &
                     Test = .false.
character(len=512) :: SW_FileOrb = 'INPORB'

public :: asymptotic_Rys, DirInt, Expert, Fake_ERIs, FMM_shortrange, force_out_of_core, force_part_c, force_part_p, G_Mode, &
          GS_Mode, IfAllOrb, iPack, IsChi, NoTab, Onenly, Primitive_Pass, PrPrt, Run_Mode, S_Mode, Short, SW_FileOrb, Test

end module Gateway_global
