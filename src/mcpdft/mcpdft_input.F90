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
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************

module mcpdft_input
    use definitions, only: iwp
    use ontop_functional, only: OTFNAL_t

    implicit none
    private

    type :: McpdftInputOptions_t
        logical :: wjob = .false.
        logical :: mspdft = .false.
        logical :: grad = .false.
        logical :: meci = .false.
        logical :: nac = .false.
        logical :: is_hdf5_wfn = .true.
        character(len=256) :: wfn_file = "JOBOLD"

        integer(kind=iwp), dimension(2) :: nac_states = 0
        type(OTFNAL_t) :: otfnal

    end type

    type(McpdftInputOptions_t) :: mcpdft_options

    public :: mcpdft_options


end module
