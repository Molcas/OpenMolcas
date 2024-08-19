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
    use definitions, only: iwp, wp
    implicit none
    private

    type :: McpdftInputOptions
        logical :: wjob = .false.
        logical :: mspdft = .false.
        logical :: grad = .false.
        logical :: meci = .false.
        logical :: nac = .false.
        integer(kind=iwp), dimension(2) :: nac_states = 0
        real(kind=wp) :: lambda = 0.0d0

        ! This should replace ksdft in the rasscf.fh
        !character, dimension(80) :: ksdft

        ! This should replace the StartOrbFile in general.fh
        logical :: is_hdf5_wfn = .true.
        !character, dimension(256) :: orbital_file
    contains
        procedure :: do_hybrid
    end type

    type(McpdftInputOptions) :: mcpdft_options = McpdftInputOptions()

    public :: mcpdft_options

contains

    logical function do_hybrid(self) result(res)
        class(McpdftInputOptions), intent(in) :: self
        res = self%lambda .gt. 0.0d0
    end function



end module
