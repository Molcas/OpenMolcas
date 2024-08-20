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

module ontop_functional
    use definitions, only: wp

    implicit none
    private

    type :: OTFNAL
        real(kind=wp) :: lambda = 0.0d0
        character(len=80) :: otxc = ""
        character(len=80) :: xc = ""
    contains
        procedure :: is_hybrid
    end type

    interface OTFNAL
        procedure :: new
    end interface

    public :: OTFNAL

contains

    type(OTFNAL) function new(otxc, lambda) result(res)
        use mcpdft_output, only: lf
        implicit none

        real(kind=wp), intent(in) :: lambda
        character(len=80), intent(in) :: otxc

        res%lambda = lambda
        res%otxc = otxc
        call upcase(res%otxc)

        if(.not. valid_otxc(res%otxc)) then
            call warningmessage(2, "Wrong on-top functional for MC-PDFT")
            write(lf,*) ' ************* ERROR **************'
            write(lf,*) ' Current on-top functionals are:   '
            write(lf,*) ' T :  translated functionals       '
            write(lf,*) ' FT:  fully translated functionals '
            write(lf,*) ' e.g. T:PBE for tPBE functional    '
            write(lf,*) ' **********************************'
            call abend()
        end if

        res%xc = get_base(res%otxc)

        if(is_hybrid_xc(res%xc)) then
            call warningmessage(2, "Hybrid functionals not supported")
            write(lf,*) ' ************* ERROR **************'
            write(lf,*) ' MC-PDFT does not translate hybrid '
            write(lf,*) ' functionals. If you want to run   '
            write(lf,*) ' hybrid MC-PDFT, use the LAMBda    '
            write(lf,*) ' keyword instead.                  '
            write(lf,*) '                                   '
            write(lf,*) ' EXAMPLE:                          '
            write(lf,*) '  tPBE0 = 75% tPBE + 25% MCSCF.    '
            write(lf,*) ' Usage:                            '
            write(lf,*) '  KSDFT=T:PBE                      '
            write(lf,*) '  LAMB =0.25                       '
            write(lf,*) ' **********************************'
            call abend()
        end if

    end function

    logical function valid_otxc(otxc)
        character(len=80), intent(in) :: otxc
        valid_otxc = otxc(1:2) == "T:" .or. otxc(1:3) == "FT:"
    end function

    logical function is_hybrid_xc(xc_base)
        character(len=80), intent(in) :: xc_base
        real(kind=wp), external :: get_exfac
        is_hybrid_xc = abs(get_exfac(xc_base)) .gt. 1.0d-8
    end function

    logical function is_hybrid(self)
        implicit none
        class(OTFNAL), intent(in) :: self
        is_hybrid =  self%lambda .gt. 0.0d0
    end function

    character(len=80) function get_base(otxc) result(xc)
        implicit none
        character(len=80), intent(in) :: otxc
        xc = otxc(index(otxc, "T:")+2:)
    end function

end module
