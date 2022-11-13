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
!***********************************************************************
!                                                                      *
! Module file for the peek/poke utilities.                             *
!                                                                      *
!***********************************************************************

module peekpoke

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: nTabDS = 32, nTabIS = 32
character(len=24) :: ds_label(nTabDS), is_label(nTabIS)
real(kind=wp) :: ds_value(nTabIS)
integer(kind=iwp) :: is_value(nTabDS)
integer(kind=iwp) :: ds_no, is_no

public :: ds_label, ds_value, ds_no, is_label, is_value, is_no

end module peekpoke
