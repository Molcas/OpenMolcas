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

module Pack_mod

use Definitions, only: wp, iwp

implicit none
private

!----------------------------------------------------------------------*
!                                                                      *
!     Entries of the packing table:                                    *
!                                                                      *
!     PkThrs : desired accuracy of packing                             *
!     isPack : Flag to indicate desired mode of action                 *
!              (isPack=.true. : packing is desired )                   *
!              (isPack=.false.: no packing desired )                   *
!                                                                      *
!----------------------------------------------------------------------*

real(kind=wp) :: PkThrs
integer(kind=iwp) :: Init_do_setup_d, Init_do_setup_e, Init_do_setup_l
logical(kind=iwp) :: isPack

public :: Init_do_setup_d, Init_do_setup_e, Init_do_setup_l, isPack, PkThrs

end module Pack_mod
