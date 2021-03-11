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
! Copyright (C) 2008, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine initialize the peek/poke utility.                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: May 2008                                                    *
!                                                                      *
!***********************************************************************
!  Init_ppu
!
!> @brief
!>   Initialize the peek/poke utility
!> @author Per-Olof Widmark
!>
!> @details
!> This routine is used to initialize the peek/poke utility.
!>
!> @param[in] Force Force initialization
!***********************************************************************

subroutine Init_ppu(Force)

use peekpoke, only: ds_no, is_no
use Definitions, only: iwp

implicit none
!----------------------------------------------------------------------*
! Dummy arguments.                                                     *
!----------------------------------------------------------------------*
logical(kind=iwp), intent(in) :: Force
!----------------------------------------------------------------------*
! Local variables.                                                     *
!----------------------------------------------------------------------*
logical(kind=iwp) :: FirstTime = .true.

!----------------------------------------------------------------------*
! Initialize the peek/poke utility                                     *
!----------------------------------------------------------------------*
if (FirstTime .or. Force) then
  ds_no = 0
  is_no = 0
end if
FirstTime = .false.
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
return

end subroutine Init_ppu
