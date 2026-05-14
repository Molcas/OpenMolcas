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
! Copyright (C) 2012, Victor P. Vysotskiy                              *
!               2025, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  GetMem
!
!> @brief
!>   A simple work space manager, used by all programs in MOLCAS
!> @author Victor P. Vysotskiy
!>
!> @details
!> This is now a simple wrapper around c_getmem, and it should not be used.
!> Use functions from the ::stdalloc module instead.
!>
!> @param[in]     NameIn Arbitrary label
!> @param[in]     KeyIn  ``RGST`` / ``EXCL`` / ``LIST`` / ``TERM``
!> @param[in]     TypeIn ``REAL`` / ``INTE`` / ``CHAR``
!> @param[in,out] iPos   Position
!> @param[in,out] Length Nr of items
!***********************************************************************

subroutine GetMem(NameIn,KeyIn,TypeIn,iPos,Length)
!***********************************************************************
!                                                                      *
! History: Victor P. Vysotskiy                                         *
!    2012: Native Molcas's Memory Allocator; Thread safety             *
!          Ignacio Fdez. Galvan                                        *
!    2025: Garble using C pointers                                     *
!    2025: Reduce to minimum expression                                *
!                                                                      *
!***********************************************************************

use, intrinsic :: iso_c_binding, only: c_null_char
use mma_module, only: c_getmem, MemStat
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: NameIn, KeyIn, TypeIn
integer(kind=iwp), intent(in) :: iPos, Length
#include "warnings.h"
integer(kind=iwp) :: irc, itmp
character(len=9) :: elbl, eopr, etyp

!----------------------------------------------------------------------*
!     Initialize the mma module the first time it is referenced        *
!----------------------------------------------------------------------*
if (.not. MemStat) call IniMem()

!----------------------------------------------------------------------*
!     prepare passed values to the C char format                       *
!----------------------------------------------------------------------*
elbl = NameIn
elbl(9:9) = c_null_char
eopr = KeyIn
eopr(9:9) = c_null_char
etyp = TypeIn
etyp(9:9) = c_null_char

!----------------------------------------------------------------------*
!     Call the C function                                              *
!----------------------------------------------------------------------*
itmp = iPos-1
iRc = c_getmem(elbl,eopr,etyp,itmp,Length)
if (iRc < 0) call Quit(_RC_MEMORY_ERROR_)

end subroutine GetMem
