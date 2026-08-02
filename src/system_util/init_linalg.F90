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
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Init_LinAlg()
! Initialization procedure for linear algebra libraries (if needed)

#ifdef _DELAYED_

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iPr
character(len=1024) :: linalg_lib
character(len=8) :: linalg_info

linalg_lib = 'Internal'
call GetEnvf('MOLCAS_LINALG',linalg_lib)
call GetEnvf('MOLCAS_LINALG_INFO',linalg_info)
iPr = 0
call UpCase(linalg_info)
linalg_info = adjustl(linalg_info)
if ((linalg_info(1:1) /= ' ') .and. &
    (linalg_info(1:3) /= 'NO ') .and. &
    (linalg_info(1:4) /= 'OFF ') .and. &
    (linalg_info(1:1) /= '0')) iPr = 1
call Initialize_BLAS(linalg_lib,iPr)

#endif

end subroutine Init_LinAlg
