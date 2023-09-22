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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************
!*******************
!** File handling **
!*******************

!IFG trivial
subroutine lendat_cvb(fileid,lenf)
! return in lenf the length of record name in file ifil
!   lenf = -1 if record does not exist

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: fileid
integer(kind=iwp) :: lenf

lenf = -1

return
! Avoid unused argument warnings
if (.false.) call Unused_real(fileid)

end subroutine lendat_cvb
