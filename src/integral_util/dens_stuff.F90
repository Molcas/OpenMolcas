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

module Dens_stuff

! TODO: This should probably be changed to arrays and pointers at some later point.

use Definitions, only: iwp

implicit none
private

integer(kind=iwp), target :: ipDDij, ipDDij2, ipDDik, ipDDik2, ipDDil, ipDDil2, ipDDjk, ipDDjk2, ipDDjl, ipDDjl2, ipDDkl, ipDDkl2, &
                             ipDij, ipDij2, ipDik, ipDik2, ipDil, ipDil2, ipDjk, ipDjk2, ipDjl, ipDjl2, ipDkl, ipDkl2, mDCRij = 1, &
                             mDCRik = 1, mDCRil = 1, mDCRjk = 1, mDCRjl = 1, mDCRkl = 1, mDij, mDik, mDil, mDjk, mDjl, mDkl

public :: ipDDij, ipDDij2, ipDDik, ipDDik2, ipDDil, ipDDil2, ipDDjk, ipDDjk2, ipDDjl, ipDDjl2, ipDDkl, ipDDkl2, ipDij, ipDij2, &
          ipDik, ipDik2, ipDil, ipDil2, ipDjk, ipDjk2, ipDjl, ipDjl2, ipDkl, ipDkl2, mDCRij, mDCRik, mDCRil, mDCRjk, mDCRjl, &
          mDCRkl, mDij, mDik, mDil, mDjk, mDjl, mDkl

end module Dens_stuff
