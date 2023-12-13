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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

function rMass(nAtom)
!***********************************************************************
!                                                                      *
! Object: to return the mass of the nucleus as a function of the       *
!         atomic number, nAtom. The mass is that one of the most       *
!         abundant isotope. In the case there is not stable isotope    *
!         we select the one with the longest lifetime.                 *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: rMass
integer(kind=iwp), intent(in) :: nAtom
real(kind=wp), external :: rMassx

rMass = rMassx(nAtom,0)

return

end function rMass
