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

module welcom

use Constants, only: Two, Four
use define_af, only: iTabMx

implicit none
private
integer, parameter :: kMax = iTabMx+6
integer ipot3(0:kmax+1)
real*8 binom(-1:kmax,-1:kmax), fiint(0:kmax,0:kmax), tetint(0:kmax,0:int(kmax/Two)+1), &
       anorm(0:kmax,0:int(kmax/Two)+1,0:int(kmax/Four)+1), fac(0:kmax)

public kMax, ipot3, binom, fiint, tetint, anorm, fac

end module welcom
