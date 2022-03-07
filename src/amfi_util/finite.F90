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
! Copyright (C) 1996,1997, Bernd Schimmelpfennig                       *
!***********************************************************************

subroutine finite()
!bs subroutine to set up parameters for finite nucleus. The
!bs s-functions are replaced by just one exponent which models the nucleus.

use AMFI_global, only: charge, Exp_finite, exponents, Lmax_occ, ncontrac, ncontrac_keep, nprimit, noccorb, occup
use Constants, only: Half

implicit none

noccorb(0) = 1
noccorb(1:lmax_occ) = 0
occup(1,0) = -charge
ncontrac_keep = ncontrac(0)
nprimit(0) = 1
ncontrac(0) = 1
exponents(1,0) = Half*Exp_finite

return

end subroutine finite
