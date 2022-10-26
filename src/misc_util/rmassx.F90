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

function rMassx(nAtom,nIso)
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

use Isotopes
implicit real*8(A-H,O-Z)
#include "real.fh"
#include "constants2.fh"
real*8 rMassx
integer nIso

if (nAtom > MaxAtomNum) then
  !write(6,*) ' Weight for this atom is not listed!'
  !write(6,*) ' Mass set to 2.6 times atom number'
  rMassx = 2.6d0*dble(nAtom)*uToau
else if (nAtom == 0) then
  !write(6,*) ' Weight for this atom is meaningless!'
  !write(6,*) ' Mass set to 0.0'
  rMassx = Zero
else if (nAtom < 0) then
  rMassx = 1.0d99*uToau
else
  isnx = nIso
  call Isotope(isnx,nAtom,rMassx)
end if

return

end function rMassx
