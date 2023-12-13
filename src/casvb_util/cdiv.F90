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

subroutine cdiv(ar,ai,br,bi,cr,ci)
! complex division, (cr,ci) = (ar,ai)/(br,bi)

use Definitions, only: wp

implicit none
real(kind=wp), intent(in) :: ar, ai, br, bi
real(kind=wp), intent(out) :: cr, ci
real(kind=wp) :: ais, ars, bis, brs, s

s = abs(br)+abs(bi)
ars = ar/s
ais = ai/s
brs = br/s
bis = bi/s
s = brs**2+bis**2
cr = (ars*brs+ais*bis)/s
ci = (ais*brs-ars*bis)/s

return

end subroutine cdiv
