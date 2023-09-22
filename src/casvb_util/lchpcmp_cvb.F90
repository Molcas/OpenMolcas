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

!IFG trivial
function lchpcmp_cvb(ltst)

use Definitions, only: iwp

implicit none
logical(kind=iwp) :: lchpcmp_cvb
logical(kind=iwp) :: ltst
integer(kind=iwp) :: itst
logical(kind=iwp), external :: chpcmp_cvb

if (ltst) then
  itst = 1
else
  itst = 0
end if
lchpcmp_cvb = chpcmp_cvb(itst)

return

end function lchpcmp_cvb
