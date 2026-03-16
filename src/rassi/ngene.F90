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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************

function NGENE(NEL,MLTPL)
! Nr of genealogical spin couplings

use definitions, only: iwp

implicit none
integer(kind=iwp) NGENE
integer(kind=iwp), intent(in) :: NEL, MLTPL
integer(kind=iwp) IS2, NU, ND
integer(kind=iwp), external :: NOVERM

NGENE = 0
if (MLTPL <= 0) return
IS2 = MLTPL-1
if (NEL < IS2) return
NU = (NEL+IS2)/2
ND = (NEL-IS2)/2
if (NU+ND /= NEL) return
NGENE = NOVERM(NEL,NU)-NOVERM(NEL,NU+1)

end function NGENE
