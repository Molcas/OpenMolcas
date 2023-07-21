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

subroutine CHO_INIMAP()
!
! Purpose: initialize integral shell pair calculation map.
!
! NB!!!!! file is assumed open (restart only)

use ChoArr, only: IntMap

implicit none
#include "cholesky.fh"
integer IOPT, NDIM, IADR

if (RSTCHO) then
  IOPT = 2
  NDIM = size(INTMAP)
  IADR = 0
  call IDAFILE(LUMAP,IOPT,INTMAP,NDIM,IADR)
else
  call IZERO(INTMAP,size(INTMAP))
end if

end subroutine CHO_INIMAP
