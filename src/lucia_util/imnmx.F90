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

function IMNMX(IVEC,NDIM,MINMAX)
! Find smallest (MINMAX=1) or largest (MINMAX=2)
! absolute value of elements in integer vector IVEC

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: IMNMX
integer(kind=iwp), intent(in) :: NDIM, IVEC(NDIM), MINMAX
integer(kind=iwp) :: IX

IX = 0
if (NDIM > 0) then
  IX = -1
  if (MINMAX == 1) IX = minval(abs(IVEC(:)))

  if (MINMAX == 2) IX = maxval(abs(IVEC(:)))

else
  ! No components : set to zero and write a warning
  IX = 0
  write(u6,*) ' Min/Max taken zero length vector set to zero'
end if

IMNMX = IX

end function IMNMX
