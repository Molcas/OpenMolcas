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

integer function norder(icoord,int_code,lmax)

implicit none
integer lmax
integer icoord(lmax), int_code(lmax)
integer nb, isite

nb = 0
do isite=1,lmax
  nb = nb+icoord(isite)*int_code(isite)
end do

norder = nb+1

return

end function norder
