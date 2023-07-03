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

integer function iUR(iR,iU)

iUR = 0
do i=0,7
  if (iand(iU,2**i) == 2**i) then
    iUR = ior(iUR,2**ieor(i,iR))
  end if
end do
!write(6,*) ' iUR=',iUR

return

end function iUR
