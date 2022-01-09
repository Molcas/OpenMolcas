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

integer function Index_Center(iCnt,iR,Index,iAtoms,nAtoms)

integer index(2,nAtoms)

98 continue
Index_Center = 0
do i=1,iAtoms
  if ((index(1,i) == iCnt) .and. (index(2,i) == iR)) then
    Index_Center = i
    Go To 99
  end if
end do

iAtoms = iAtoms+1
index(1,iAtoms) = iCnt
index(2,iAtoms) = iR
Go To 98

99 return

end function Index_Center
