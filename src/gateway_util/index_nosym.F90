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

integer function Index_NoSym(iCntr,iCmp,iCnt,iAng,iR,Index,iBas,nBas)

integer index(5,nBas)

98 continue
Index_NoSym = 0
do i=1,iBas
  if ((index(1,i) == iCntr) .and. (index(2,i) == iCmp) .and. (index(3,i) == iCnt) .and. (index(4,i) == iAng) .and. &
      (index(5,i) == iR)) then
    Index_NoSym = i
    Go To 99
  end if
end do

iBas = iBas+1
index(1,iBas) = iCntr
index(2,iBas) = iCmp
index(3,iBas) = iCnt
index(4,iBas) = iAng
index(5,iBas) = iR
Go To 98

99 return

end function Index_NoSym
