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

subroutine PRTSTR(ISTR,NEL,NSTR)
! Print NSTR strings each containing NEL electrons

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NEL, NSTR, ISTR(NEL,NSTR)
integer(kind=iwp) :: IEL, JSTR

do JSTR=1,NSTR
  write(u6,100) ' String ',JSTR,' : ',(ISTR(IEL,JSTR),IEL=1,NEL)
end do

return
100 format('0',A,I6,A,4X,10(2X,I4),/,(1X,19X,10(2X,I4)))

end subroutine PRTSTR
