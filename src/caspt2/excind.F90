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

subroutine EXCIND(IAS,INS,ISYM,ICASE,IP,IQ,IR,IS)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: IAS, INS, ISYM, ICASE
integer(kind=iwp), intent(out) :: IP, IQ, IR, IS
integer(kind=iwp) :: IR1, IR2

!PAM99 New call sequence for ASIND
call ASIND(IAS,ISYM,ICASE,IR,IP,IR1)
call NSIND(INS,ISYM,ICASE,IS,IQ,IR2)

end subroutine EXCIND
