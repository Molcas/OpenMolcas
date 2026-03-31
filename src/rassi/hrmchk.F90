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

subroutine HRMCHK(NDIMEN,ARRRE,ARRIM,ERRRE,ERRIM)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NDIMEN
real(kind=wp) :: ARRRE(NDIMEN,NDIMEN), ARRIM(NDIMEN,NDIMEN), ERRRE, ERRIM
integer(kind=iwp) :: I, J

ERRRE = Zero
ERRIM = Zero
do I=1,NDIMEN
  do J=1,I-1
    ERRRE = max(ERRRE,abs(ARRRE(I,J)-ARRRE(J,I)))
    ERRIM = max(ERRIM,abs(ARRIM(I,J)+ARRIM(J,I)))
  end do
end do

end subroutine HRMCHK
