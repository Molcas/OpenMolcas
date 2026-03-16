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

subroutine ADARASSI(N,A,D,DROT)

use definitions, only: iwp, wp
use constants, only: Zero, One

implicit none
integer(kind=iwp), intent(in) :: N
complex(kind=wp), intent(in) :: A(N,N), D(N,N)
complex(kind=wp), intent(out) :: DROT(N,N)
integer(kind=iwp) I, J
complex(kind=wp) TEMP(N,N)

! initialization
do I=1,N
  do J=1,N
    DROT(I,J) = (Zero,Zero)
    TEMP(I,J) = (Zero,Zero)
  end do
end do

! actual multiplication
call ZGEMM('C','N',N,N,N,(One,Zero),A,N,D,N,(Zero,Zero),TEMP,N)
call ZGEMM('N','N',N,N,N,(One,Zero),TEMP,N,A,N,(Zero,Zero),DROT,N)

end subroutine ADARASSI
