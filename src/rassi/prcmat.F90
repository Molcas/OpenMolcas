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

subroutine PRCMAT(NSS,XMATR,XMATI)
! Write out matrix elements over states as a complex matrix
! in square format

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSS
real(kind=wp), intent(in) :: XMATR(NSS,NSS), XMATI(NSS,NSS)
integer(kind=iwp) :: ISS, JEND, JSS, JSTA

do JSTA=1,NSS,2
  JEND = min(NSS,JSTA+1)
  write(u6,*)
  write(u6,'(1X,A8,12X,I3,35X,I3)') ' STATE  ',(JSS,JSS=JSTA,JEND)
  do ISS=1,NSS
    write(u6,'(1X,I4,2x,2(A1,F10.6,A1,F10.6,A1,3x))') ISS,('(',XMATR(ISS,JSS),',',XMATI(ISS,JSS),')',JSS=JSTA,JEND)
  end do
end do

end subroutine PRCMAT
