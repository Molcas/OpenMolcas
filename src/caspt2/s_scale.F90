!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine S_SCALE(NAS,SCA,S,iLo,iHi,jLo,jHi,LDS)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NAS, iLo, iHi, jLo, jHi, LDS
real(kind=wp), intent(in) :: SCA(NAS)
real(kind=wp), intent(inout) :: S(LDS,*)
integer(kind=iwp) :: J

do J=jLo,jHi
  S(1:iHi-iLo+1,J-jLo+1) = SCA(iLo:iHi)*SCA(J)*S(1:iHi-iLo+1,J-jLo+1)
end do

end subroutine S_SCALE
