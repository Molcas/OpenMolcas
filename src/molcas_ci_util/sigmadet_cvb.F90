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

subroutine SIGMADET_CVB(C,HC,IREFSM,NCI)

use Definitions, only: wp, iwp
use Lucia_Interface, only: Lucia_Util

implicit none
integer(kind=iwp), intent(in) :: IREFSM, NCI
real(kind=wp), intent(in) :: C(NCI)
real(kind=wp), intent(out) :: HC(NCI)

! Export arguments to be used in sigma_master_cvb

! Call the sigma routine
call LUCIA_UTIL('SIGMA_CVB',         &
                CI_Vector=C, &
                SIGMA_Vector=HC, &
                iSym=iREFSM)

end subroutine SIGMADET_CVB
