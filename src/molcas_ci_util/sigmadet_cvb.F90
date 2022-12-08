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

implicit none
integer(kind=iwp), intent(in) :: IREFSM, NCI
real(kind=wp), intent(in) :: C(NCI)
real(kind=wp), intent(out) :: HC(NCI)
integer(kind=iwp) :: IDUMMY
real(kind=wp) :: DUMMY(1)
#include "WrkSpc.fh"
#include "rasscf_lucia.fh"

! Export arguments to be used in sigma_master_cvb

C_POINTER = KCI_POINTER
call DCOPY_(NCI,C,1,WORK(C_POINTER),1)
! Call the sigma routine
call LUCIA_UTIL('SIGMA_CVB',IREFSM,IDUMMY,DUMMY)
call DCOPY_(NCI,WORK(KSIGMA_POINTER),1,HC,1)

return

end subroutine SIGMADET_CVB
