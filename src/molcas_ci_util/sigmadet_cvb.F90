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

subroutine SIGMADET_CVB(C,HC,IREFSM_CASVB,NCI)

use CI_Library, only: CI_Lib_Sigma_Det
use wadr, only: FMO, TUVX
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IREFSM_CASVB, NCI
real(kind=wp), intent(_IN_) :: C(NCI)
real(kind=wp), intent(out) :: HC(NCI)
integer(kind=iwp), parameter :: iState=1

call CI_Lib_Sigma_Det(iState,IREFSM_CASVB,size(FMO),FMO,size(TUVX),TUVX,nCI,C,HC)

end subroutine SIGMADET_CVB
