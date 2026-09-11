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

use ci_interfaces, only: Mk_H_Psi_Det
use wadr, only: FMO, TUVX
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IREFSM_CASVB, NCI
real(kind=wp), intent(_IN_) :: C(NCI)
real(kind=wp), intent(out) :: HC(NCI)

Call Mk_H_Psi_Det(IREFSM_CASVB,nCI,C,HC,nTU=size(FMO),TU=FMO,nTUVX=size(TUVX),TUVX=TUVX)

end subroutine SIGMADET_CVB
