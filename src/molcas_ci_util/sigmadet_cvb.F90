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

use wadr, only: FMO, TUVX
use CandS, only: ICSM, ISSM
use Lucia_Interface, only: Lucia_Util
use Lucia_Data, only: IREFSM, SIGMA_ON_DISK
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IREFSM_CASVB, NCI
real(kind=wp), intent(_IN_) :: C(NCI)
real(kind=wp), intent(out) :: HC(NCI)

integer(kind=iwp) :: IREFSM_SAVE

IREFSM_SAVE = IREFSM
IREFSM=IREFSM_CASVB

ICSM=IREFSM_CASVB
ISSM=IREFSM_CASVB

! Call the sigma routine
call LUCIA_UTIL('SIGMA', &
                CI_Vector=C, &
                SIGMA_Vector=HC, &
                nTU=size(FMO),TU=FMO, &
                nTUVX=size(TUVX),TUVX=TUVX)

SIGMA_ON_DISK = .True.

call Lucia_Util('Diag', &
                nTU=size(FMO),TU=FMO, &
                nTUVX=size(TUVX),TUVX=TUVX)


IREFSM=IREFSM_SAVE
ICSM=IREFSM
ISSM=IREFSM

end subroutine SIGMADET_CVB
