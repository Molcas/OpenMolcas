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

!#define _DEBUGPRINT_
subroutine CIDIA(NCONF,CSFDIA,LUDAVID)
! PURPOSE: - COMPUTE DIAGONAL ELEMENTS OF THE CI-MATRIX
!            THE DETERMINANT BASIS
!          - TRANSLATE FORM DET => CSF BASIS
!
! CALLING PARAMETERS:
! NCONF :  NO. OF CSF
! CSFDIA:  DIAGONAL OF CI MATRIX IN CSF BASIS

use CI_Library, only: CI_Diag_Cnf_Avg, CI_Lib_Diag
use CI_Info, only: CI_St
use wadr, only: FMO, TUVX
use rasscf_global, only: NSEL
use timers, only: TimeHDiag
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NCONF, LUDAVID
real(kind=wp), intent(out) :: CSFDIA(NCONF)
real(kind=wp) :: dum1, dum2, dum3, Time(2)
integer(kind=iwp), parameter :: iState = 1

call Timing(Time(1),dum1,dum2,dum3)

! the configuration average of the CI library (the CI state of RASSCF; the
! P space of NSEL CSFs for the checks of MOLCAS_CI_VERIFY)
call CI_Lib_Diag(iState,CI_Diag_Cnf_Avg,CI_St(iState)%STSym,size(FMO),FMO,size(TUVX),TUVX,Zero,nConf,CSFDIA,nPMax=NSEL)

! PRINT CI-DIAGONAL

#ifdef _DEBUGPRINT_
call RecPrt('CI-DIAGONAL (max.200 elements)',' ',CSFDIA,1,nConf)
#endif

! SAVE THE CI_DIAGONAL ON TAPE

call Save_H_diag(nConf,CSFDIA,LUDAVID)

call Timing(Time(2),dum1,dum2,dum3)
TimeHDiag = TimeHDiag+Time(2)-Time(1)

end subroutine CIDIA
