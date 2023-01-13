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

subroutine CIDIA_CI_UTIL(NCONF,IREFSM,CSFDIA,LUDAVID)
! PURPOSE: - COMPUTE DIAGONAL ELEMENTS OF THE CI-MATRIX
!            THE DETERMINANT BASIS
!          - TRANSLATE FORM DET => CSF BASIS
!
! CALLING PARAMETERS:
! NCONF :  NO. OF CSF
! IREFSM:  REFERENCE SYMMETRY
! CSFDIA:  DIAGONAL OF CI MATRIX IN CSF BASIS

use csfbas, only: CTS
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NCONF, IREFSM, LUDAVID
real(kind=wp), intent(out) :: CSFDIA(NCONF)
integer(kind=iwp) :: iDummy, IPRINT, IPRL, IPRLEV
real(kind=wp) :: Dummy(1), eCore_Hex
real(kind=wp), allocatable :: DDIA(:)
real(kind=wp), external :: Get_eCore
#include "ciinfo.fh"
#include "spinfo.fh"
#include "WrkSpc.fh"
#include "timers.fh"
#include "output_ras.fh"

call Timing(Tissot_1,Swatch,Swatch,Swatch)
IPRLEV = IPRLOC(3)

! COMPUTE CI DIAGONAL IN DETERMINANT BASIS

call Lucia_Util('Diag',iDummy,iDummy,Dummy)

call mma_allocate(DDIA,NDET,label='DETDIA')
call get_diag(DDIA,ndet)

! TRANSFORM CI DIAGONAL FROM DET TO CSF BASIS

IPRINT = 0
if (IPRLEV == INSANE) IPRINT = 40
call CSDIAG_CI_UTIL(NCONF,ndet,CSFDIA,DDIA,NCNFTP(1,IREFSM),NTYP,CTS,NDTFTP,NCSFTP,IPRINT)
eCore_Hex = Get_eCore()
CSFDIA(:) = CSFDIA(:)+eCore_Hex

! DEALLOCATE LOCAL MEMORY

call mma_deallocate(DDIA)

! PRINT CI-DIAGONAL

if (IPRLEV >= DEBUG) then
  IPRL = NCONF
  if (IPRLEV < INSANE) IPRL = min(IPRL,200)
  call dVcPrt('CI-DIAGONAL (max.200 elemwnts)',' ',CSFDIA,NCONF)
end if

! SAVE THE CI_DIAGONAL ON TAPE

call Save_H_diag(nConf,CSFDIA,LUDAVID)

call Timing(Tissot_2,Swatch,Swatch,Swatch)
Tissot_2 = Tissot_2-Tissot_1
Tissot_3 = Tissot_3+Tissot_2

return

end subroutine CIDIA_CI_UTIL
