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

subroutine CIDIA_CI_UTIL(NORB,NCONF,IREFSM,CSFDIA,G,LUDAVID)
! PURPOSE: - COMPUTE DIAGONAL ELEMENTS OF THE CI-MATRIX
!            THE DETERMINANT BASIS
!          - TRANSLATE FORM DET => CSF BASIS
!
! CALLING PARAMETERS:
! NORB  :  NO. OF ACTIVE ORBITALS
! NCONF :  NO. OF CSF
! IREFSM:  REFERENCE SYMMETRY
! CSFDIA:  DIAGONAL OF CI MATRIX IN CSF BASIS
! G     :  MODIFIED ONE ELECTRON HAMILTONIAN INCLUDING CORE ELECTR.

use Constants, only: One
use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp), intent(in) :: NORB, NCONF, IREFSM, LUDAVID
real(kind=wp), intent(out) :: CSFDIA(*)
real(kind=wp), intent(in) :: G(*)
integer(kind=iwp) :: I, iDummy, II, IPRINT, IPRL, IPRLEV, KDDIA, KH1DIA, KSCR, KX
real(kind=wp) :: Dummy(1), eCore_Hex
real(kind=r8), external :: Get_eCore
#include "ciinfo.fh"
#include "spinfo.fh"
#include "csfbas.fh"
#include "WrkSpc.fh"
#include "timers.fh"
#include "output_ras.fh"

call Timing(Tissot_1,Swatch,Swatch,Swatch)
IPRLEV = IPRLOC(3)

! ALLOCATE LOCAL MEMORY

call GETMEM('XA','ALLO','REAL',KX,NORB)
call GETMEM('SCR','ALLO','REAL',KSCR,2*NORB)
call GETMEM('H1DIA','ALLO','REAL',KH1DIA,NORB)

! SELECT DIAGONAL ONEBODY INTEGRALS

II = 0
do I=1,NORB
  II = II+I
  Work(KH1DIA-1+I) = G(II)
end do

! COMPUTE CI DIAGONAL IN DETERMINANT BASIS

call Lucia_Util('Diag',iDummy,iDummy,Dummy)

call GETMEM('DETDIA','ALLO','REAL',KDDIA,NDET)
call get_diag(work(kddia),ndet)

! TRANSFORM CI DIAGONAL FROM DET TO CSF BASIS

IPRINT = 0
if (IPRLEV == INSANE) IPRINT = 40
call CSDIAG_CI_UTIL(CSFDIA,Work(KDDIA),NCNFTP(1,IREFSM),NTYP,iWork(KICTS(1)),NDTFTP,NCSFTP,IPRINT)
eCore_Hex = Get_eCore()
call DAXPY_(NCONF,eCore_Hex,[One],0,CSFDIA,1)

! DEALLOCATE LOCAL MEMORY

call GETMEM('XA','FREE','REAL',KX,NORB)
call GETMEM('SCR','FREE','REAL',KSCR,2*NORB)
call GETMEM('H1DIA','FREE','REAL',KH1DIA,NORB)
call GETMEM('DETDIA','FREE','REAL',KDDIA,NDET)

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
