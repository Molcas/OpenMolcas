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
! Copyright (C) Yingjin Ma                                             *
!***********************************************************************

subroutine CSDTVC_dmrg(CSFVEC,DETVEC,DTOCMT,ICTSDT,IREFSM,ICOPY)
! DETERMINANT TO CSF TRANSFORMATION
!
! ICOPY /= 0 : Copy output into input
!              so input becomes output while
!              output remains output
! Modified version for DMRG only -- yma

use MCLR_Data, only: NCNATS, NCPCNT, NCSASM, NDPCNT, NDTASM, NTYP
use Constants, only: Zero, One
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: CSFVEC(*)
real(kind=wp), intent(inout) :: DETVEC(*)
real(kind=wp), intent(in) :: DTOCMT(*)
integer(kind=iwp), intent(in) :: ICTSDT(*), IREFSM, ICOPY
integer(kind=iwp) :: ICNF, ICSF, IDET, IOFFCD, IOFFCS, IOFFDT, ITYP, NCSF, NDET

NDET = NDTASM(IREFSM)
NCSF = NCSASM(IREFSM)

! ===================================
!  Determinant to csf transformation
! ===================================

! To CSF ordering

call GATVCS(CSFVEC,DETVEC,ICTSDT,NDET)
DETVEC(1:NDET) = CSFVEC(1:NDET)
! Multiply with CIND expansion matrix
IOFFCS = 1
IOFFDT = 1
IOFFCD = 1
do ITYP=1,NTYP
  IDET = NDPCNT(ITYP)
  ICSF = NCPCNT(ITYP)
  ICNF = NCNATS(ITYP,IREFSM)
  if (ITYP > 1) then
    IOFFCS = IOFFCS+NCNATS(ITYP-1,IREFSM)*NCPCNT(ITYP-1)
    IOFFDT = IOFFDT+NCNATS(ITYP-1,IREFSM)*NDPCNT(ITYP-1)
    IOFFCD = IOFFCD+NDPCNT(ITYP-1)*NCPCNT(ITYP-1)
  end if
  if (IDET*ICNF*ICSF > 0) call DGEMM_('T','N',ICSF,ICNF,IDET,ONE,DTOCMT(IOFFCD),IDET,DETVEC(IOFFDT),IDET,ZERO,CSFVEC(IOFFCS),ICSF)
end do
if (ICOPY /= 0) DETVEC(1:NCSF) = CSFVEC(1:NCSF)

end subroutine CSDTVC_dmrg
