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

subroutine CSDTVC_dmrg(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,IREFSM,ICOPY)
! IWAY = 1 : CSF to DETERMINANT TRANSFORMATION
! IWAY = 2 : DETERMINANT TO CSF TRANSFORMATION
!
! ICOPY /= 0 : Copy output into input
!              so input becomes output while
!              output remains output
! Modified version for DMRG only -- yma

use Constants, only: Zero, One
use MCLR_Data, only: NTYP, NCNATS, NCPCNT, NCSASM, NDPCNT, NDTASM

implicit none
real*8 CSFVEC(*), DETVEC(*)
integer IWAY
real*8 DTOCMT(*)
integer ICTSDT(*)
integer IREFSM, ICOPY
integer IOFFCS, IOFFDT, IOFFCD, NDET, NCSF, ITYP, IDET, ICSF, ICNF

IOFFCS = 0 ! dummy initialize
IOFFDT = 0 ! dummy initialize
IOFFCD = 0 ! dummy initialize

NDET = NDTASM(IREFSM)
NCSF = NCSASM(IREFSM)

if (IWAY == 1) then

  ! ===========================
  !  CSF to DET transformation
  ! ===========================

  DETVEC(1:NDET) = Zero
  ! Multiply with  expansion matrix
  do ITYP=1,NTYP
    IDET = NDPCNT(ITYP)
    ICSF = NCPCNT(ITYP)
    ICNF = NCNATS(ITYP,IREFSM)
    if (ITYP == 1) then
      IOFFCS = 1
      IOFFDT = 1
      IOFFCD = 1
    else
      IOFFCS = IOFFCS+NCNATS(ITYP-1,IREFSM)*NCPCNT(ITYP-1)
      IOFFDT = IOFFDT+NCNATS(ITYP-1,IREFSM)*NDPCNT(ITYP-1)
      IOFFCD = IOFFCD+NDPCNT(ITYP-1)*NCPCNT(ITYP-1)
    end if
    if (IDET*ICNF*ICSF > 0) call DGEMM_('N','N',IDET,ICNF,ICSF,ONE,DTOCMT(IOFFCD),IDET,CSFVEC(IOFFCS),ICSF,ZERO,DETVEC(IOFFDT),IDET)
  end do
  ! Sign changes
  CSFVEC(1:NDET) = DETVEC(1:NDET)
  ! Change to string ordering
  call SCAVCS(DETVEC,CSFVEC,ICTSDT,NDET)
  if (ICOPY /= 0) CSFVEC(1:NDET) = DETVEC(1:NDET)
else

  ! ===================================
  !  Determinant to csf transformation
  ! ===================================

  ! To CSF ordering

  call GATVCS(CSFVEC,DETVEC,ICTSDT,NDET)
  DETVEC(1:NDET) = CSFVEC(1:NDET)
  ! Multiply with CIND expansion matrix
  do ITYP=1,NTYP
    IDET = NDPCNT(ITYP)
    ICSF = NCPCNT(ITYP)
    ICNF = NCNATS(ITYP,IREFSM)
    if (ITYP == 1) then
      IOFFCS = 1
      IOFFDT = 1
      IOFFCD = 1
    else
      IOFFCS = IOFFCS+NCNATS(ITYP-1,IREFSM)*NCPCNT(ITYP-1)
      IOFFDT = IOFFDT+NCNATS(ITYP-1,IREFSM)*NDPCNT(ITYP-1)
      IOFFCD = IOFFCD+NDPCNT(ITYP-1)*NCPCNT(ITYP-1)
    end if
    if (IDET*ICNF*ICSF > 0) call DGEMM_('T','N',ICSF,ICNF,IDET,ONE,DTOCMT(IOFFCD),IDET,DETVEC(IOFFDT),IDET,ZERO,CSFVEC(IOFFCS),ICSF)
  end do
  if (ICOPY /= 0) DETVEC(1:NCSF) = CSFVEC(1:NCSF)
end if

end subroutine CSDTVC_dmrg
