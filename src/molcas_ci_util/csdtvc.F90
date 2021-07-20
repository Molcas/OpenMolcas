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

subroutine CSDTVC(CSFVEC,DETVEC,IWAY,DTOCMT,ICTSDT,IREFSM,ICOPY)
! PURPOSE: TRANSFORM FROM DETERMINANT TO CSF BASIS AND VICE VERSA
!          IWAY = 1 : CSF TO DETERMINANT TRANSFORMATION
!          IWAY = 2 : DETERMINANT TO CSF TRANSFORMATION

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "ciinfo.fh"
#include "spinfo.fh"
integer(kind=iwp), intent(in) :: IWAY, ICTSDT(*), IREFSM, ICOPY
real(kind=wp), intent(inout) :: CSFVEC(NDTASM(IREFSM)), DETVEC(NDTASM(IREFSM))
real(kind=wp), intent(in) :: DTOCMT(*)
integer(kind=iwp) :: ICNF, ICSF, IDET, IOFFCD, IOFFCS, IOFFDT, ITYP, NCSF, NTEST

! To avoid compiler complaints

IOFFCS = 0
IOFFDT = 0
IOFFCD = 0

NTEST = 00

NDET = NDTASM(IREFSM)
NCSF = NCSASM(IREFSM)
if (NTEST >= 100) then
  write(u6,*) '======================================='
  write(u6,*) '         CSDTVC INFORMATIONS'
  if (IWAY == 1) then
    write(u6,*) '   CSF TO DETERMINANT TRANSFORMATION'
  else
    write(u6,*) '   DETERMINANT TO CSF TRANSFORMATION'
  end if
  write(u6,*) '======================================='
  write(u6,*)
  write(u6,*) '  NDET                  = ',NDET
  write(u6,*) '  NCSF                  = ',NCSF
  write(u6,*)
end if

! CSF ==> DET TRANSFORMATION

if (IWAY == 1) then
  if (NTEST >= 100) then
    write(u6,*) '   INPUT CSF VECTOR:'
    call WRTMAT(CSFVEC,1,NCSF,1,NCSF)
    write(u6,*)
  end if
  DETVEC(:) = Zero
  do ITYP=1,NTYP
    IDET = NDTFTP(ITYP)
    ICSF = NCSFTP(ITYP)
    ICNF = NCNFTP(ITYP,IREFSM)
    if (ITYP == 1) then
      IOFFCS = 1
      IOFFDT = 1
      IOFFCD = 1
    else
      IOFFCS = IOFFCS+NCNFTP(ITYP-1,IREFSM)*NCSFTP(ITYP-1)
      IOFFDT = IOFFDT+NCNFTP(ITYP-1,IREFSM)*NDTFTP(ITYP-1)
      IOFFCD = IOFFCD+NDTFTP(ITYP-1)*NCSFTP(ITYP-1)
    end if
    if ((IDET*ICNF*ICSF) > 0) call MATML4(DETVEC(IOFFDT),DTOCMT(IOFFCD),CSFVEC(IOFFCS),IDET,ICNF,IDET,ICSF,ICSF,ICNF,0)
  end do
  call Sort_Cdet(nDet,ICTSDT,DetVec)
  if (ICOPY /= 0) CSFVEC(:) = DETVEC(:)
  if (NTEST >= 100) then
    write(u6,*) '   OUTPUT DET VECTOR:'
    call WRTMAT(DETVEC,1,NDET,1,NDET)
    write(u6,*)
  end if

else

  ! DET ==> CSF TRANSFORMATION

  if (NTEST >= 100) then
    write(u6,*) '   INPUT DET VECTOR:'
    call WRTMAT(DETVEC,1,NDET,1,NDET)
    write(u6,*)
  end if
  call GATVCS(CSFVEC,DETVEC,ICTSDT,NDET)
  if (NTEST >= 100) then
    write(u6,*) ' ICTSDT reorder array'
    call IWRTMA(ICTSDT,1,100,1,100)
  end if
  DETVEC(:) = CSFVEC(:)
  do ITYP=1,NTYP
    IDET = NDTFTP(ITYP)
    ICSF = NCSFTP(ITYP)
    ICNF = NCNFTP(ITYP,IREFSM)
    if (ITYP == 1) then
      IOFFCS = 1
      IOFFDT = 1
      IOFFCD = 1
    else
      IOFFCS = IOFFCS+NCNFTP(ITYP-1,IREFSM)*NCSFTP(ITYP-1)
      IOFFDT = IOFFDT+NCNFTP(ITYP-1,IREFSM)*NDTFTP(ITYP-1)
      IOFFCD = IOFFCD+NDTFTP(ITYP-1)*NCSFTP(ITYP-1)
    end if
    if ((IDET*ICNF*ICSF) > 0) call MATML4(CSFVEC(IOFFCS),DTOCMT(IOFFCD),DETVEC(IOFFDT),ICSF,ICNF,IDET,ICSF,IDET,ICNF,1)
  end do
  if (ICOPY /= 0) DETVEC(1:NCSF) = CSFVEC(1:NCSF)
  if (NTEST >= 100) then
    write(u6,*) '   OUTPUT CSF VECTOR:'
    call WRTMAT(CSFVEC,1,NCSF,1,NCSF)
    write(u6,*)
  end if
end if

return

end subroutine CSDTVC
