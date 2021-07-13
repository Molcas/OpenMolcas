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

implicit real*8(A-H,O-Z)
#include "ciinfo.fh"
#include "spinfo.fh"
dimension CSFVEC(*), DETVEC(*)
dimension DTOCMT(*), ICTSDT(*)

! To avoid compiler complaints

IOFFCS = 0
IOFFDT = 0
IOFFCD = 0

NTEST = 00

NDET = NDTASM(IREFSM)
NCSF = NCSASM(IREFSM)
if (NTEST >= 100) then
  write(6,*) '======================================='
  write(6,*) '         CSDTVC INFORMATIONS'
  if (IWAY == 1) then
    write(6,*) '   CSF TO DETERMINANT TRANSFORMATION'
  else
    write(6,*) '   DETERMINANT TO CSF TRANSFORMATION'
  end if
  write(6,*) '======================================='
  write(6,*)
  write(6,*) '  NDET                  = ',NDET
  write(6,*) '  NCSF                  = ',NCSF
  write(6,*)
end if

! CSF ==> DET TRANSFORMATION

if (IWAY == 1) then
  if (NTEST >= 100) then
    write(6,*) '   INPUT CSF VECTOR:'
    call WRTMAT(CSFVEC,1,NCSF,1,NCSF)
    write(6,*)
  end if
  call DCOPY_(NDET,[0.0d0],0,DETVEC,1)
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
  if (ICOPY /= 0) call DCOPY_(NDET,DETVEC,1,CSFVEC,1)
  if (NTEST >= 100) then
    write(6,*) '   OUTPUT DET VECTOR:'
    call WRTMAT(DETVEC,1,NDET,1,NDET)
    write(6,*)
  end if

else

  ! DET ==> CSF TRANSFORMATION

  if (NTEST >= 100) then
    write(6,*) '   INPUT DET VECTOR:'
    call WRTMAT(DETVEC,1,NDET,1,NDET)
    write(6,*)
  end if
  call GATVCS(CSFVEC,DETVEC,ICTSDT,NDET)
  if (NTEST >= 100) then
    write(6,*) ' ICTSDT reorder array '
    call IWRTMA(ICTSDT,1,100,1,100)
  end if
  call DCOPY_(NDET,CSFVEC,1,DETVEC,1)
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
  if (ICOPY /= 0) call DCOPY_(NCSF,CSFVEC,1,DETVEC,1)
  if (NTEST >= 100) then
    write(6,*) '   OUTPUT CSF VECTOR:'
    call WRTMAT(CSFVEC,1,NCSF,1,NCSF)
    write(6,*)
  end if
end if

return

end subroutine CSDTVC
