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
! Copyright (C) Jeppe Olsen                                            *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine STRTYP(MS2,NACTEL,MNRS10,MXRS30)
! construct input common blocks /STRINP/
! from /LUCINP/ and /ORBINP/
!
!========
! Output
!========
!   Common block /STRINP/
!
! Jeppe Olsen,  Dec.24, Almaden
!               Last Revision March 31
!
! Where INTXC is internal excitation level 1 => no int exc
!                                          2 => single int exc
!                                          3 => double int exc
! DELTA : number of electrons in string - reference + 5
! ISTTP = 0 => zero order space
! ISTTP = 1 => reference space, no internal excitations
! ISTTP = 2 => reference space, single internal excitations
! ISTTP = 3 => reference space, single internal excitations

use Str_Info, only: ISTAC, IAZTP, IATPM1, IATPM2, IBZTP, IBTPM1, IBTPM2, NSTTYP, NSTTYP_MAX, MNRS1, MXRS1, MNRS3, MXRS3, NELEC, &
                    ISTTP, IZORR, IUNIQMP, IUNIQTP
#ifdef _DEBUGPRINT_
use Str_Info, only: IARTP, IBRTP
#endif
use MCLR_Data, only: NORB1, nORB3

implicit none
integer MS2, NACTEL, MNRS10, MXRS30
! Local variables
logical, external :: Reduce_Prt
integer NAEL, NBEL, IPL, MXRS10, MNRS30, ITYPE, ITYP
integer, external :: iPrintLevel

ISTAC(:,:) = 0
! Number of alpha and beta electrons
NAEL = (MS2+NACTEL)/2
NBEL = (NACTEL-MS2)/2
if (NAEL+NBEL /= NACTEL) then
  write(6,*) 'STRTYP: NAEL + NBEL /= NACTEL'
  write(6,*) 'NAEL,NBEL,NACTEL=',NAEL,NBEL,NACTEL
  !*********************************************************************
  ! The argument iPL was missing so I inserted this piece inside the
  ! stars to calculate it the same way as in the start of mclr.f
  ! //Jonas B
  iPL = iPrintLevel(-1)
  if (Reduce_Prt() .and. (iPL < 3)) iPL = iPL-1
  !*********************************************************************
  call PrInp_MCLR(iPL)
  call Abend()
end if
! Default on alternatice RAS limits
MXRS10 = max(NACTEL,2*NORB1)
MNRS30 = 0
! ===========================
! Strings in zero order space
! ===========================
! Type : alpha-strings
ITYPE = 1
IAZTP = ITYPE
NELEC(ITYPE) = NAEL
MNRS1(ITYPE) = max(0,MNRS10-min(NBEL,NORB1))
MXRS1(ITYPE) = min(NAEL,NORB1,MXRS10)
MNRS3(ITYPE) = max(0,MNRS30-min(NBEL,NORB3))
MXRS3(ITYPE) = min(NAEL,NORB3,MXRS30)

IZORR(ITYPE) = 1
ISTTP(ITYPE) = 0
! Type : single annihilated alphastrings
if (NAEL >= 1) then
  ITYPE = ITYPE+1
  NELEC(ITYPE) = NAEL-1
  MNRS1(ITYPE) = max(0,MNRS1(1)-1)
  MXRS1(ITYPE) = min(NAEL-1,MXRS1(1))
  MNRS3(ITYPE) = max(0,MNRS3(1)-1)
  MXRS3(ITYPE) = min(NAEL-1,MXRS3(1))
  IZORR(ITYPE) = 1
  ISTTP(ITYPE) = 0
  IATPM1 = ITYPE
end if
! Type : double annihilated alphastrings
if (NAEL >= 2) then
  ITYPE = ITYPE+1
  NELEC(ITYPE) = NAEL-2
  MNRS1(ITYPE) = max(0,MNRS1(1)-2)
  MXRS1(ITYPE) = min(NAEL-2,MXRS1(1))
  MNRS3(ITYPE) = max(0,MNRS3(1)-2)
  MXRS3(ITYPE) = min(NAEL-2,MXRS3(1))
  IZORR(ITYPE) = 1
  ISTTP(ITYPE) = 0
  IATPM2 = ITYPE
end if
! Type : beta strings
if (NAEL == NBEL) then
  IBZTP = IAZTP
  IBTPM1 = IATPM1
  IBTPM2 = IATPM2
else
  ITYPE = ITYPE+1
  IBZTP = ITYPE
  NELEC(ITYPE) = NBEL
  MNRS1(ITYPE) = max(0,MNRS10-min(NAEL,NORB1))
  MXRS1(ITYPE) = min(NBEL,NORB1,MXRS10)
  MNRS3(ITYPE) = max(0,MNRS30-min(NAEL,NORB3))
  MXRS3(ITYPE) = min(NBEL,NORB3,MXRS30)
  IZORR(ITYPE) = 1
  ISTTP(ITYPE) = 0
  ! Type : single annihilated betastrings
  if (NBEL >= 1) then
    ITYPE = ITYPE+1
    NELEC(ITYPE) = NBEL-1
    MNRS1(ITYPE) = max(0,MNRS1(IBZTP)-1)
    MXRS1(ITYPE) = min(NBEL-1,MXRS1(IBZTP))
    MNRS3(ITYPE) = max(0,MNRS3(IBZTP)-1)
    MXRS3(ITYPE) = min(NBEL-1,MXRS3(IBZTP))
    IZORR(ITYPE) = 1
    ISTTP(ITYPE) = 0
    IBTPM1 = ITYPE
  end if
  ! Type : double annihilated alphastrings
  if (NBEL >= 2) then
    ITYPE = ITYPE+1
    NELEC(ITYPE) = NBEL-2
    MNRS1(ITYPE) = max(0,MNRS1(IBZTP)-2)
    MXRS1(ITYPE) = min(NBEL-2,MXRS1(IBZTP))
    MNRS3(ITYPE) = max(0,MNRS3(IBZTP)-2)
    MXRS3(ITYPE) = min(NBEL-2,MXRS3(IBZTP))
    IZORR(ITYPE) = 1
    ISTTP(ITYPE) = 0
    IBTPM2 = ITYPE
  end if
end if

NSTTYP = ITYPE
if (NSTTYP > NSTTYP_Max) then
  write(6,*) 'STRTYP: NSTTYP>NSTTYP_Max'
  write(6,*) 'STRTYP: NSTTYP=',NSTTYP
  call Abend()
end if
#ifdef _DEBUGPRINT_
write(6,*) ' Information about string types generated'
write(6,*) ' ========================================'
write(6,*)
write(6,'(A,I3)') ' Number of types generated ',NSTTYP
write(6,*)
write(6,'(A)') ' ==========================================='
write(6,'(A)') '  Type  NELEC MNRS1 MXRS1 MNRS3 MXRS3 ISTTP'
write(6,'(A)') ' ==========================================='
do ITYP=1,NSTTYP
  write(6,'(7I6)') ITYP,NELEC(ITYP),MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP),ISTTP(ITYP)
end do
write(6,*) ' IARTP IBRTP'
call IWRTMA(IARTP,3,7,3,10)
call IWRTMA(IBRTP,3,7,3,10)
#endif

!EAW
do ITYP=1,NSTTYP
  IUNIQMP(ITYP) = ITYP
  IUNIQTP(ITYP) = ITYP
end do
!EAW

end subroutine STRTYP
