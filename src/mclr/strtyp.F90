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

use Str_Info, only: IATPM1, IATPM2, IAZTP, IBTPM1, IBTPM2, IBZTP, ISTAC, IUNIQMP, IUNIQTP, MNRS1, MNRS3, MXRS1, MXRS3, NELEC, NSTTYP
use MCLR_Data, only: NORB1, NORB3
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: MS2, NACTEL, MNRS10, MXRS30
integer(kind=iwp) :: IPL, ITYP, ITYPE, MNRS30, MXRS10, NAEL, NBEL, NSTTYP_Max
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt

ISTAC(:,:) = 0
! Number of alpha and beta electrons
NAEL = (MS2+NACTEL)/2
NBEL = (NACTEL-MS2)/2
if (NAEL+NBEL /= NACTEL) then
  write(u6,*) 'STRTYP: NAEL + NBEL /= NACTEL'
  write(u6,*) 'NAEL,NBEL,NACTEL=',NAEL,NBEL,NACTEL
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

! Type : single annihilated alphastrings
if (NAEL >= 1) then
  ITYPE = ITYPE+1
  NELEC(ITYPE) = NAEL-1
  MNRS1(ITYPE) = max(0,MNRS1(1)-1)
  MXRS1(ITYPE) = min(NAEL-1,MXRS1(1))
  MNRS3(ITYPE) = max(0,MNRS3(1)-1)
  MXRS3(ITYPE) = min(NAEL-1,MXRS3(1))
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
  ! Type : single annihilated betastrings
  if (NBEL >= 1) then
    ITYPE = ITYPE+1
    NELEC(ITYPE) = NBEL-1
    MNRS1(ITYPE) = max(0,MNRS1(IBZTP)-1)
    MXRS1(ITYPE) = min(NBEL-1,MXRS1(IBZTP))
    MNRS3(ITYPE) = max(0,MNRS3(IBZTP)-1)
    MXRS3(ITYPE) = min(NBEL-1,MXRS3(IBZTP))
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
    IBTPM2 = ITYPE
  end if
end if

NSTTYP_Max = size(NELEC)
NSTTYP = ITYPE
if (NSTTYP > NSTTYP_Max) then
  write(u6,*) 'STRTYP: NSTTYP>NSTTYP_Max'
  write(u6,*) 'STRTYP: NSTTYP=',NSTTYP
  call Abend()
end if
#ifdef _DEBUGPRINT_
write(u6,*) ' Information about string types generated'
write(u6,*) ' ========================================'
write(u6,*)
write(u6,'(A,I3)') ' Number of types generated ',NSTTYP
write(u6,*)
write(u6,'(A)') ' ====================================='
write(u6,'(A)') '  Type  NELEC MNRS1 MXRS1 MNRS3 MXRS3'
write(u6,'(A)') ' ====================================='
do ITYP=1,NSTTYP
  write(u6,'(6I6)') ITYP,NELEC(ITYP),MNRS1(ITYP),MXRS1(ITYP),MNRS3(ITYP),MXRS3(ITYP)
end do
#endif

!EAW
IUNIQMP(1:NSTTYP) = [(ITYP,ITYP=1,NSTTYP)]
IUNIQTP(1:NSTTYP) = [(ITYP,ITYP=1,NSTTYP)]
!EAW

end subroutine STRTYP
