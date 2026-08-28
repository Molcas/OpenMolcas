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
! Copyright (C) 1991, Jeppe Olsen                                      *
!***********************************************************************

subroutine ADADST(IOBTP,IOBSM,IOBOFF,NIOB,JOBTP,JOBSM,JOBOFF,NJOB,IJORD,ICLS,ISM,IGRP,KMIN,KMAX,I1,XI1S,NK,NKDIM,IEND)
! Obtain mappings
! a+IORB a+ JORB !KSTR> = +/-!ISTR>
! Where IORB belongs to orbitals IOBTP,IOBSM
! and JORB belongs to JOBTP,JOBSM
! In the form
! I1(KSTR) =  ISTR if a+IORB a+ JORB !KSTR> = +/-!ISTR>, ISTR is in
! ICLS,ISM,IGRP.
! (numbering relative to TS start)
!
! Above +/- is stored in XI1S
! Number of K strings checked is returned in NK
! Only Kstrings with relative numbers from KMIN to KMAX are included
!
! If IEND /= 0 last string has been checked
!
! Jeppe Olsen, Winter of 1991

use Str_Info, only: ISTAC, IUNIQMP, NOCTYP, STR
use MCLR_Data, only: NACOB
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IOBTP, IOBSM, IOBOFF, NIOB, JOBTP, JOBSM, JOBOFF, NJOB, IJORD, ICLS, ISM, IGRP, KMIN, KMAX, NKDIM
integer(kind=iwp), intent(_OUT_) :: I1(NKDIM,*)
real(kind=wp), intent(_OUT_) :: XI1S(NKDIM,*)
integer(kind=iwp), intent(out) :: NK, IEND
integer(kind=iwp) :: I1MPF, I2MPF, JGRP, KGRP, L1MP, L2MP

JGRP = IGRP+1
if (IUNIQMP(JGRP) /= JGRP) then
  JGRP = -IUNIQMP(JGRP)
  !write(u6,*) ' Unique string group for mappings ',JGRP
end if
! Are the creation arrays full or in compact form
! N-1 => N
if ((ISTAC(JGRP,1) /= 0) .and. (ISTAC(JGRP,2) /= 0)) then
  I1MPF = 1
  L1MP = NACOB
else
  I1MPF = 0
  L1MP = 0
end if

KGRP = IGRP+2
if (IUNIQMP(KGRP) /= KGRP) then
  KGRP = -IUNIQMP(KGRP)
  !write(u6,*) ' Unique string group for mappings ',KGRP
end if
! N-2 => N-1
if ((ISTAC(KGRP,1) /= 0) .and. (ISTAC(KGRP,2) /= 0)) then
  I2MPF = 1
  L2MP = NACOB
else
  I2MPF = 0
  L2MP = 0
end if

call ADADS1(NK,I1,XI1S,IOBSM,IOBTP,IOBOFF,NIOB,JOBSM,JOBTP,JOBOFF,NJOB,IJORD,NKDIM,ICLS,ISM,Str(KGRP)%STSTM(:,1), &
            Str(KGRP)%STSTM(:,2),I2MPF,L2MP,Str(KGRP)%STSTMI,Str(KGRP)%STSTMN,Str(JGRP)%STSTM(:,1),Str(JGRP)%STSTM(:,2),I1MPF, &
            L1MP,Str(JGRP)%STSTMI,Str(JGRP)%STSTMN,Str(IGRP)%EL1,Str(IGRP)%EL3,Str(IGRP+2)%EL1,Str(IGRP+2)%EL3,Str(IGRP)%ISTSO, &
            Str(IGRP+2)%ISTSO,Str(IGRP+2)%NSTSO,NOCTYP(IGRP),NOCTYP(IGRP+2),NACOB,KMAX,KMIN,IEND)

end subroutine ADADST
