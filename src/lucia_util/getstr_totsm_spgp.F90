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
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

subroutine GETSTR_TOTSM_SPGP(ISTRTP,ISPGRP,ISPGRPSM,NEL,NSTR,ISTR,NORBT,IDOREO,IZ,IREO)
! Obtain all super-strings of given total symmetry and given
! occupation in each GAS space
!
! If IDOREO /= 0 THEN reordering array : lexical => actual order is obtained
!
! Nomenclature of the day : superstring : string in complete
!                           orbital space, product of strings in
!                           each GAS space
! =====
! Input
! =====
!
! ISTRTP  : Type of of superstrings ( alpha => 1, beta => 2 )
! ISPGRP :  supergroup number, (relative to start of this type )
! ISPGRPSM : Total symmetry of superstrings
! NEL : Number of electrons
! IZ  : Reverse lexical ordering array for this supergroup
!
! ======
! Output
! ======
!
! NSTR : Number of superstrings generated
! ISTR : Occupation of superstring
! IREO : Reorder array (if IDOREO /= 0)
!
! Jeppe Olsen, July 1995

use strbas, only: NSTSGP, ISTSGP
use lucia_data, only: NGAS
use lucia_data, only: IBSPGPFTP, ISPGPFTP, NELFGP
use lucia_data, only: MXPNGAS, MXPNSMST
use csm_data, only: NSMST
use Definitions, only: u6

implicit none
integer ISTRTP, ISPGRP, ISPGRPSM, NEL, NSTR, NORBT, IDOREO
integer IZ(NORBT,NEL)
! output
integer ISTR(*), IREO(*)
! Local scratch
integer NELFGS(MXPNGAS), ISMFGS(MXPNGAS), ITPFGS(MXPNGAS)
integer maxval(MXPNGAS), minval(MXPNGAS)
integer NNSTSGP(MXPNSMST,MXPNGAS)
integer IISTSGP(MXPNSMST,MXPNGAS)
integer NTEST, ISPGRPA, NGASL, IGAS, ISMST, MAXLEX, IFIRST, ISTRBS, NONEW, ISTSMM1, JSTSMM1, ISMGSN, JSTR, LEX, IEL

NTEST = 0
if (NTEST >= 100) then
  write(u6,*)
  write(u6,*) ' ============================'
  write(u6,*) ' Welcome to GETSTR_TOTSM_SPGP'
  write(u6,*) ' ============================'
  write(u6,*)
  write(u6,'(A,3I3)') ' Strings to be obtained : Type, supergroup, symmetry ',ISTRTP,ISPGRP,ISPGRPSM
  write(u6,*)
end if
! Absolute number of this supergroup
ISPGRPA = IBSPGPFTP(ISTRTP)-1+ISPGRP
! Occupation per gasspace
! Largest occupied space
NGASL = 0
! Largest and lowest symmetries active in each GAS space
do IGAS=1,NGAS
  ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPA)
  NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
  if (NELFGS(IGAS) > 0) NGASL = IGAS
end do
if (NGASL == 0) NGASL = 1
! Number of strings per GAS space and offsets for strings of given sym
do IGAS=1,NGAS
  call ICOPVE2(NSTSGP(1)%I,(ITPFGS(IGAS)-1)*NSMST+1,NSMST,NNSTSGP(1,IGAS))
  call ICOPVE2(ISTSGP(1)%I,(ITPFGS(IGAS)-1)*NSMST+1,NSMST,IISTSGP(1,IGAS))
end do

do IGAS=1,NGAS
  do ISMST=1,NSMST
    if (NNSTSGP(ISMST,IGAS) > 0) maxval(IGAS) = ISMST
  end do
  do ISMST=NSMST,1,-1
    if (NNSTSGP(ISMST,IGAS) > 0) minval(IGAS) = ISMST
  end do
end do
! Largest and lowest active symmetries for each GAS space
if (NTEST >= 200) then
  write(u6,*) ' Type of each GAS space'
  call IWRTMA(ITPFGS,1,NGAS,1,NGAS)
  write(u6,*) ' Number of elecs per GAS space'
  call IWRTMA(NELFGS,1,NGAS,1,NGAS)
end if

! Loop over symmetries of each GAS

MAXLEX = 0
IFIRST = 1
ISTRBS = 1
1000 continue
if (IFIRST == 1) then
  do IGAS=1,NGASL-1
    ISMFGS(IGAS) = minval(IGAS)
  end do
else
  ! Next distribution of symmetries in NGAS -1
  !     NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
  !call NXTNUM2(ISMFGS,NGASL-1,1,MAXVAL,NONEW)
  !call NXTNUM3(IOCA,NGAS,IGSMIN,IGSMAX,NONEW)
  call NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
  if (NONEW /= 0) goto 1001
end if
IFIRST = 0
if (NTEST >= 200) then
  write(u6,*) ' next symmetry of NGASL-1 spaces'
  call IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
end if
! Symmetry of NGASL -1 spaces given, symmetry of total space
ISTSMM1 = 1
do IGAS=1,NGASL-1
  !    SYMCOM(ITASK,IOBJ,I1,I2,I12)
  call SYMCOM(3,1,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
  ISTSMM1 = JSTSMM1
  !write(u6,*) ' ISTSMM1 : ',ISTSMM1
end do
! required sym of SPACE NGASL
call SYMCOM(2,1,ISTSMM1,ISMGSN,ISPGRPSM)
ISMFGS(NGASL) = ISMGSN

do IGAS=NGASL+1,NGAS
  ISMFGS(IGAS) = 1
end do
if (NTEST >= 200) then
  write(u6,*) ' Next symmetry distribution'
  call IWRTMA(ISMFGS,1,NGAS,1,NGAS)
end if
! Obtain all strings of this symmetry
call GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTR(1+NEL*(ISTRBS-1)),NSTR,NEL,NNSTSGP,IISTSGP)
! Reorder Info : Lexical => actual number
if (IDOREO /= 0) then
  ! Lexical number of NEL electrons
  ! Can be made smart by using common factor for first NGAS-1 spaces
  do JSTR=ISTRBS,ISTRBS+NSTR-1
    LEX = 1
    do IEL=1,NEL
      LEX = LEX+IZ(ISTR(IEL+NEL*(JSTR-1)),IEL)
    end do
    !write(u6,*) ' string'
    !call IWRTMA(ISTR(1,JSTR),1,NEL,1,NEL)
    !write(u6,*) ' JSTR and LEX ',JSTR,LEX

    MAXLEX = max(MAXLEX,LEX)
    IREO(LEX) = JSTR
  end do
end if

ISTRBS = ISTRBS+NSTR
! ready for next symmetry distribution
if (NGAS-1 /= 0) goto 1000
1001 continue
! End of loop over symmetry distributions
NSTR = ISTRBS-1

if (NTEST >= 100) then
  write(u6,*) ' Number of strings generated ',NSTR
  write(u6,*)
  write(u6,*) ' Strings :'
  write(u6,*)
  call PRTSTR(ISTR,NEL,NSTR)

  if (IDOREO /= 0) then
    write(u6,*) 'Largest Lexical number obtained ',MAXLEX
    write(u6,*) ' Reorder array'
    call IWRTMA(IREO,1,NSTR,1,NSTR)
  end if
end if

end subroutine GETSTR_TOTSM_SPGP
