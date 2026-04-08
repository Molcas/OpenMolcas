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
! Copyright (C) 1995,1997, Jeppe Olsen                                 *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,NORBT,IDOREO,IZ,IREO)
! Obtain all super-strings of given total symmetry and given
! occupation in each GAS space
!
! If IDOREO /= 0 THEN reordering array : lexical => actual order is obtained
!
! Nomenclature of the day : superstring : string in complete
!                           orbital space, product of strings in
!                           each GAS space
!
! Compared to GETSTR2_TOTSM_SPGP : Based upon IGRP(NIGRP)
!                                  (Just a few changes in the beginning)
!
! =====
! Input
! =====
!
! IGRP :  supergroup, here as an array of GAS space
! NIGRP : Number of active groups
! ISPGRPSM : Total symmetry of superstrings
! NEL : Number of electrons
! IZ  : Reverse lexical ordering array for this supergroup (IF IDOREO /= 0)
!
! ======
! Output
! ======
!
! NSTR : Number of superstrings generated
! ISTR : Occupation of superstring
! IREO : Reorder array (if IDOREO /= 0)
!
! Jeppe Olsen, Written  July 1995
!              Version of Dec 1997

use Symmetry_Info, only: Mul
use lucia_data, only: ISTSGP, MXPNGAS, MXPNSMST, NELFGP, NGAS, NIRREP, NSTSGP
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NIGRP, IGRP(NIGRP), ISPGRPSM, NEL, NORBT, IDOREO, IZ(NORBT,NEL)
integer(kind=iwp), intent(out) :: NSTR
integer(kind=iwp), intent(inout) :: ISTR(*)
integer(kind=iwp), intent(_OUT_) :: IREO(*)
integer(kind=iwp) :: IEL, IFIRST, IGAS, IISTSGP(MXPNSMST,MXPNGAS), ISMFGS(MXPNGAS), ISMST, ISTRBS, ISTSMM1, ITPFGS(MXPNGAS), JSTR, &
                     LEX, MAXLEX, MNVAL(MXPNGAS), MXVAL(MXPNGAS), NELFGS(MXPNGAS), NGASL, NNSTSGP(MXPNSMST,MXPNGAS), NONEW
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I

write(u6,*)
write(u6,*) ' ============================'
write(u6,*) ' Welcome to GETSTR_TOTSM_SPGP'
write(u6,*) ' ============================'
write(u6,*)
write(u6,'(A)') ' Strings to be obtained :'
write(u6,'(A)') ' ************************'
write(u6,'(A)')
write(u6,'(A,I2)') '   Symmetry : ',ISPGRPSM
write(u6,'(A,16I3)') ' Groups : ',(IGRP(I),I=1,NIGRP)
write(u6,*) ' NEL = ',NEL
if (IDOREO /= 0) then
  write(u6,*)
  write(u6,*) ' ============='
  write(u6,*) ' The Z array :'
  write(u6,*) ' ============='
  write(u6,*)
  write(u6,*) ' NORBT,NEL = ',NORBT,NEL
  call IWRTMA(IZ,NORBT,NEL,NORBT,NEL)
end if
#endif
! Absolute number of this supergroup
! Occupation per gasspace
! Largest occupied space
NGASL = 0
! Largest and lowest symmetries active in each GAS space
do IGAS=1,NGAS
  ITPFGS(IGAS) = IGRP(IGAS)
  NELFGS(IGAS) = NELFGP(IGRP(IGAS))
  if (NELFGS(IGAS) > 0) NGASL = IGAS
end do
if (NGASL == 0) NGASL = 1
! Number of strings per GAS space and offsets for strings of given sym
do IGAS=1,NGAS
  IISTSGP(1:NIRREP,IGAS) = ISTSGP((ITPFGS(IGAS)-1)*NIRREP+1:ITPFGS(IGAS)*NIRREP)
  NNSTSGP(1:NIRREP,IGAS) = NSTSGP((ITPFGS(IGAS)-1)*NIRREP+1:ITPFGS(IGAS)*NIRREP)
end do

do IGAS=1,NGAS
  do ISMST=1,NIRREP
    if (NNSTSGP(ISMST,IGAS) > 0) MXVAL(IGAS) = ISMST
  end do
  do ISMST=NIRREP,1,-1
    if (NNSTSGP(ISMST,IGAS) > 0) MNVAL(IGAS) = ISMST
  end do
end do
! Largest and lowest active symmetries for each GAS space
#ifdef _DEBUGPRINT_
write(u6,*) ' Type of each GAS space'
call IWRTMA(ITPFGS,1,NGAS,1,NGAS)
write(u6,*) ' Number of elecs per GAS space'
call IWRTMA(NELFGS,1,NGAS,1,NGAS)
#endif

! Loop over symmetries of each GAS

MAXLEX = 0
IFIRST = 1
ISTRBS = 1
do
  if (IFIRST == 1) then
    ISMFGS(1:NGASL-1) = MNVAL(1:NGASL-1)
  else
    ! Next distribution of symmetries in NGAS -1
    call NXTNUM3(ISMFGS,NGASL-1,MNVAL,MXVAL,NONEW)
    if (NONEW /= 0) exit
  end if
  IFIRST = 0
# ifdef _DEBUGPRINT_
  write(u6,*) ' next symmetry of NGASL-1 spaces'
  call IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
# endif
  ! Symmetry of NGASL -1 spaces given, symmetry of total space
  ISTSMM1 = 1
  do IGAS=1,NGASL-1
    ISTSMM1 = Mul(ISTSMM1,ISMFGS(IGAS))
  end do
  ! required sym of SPACE NGASL
  ISMFGS(NGASL) = Mul(ISTSMM1,ISPGRPSM)

  ISMFGS(NGASL+1:NGAS) = 1
# ifdef _DEBUGPRINT_
  write(u6,*) ' Next symmetry distribution'
  call IWRTMA(ISMFGS,1,NGAS,1,NGAS)
# endif
  ! Obtain all strings of this symmetry
  !call GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTR(1,ISTRBS),NSTR,NEL,NNSTSGP,IISTSGP)
  call GETSTRN_GASSM_SPGP(ISMFGS,ITPFGS,ISTR(1+NEL*(ISTRBS-1)),NSTR,NEL,NNSTSGP,IISTSGP)
  ! Reorder Info : Lexical => actual number
  if (IDOREO /= 0) then
    ! Lexical number of NEL electrons
    ! Can be made smart by using common factor for first NGAS-1 spaces
    do JSTR=ISTRBS,ISTRBS+NSTR-1
      LEX = 1
      do IEL=1,NEL
        !LEX = LEX+IZ(ISTR(IEL,JSTR)),IEL)
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
  if (NGAS == 1) exit
end do
! End of loop over symmetry distributions
NSTR = ISTRBS-1

#ifdef _DEBUGPRINT_
write(u6,*) ' NEL(b) = ',NEL
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
#endif

end subroutine GETSTR2_TOTSM_SPGP
