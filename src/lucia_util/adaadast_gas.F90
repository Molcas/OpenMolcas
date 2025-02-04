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

subroutine ADAADAST_GAS(IOB,IOBSM,IOBTP,NIOB,IAC,JOB,JOBSM,JOBTP,NJOB,JAC,ISPGP,ISM,ITP,KMIN,KMAX,I1,XI1S,LI1,NK,IEND,IFRST,KFRST, &
                        I12,K12,SCLFAC)
! Obtain two-operator mappings
! a+/a IORB a+/a JORB !KSTR> = +/-!ISTR>
!
! Whether creation- or annihilation operators are in use depends
! upon IAC, JAC : 1=> Annihilation,
!                 2=> Creation
!
! In the form
! I1(KSTR) =  ISTR if a+/a IORB a+/a JORB !KSTR> = +/-!ISTR>, ISTR is in
! ISPGP,ISM,IGRP.
! (numbering relative to TS start)
! Only excitations IOB >= JOB are included
! The orbitals are in GROUP-SYM IOBTP,IOBSM, JOBTP,JOBSM respectively,
! and IOB (JOB) is the first orbital to be used, and the number of orbitals
! to be checked is NIOB ( NJOB).
!
! Only orbital pairs IOB > JOB are included (if the types are identical)
!
! The output is given in I1(KSTR,I,J) = I1 ((KSTR,(J-1)*NIOB + I)
!
! Above +/- is stored in XI1S
! Number of K strings checked is returned in NK
! Only Kstrings with relative numbers from KMIN to KMAX are included
!
! If IEND /= 0 last string has been checked
!
! Jeppe Olsen, August of 95   ( adadst)
!              November 1997 : annihilation added

use HIDSCR, only: ZSCR, ZOCSTR => OCSTR, REO, Z
use lucia_data, only: NGAS
use lucia_data, only: IBGPSTR, IBSPGPFTP, ISPGPFTP, NELFGP, NELFSPGP, NELFTP, NGPSTR
use lucia_data, only: NELIS, NSTRKS
use lucia_data, only: IOBPTS, NOBPT, NOCOB
use lucia_data, only: MXPNGAS
use Definitions, only: u6

implicit none
integer IOB, IOBSM, IOBTP, NIOB, IAC, JOB, JOBSM, JOBTP, NJOB, JAC, ISPGP, ISM, ITP, KMIN, KMAX, LI1, NK, IEND, IFRST, KFRST, I12, &
        K12
real*8 SCLFAC
integer I1(*)
real*8 XI1S(*)
integer KGRP(MXPNGAS)
integer IDUM(1)
integer, save :: NSTRI_
integer IIGRP, JJGRP, NTEST, K1SM, KSM, ISPGPABS, IACADJ, IDELTA, JACADJ, JDELTA, IEL, JEL, ITRIVIAL, IGRP, JGRP, NTEST2, NELI, &
        NSTRI, NELK, NSTRK, IIOB, JJOB, IZERO

! Some dummy initializations
IIGRP = 0 ! jwk-cleanup
JJGRP = 0 ! jwk-cleanup

NTEST = 0
if (NTEST >= 100) then
  write(u6,*)
  write(u6,*) ' ======================'
  write(u6,*) ' ADAADST_GAS in service'
  write(u6,*) ' ======================'
  write(u6,*)
  write(u6,*) ' IOB,IOBSM,IOBTP,IAC ',IOB,IOBSM,IOBTP,IAC
  write(u6,*) ' JOB,JOBSM,JOBTP,JAC ',JOB,JOBSM,JOBTP,JAC
  write(u6,*) ' I12, K12 ',I12,K12
  write(u6,*) ' IFRST,KFRST',IFRST,KFRST
end if

! Internal affairs

if ((I12 > size(Z,2)) .or. (K12 > size(ZOCSTR,2))) then
  write(u6,*) ' ADST_GAS : Illegal value of I12 or K12 ',I12,K12
  !stop ' ADST_GAS : Illegal value of I12 or K12'
  call SYSABENDMSG('lucia_util/adst_gas','Internal error','')
  return
end if

! Supergroup and symmetry of K strings

call SYMCOM(2,0,IOBSM,K1SM,ISM)
call SYMCOM(2,0,JOBSM,KSM,K1SM)
if (NTEST >= 100) write(u6,*) ' K1SM,KSM : ',K1SM,KSM
ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
IACADJ = 2
IDELTA = -1
if (IAC == 2) then
  IACADJ = 1
  IDELTA = 1
end if
JACADJ = 2
JDELTA = -1
if (JAC == 2) then
  JACADJ = 1
  JDELTA = 1
end if
if (NTEST >= 100) then
  write(u6,*) ' IACADJ, JACADJ',IACADJ,JACADJ
  write(u6,*) ' IDELTA, JDELTA',IDELTA,JDELTA
end if
! Occupation of K-strings
if (IOBTP == JOBTP) then
  IEL = NELFSPGP(IOBTP,ISPGPABS)-IDELTA-JDELTA
  JEL = IEL
else
  IEL = NELFSPGP(IOBTP,ISPGPABS)-IDELTA
  JEL = NELFSPGP(JOBTP,ISPGPABS)-JDELTA
end if
if (NTEST >= 100) write(u6,*) ' IEL, JEL',IEL,JEL
! Trivial zero ? (Nice, then mission is complete )
ITRIVIAL = 0
if ((IEL < 0) .or. (JEL < 0) .or. (IEL > NOBPT(IOBTP)) .or. (JEL > NOBPT(JOBTP))) then
  ! No strings with this number of elecs - be happy : No work
  NK = 0
  if (NTEST >= 100) write(u6,*) ' Trivial zero excitations'
  ITRIVIAL = 1
  !return
else
  ! Find group with IEL electrons in IOBTP, JEL in JOBTP
  IIGRP = 0
  do IGRP=IBGPSTR(IOBTP),IBGPSTR(IOBTP)+NGPSTR(IOBTP)-1
    if (NELFGP(IGRP) == IEL) IIGRP = IGRP
  end do
  JJGRP = 0
  do JGRP=IBGPSTR(JOBTP),IBGPSTR(JOBTP)+NGPSTR(JOBTP)-1
    if (NELFGP(JGRP) == JEL) JJGRP = JGRP
  end do
  !write(u6,*) ' ADAADA : IIGRP, JJGRP',IIGRP,JJGRP

  if ((IIGRP == 0) .or. (JJGRP == 0)) then
    write(u6,*) ' ADAADAST : cul de sac, active K groups not found'
    write(u6,*) ' Active GAS spaces  ',IOBTP,JOBTP
    write(u6,*) ' Number of electrons',IEL,JEL
    !stop ' ADAADAST : cul de sac, active K groups not found'
    call SYSABENDMSG('lucia_util/adaadast_gas','Internal error','')
  end if

end if
! Groups defining Kstrings
if (ITRIVIAL /= 1) then
  call ICOPVE(ISPGPFTP(1,ISPGPABS),KGRP,NGAS)
  KGRP(IOBTP) = IIGRP
  KGRP(JOBTP) = JJGRP
  if (NTEST >= 100) then
    write(u6,*) ' Groups in KGRP'
    call IWRTMA(KGRP,1,NGAS,1,NGAS)
  end if
end if

! In ADADS1_GAS we need : Occupation of KSTRINGS
!                         lexical => Actual order for I strings
! Generate if required

if (IFRST /= 0) then
  ! Generate information about I strings
  ! Arc weights for ISPGP
  NTEST2 = NTEST
  call WEIGHT_SPGP(Z(:,I12),NGAS,NELFSPGP(1,ISPGPABS),NOBPT,ZSCR,NTEST2)
  NELI = NELFTP(ITP)
  NELIS(I12) = NELI
  ! Reorder array for I strings
  call GETSTR_TOTSM_SPGP(ITP,ISPGP,ISM,NELI,NSTRI,ZOCSTR(:,K12),NOCOB,1,Z(:,I12),REO(:,I12))
  if (NTEST >= 1000) then
    write(u6,*) ' Info on I strings generated'
    write(u6,*) ' NSTRI = ',NSTRI
    write(u6,*) ' REORDER array'
    call IWRTMA(REO(:,I12),1,NSTRI,1,NSTRI)
  end if
  NSTRI_ = NSTRI

end if
if (NTEST >= 1000) then
  write(u6,*) ' REORDER array for I STRINGS'
  call IWRTMA(REO(:,I12),1,NSTRI,1,NSTRI)
end if

if (ITRIVIAL == 1) return
NELK = NELIS(I12)
if (IAC == 1) then
  NELK = NELK+1
else
  NELK = NELK-1
end if
if (JAC == 1) then
  NELK = NELK+1
else
  NELK = NELK-1
end if
if (NTEST >= 100) write(u6,*) ' NELK = ',NELK
if (KFRST /= 0) then
  ! Generate occupation of K STRINGS
  IDUM(1) = 0
  call GETSTR2_TOTSM_SPGP(KGRP,NGAS,KSM,NELK,NSTRK,ZOCSTR(:,K12),NOCOB,0,IDUM,IDUM)
  !    GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,NORBT,IDOREO,IZ,IREO)
  NSTRKS(K12) = NSTRK
  if (NTEST >= 1000) then
    write(u6,*) ' K strings generated'
    write(u6,*) ' Reorder array after generation of K strings'
    call IWRTMA(REO(:,I12),1,NSTRI,1,NSTRI)
  end if
end if

NSTRK = NSTRKS(K12)

IIOB = IOBPTS(IOBTP,IOBSM)+IOB-1
JJOB = IOBPTS(JOBTP,JOBSM)+JOB-1

IZERO = 0
call ISETVC(I1,IZERO,LI1*NIOB*NJOB)
!OLD call SETVEC(XI1S,Zero,LI1*NIOB*NJOB)

call ADAADAS1_GAS(NK,I1,XI1S,LI1,IIOB,NIOB,IAC,JJOB,NJOB,JAC,ZOCSTR(:,K12),NELK,NSTRK,REO(:,I12),Z(:,I12),NOCOB,KMAX,KMIN,IEND, &
                  SCLFAC,NSTRI_)

if (NTEST >= 1000) then
  write(u6,*) ' Reorder array after ADAADAS1'
  call IWRTMA(REO(:,I12),1,NSTRI,1,NSTRI)
end if

end subroutine ADAADAST_GAS
