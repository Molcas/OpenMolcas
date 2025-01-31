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

subroutine ADADST_GAS(IOB,IOBSM,IOBTP,NIOB,JOB,JOBSM,JOBTP,NJOB,ISPGP,ISM,ITP,KMIN,KMAX,I1,XI1S,LI1,NK,IEND,IFRST,KFRST,I12,K12, &
                      SCLFAC)
! Obtain mappings
! a+IORB a+ JORB !KSTR> = +/-!ISTR>
! In the form
! I1(KSTR) =  ISTR if a+IORB a+ JORB !KSTR> = +/-!ISTR>, ISTR is in
! ISPGP,ISM,IGRP.
! (numbering relative to TS start)
! Only excitations IOB >= JOB are included
! The orbitals are in GROUP-SYM IOBTP,IOBSM, JOBTP,JOBSM respectively,
! and IOB (JOB) is the first orbital to be used, and the number of orbitals
! to be checked is NIOB ( NJOB).
!
! Only orbital pairs IOB > JOB are included
!
! The output is given in I1(KSTR,I,J) = I1 ((KSTR,(J-1)*NIOB + I)
!
! Above +/- is stored in XI1S
! Number of K strings checked is returned in NK
! Only Kstrings with relative numbers from KMIN to KMAX are included
!
! If IEND /= 0 last string has been checked
!
! Jeppe Olsen, August of 95

use HIDSCR, only: ZSCR, ZOCSTR => OCSTR, REO, Z
use lucia_data, only: NGAS
use lucia_data, only: IBSPGPFTP, NELFSPGP, NELFTP
use lucia_data, only: NELIS, NSTRKS
use lucia_data, only: IOBPTS, NOBPT, NOCOB

implicit none
integer IOB, IOBSM, IOBTP, NIOB, JOB, JOBSM, JOBTP, NJOB, ISPGP, ISM, ITP, KMIN, KMAX, LI1, NK, IEND, IFRST, KFRST, I12, K12
real*8 SCLFAC
integer I1(*)
real*8 XI1S(*)
integer IDUM_ARR(1)
integer NTEST, ISPGPABS, K1SM, K1SPGPABS, KSM, KSPGPABS, NTEST2, NELI, NELK, NSTRK, IIOB, JJOB, NSTRI

NTEST = 0
if (NTEST >= 100) then
  write(6,*)
  write(6,*) ' ====================='
  write(6,*) ' ADADST_GAS in service'
  write(6,*) ' ====================='
  write(6,*)
  write(6,*) ' IOB,IOBSM,IOBTP ',IOB,IOBSM,IOBTP
  write(6,*) ' JOB,JOBSM,JOBTP ',JOB,JOBSM,JOBTP
end if

!if (SCLFAC /= 1.0D0) then
!  write(6,*) 'Problemo, ADADST'
!  write(6,*) ' SCLFAC = ',SCLFAC
!end if

! Internal affairs

if ((I12 > size(Z,2)) .or. (K12 > size(ZOCSTR,2))) then
  write(6,*) ' ADST_GAS : Illegal value of K12 = ',K12
  write(6,*) ' ADST_GAS : Illegal value of I12 = ',I12
  !stop ' ADST_GAS : Illegal value of I12'
  call SYSABENDMSG('lucia_util/adst_gas','Internal error','')
  return
end if

! Supergroup and symmetry of K strings

ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
call NEWTYP(ISPGPABS,1,IOBTP,K1SPGPABS)
call NEWTYP(K1SPGPABS,1,JOBTP,KSPGPABS)
call SYMCOM(2,0,IOBSM,K1SM,ISM)
call SYMCOM(2,0,JOBSM,KSM,K1SM)
if (NTEST >= 100) write(6,*) ' K1SM,K1SPGPABS,KSM,KSPGPABS : ',K1SM,K1SPGPABS,KSM,KSPGPABS
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
end if
NELK = NELIS(I12)-2
if (KFRST /= 0) then
  ! Generate occupation of K STRINGS
  call GETSTR_TOTSM_SPGP(1,KSPGPABS,KSM,NELK,NSTRK,ZOCSTR(:,K12),NOCOB,0,IDUM_ARR,IDUM_ARR)
  NSTRKS(K12) = NSTRK
end if

NSTRK = NSTRKS(K12)

IIOB = IOBPTS(IOBTP,IOBSM)+IOB-1
JJOB = IOBPTS(JOBTP,JOBSM)+JOB-1
call ADADS1_GAS(NK,I1,XI1S,LI1,IIOB,NIOB,JJOB,NJOB,ZOCSTR(:,K12),NELK,NSTRK,REO(:,I12),Z(:,I12),NOCOB,KMAX,KMIN,IEND,SCLFAC)

end subroutine ADADST_GAS
