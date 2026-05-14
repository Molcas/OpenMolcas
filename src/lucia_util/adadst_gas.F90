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

!#define _DEBUGPRINT_
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
! to be checked is NIOB (NJOB).
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

use Symmetry_Info, only: Mul
use lucia_data, only: IBSPGPFTP, IOBPTS, NELFSPGP, NELFTP, NELIS, NGAS, NOBPT, NOCOB, NSTRKS, OCSTR, REO, Z, ZSCR
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IOB, IOBSM, IOBTP, NIOB, JOB, JOBSM, JOBTP, NJOB, ISPGP, ISM, ITP, KMIN, KMAX, LI1, IFRST, KFRST, &
                                 I12, K12
integer(kind=iwp), intent(out) :: I1(LI1,NIOB*NJOB), NK, IEND
real(kind=wp), intent(out) :: XI1S(LI1,NIOB*NJOB)
real(kind=wp), intent(in) :: SCLFAC
integer(kind=iwp) :: IDUM_ARR(1), IIOB, ISPGPABS, JJOB, K1SM, K1SPGPABS, KSM, KSPGPABS, NELI, NELK, NSTRI, NSTRK

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ====================='
write(u6,*) ' ADADST_GAS in service'
write(u6,*) ' ====================='
write(u6,*)
write(u6,*) ' IOB,IOBSM,IOBTP ',IOB,IOBSM,IOBTP
write(u6,*) ' JOB,JOBSM,JOBTP ',JOB,JOBSM,JOBTP
#endif

!if (SCLFAC /= One) then
!  write(u6,*) 'Problemo, ADADST'
!  write(u6,*) ' SCLFAC = ',SCLFAC
!end if

! Internal affairs

if ((I12 > size(Z,2)) .or. (K12 > size(OCSTR,2))) then
  write(u6,*) ' ADST_GAS : Illegal value of K12 = ',K12
  write(u6,*) ' ADST_GAS : Illegal value of I12 = ',I12
  !stop ' ADST_GAS : Illegal value of I12'
  call SYSABENDMSG('lucia_util/adst_gas','Internal error','')
  return
end if

! Supergroup and symmetry of K strings

ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
call NEWTYP(ISPGPABS,1,IOBTP,K1SPGPABS)
call NEWTYP(K1SPGPABS,1,JOBTP,KSPGPABS)
K1SM = Mul(IOBSM,ISM)
KSM = Mul(JOBSM,K1SM)
#ifdef _DEBUGPRINT_
write(u6,*) ' K1SM,K1SPGPABS,KSM,KSPGPABS : ',K1SM,K1SPGPABS,KSM,KSPGPABS
#endif
! In ADADS1_GAS we need : Occupation of KSTRINGS
!                         lexical => Actual order for I strings
! Generate if required

if (IFRST /= 0) then
  ! Generate information about I strings
  ! Arc weights for ISPGP
  call WEIGHT_SPGP(Z(:,I12),NGAS,NELFSPGP(1,ISPGPABS),NOBPT,ZSCR)
  NELI = NELFTP(ITP)
  NELIS(I12) = NELI
  ! Reorder array for I strings
  call GETSTR_TOTSM_SPGP(ITP,ISPGP,ISM,NELI,NSTRI,OCSTR(:,K12),NOCOB,1,Z(:,I12),REO(:,I12))
end if
NELK = NELIS(I12)-2
if (KFRST /= 0) then
  ! Generate occupation of K STRINGS
  call GETSTR_TOTSM_SPGP(1,KSPGPABS,KSM,NELK,NSTRK,OCSTR(:,K12),NOCOB,0,IDUM_ARR,IDUM_ARR)
  NSTRKS(K12) = NSTRK
end if

NSTRK = NSTRKS(K12)

IIOB = IOBPTS(IOBTP,IOBSM)+IOB-1
JJOB = IOBPTS(JOBTP,JOBSM)+JOB-1
call ADADS1_GAS(NK,I1,XI1S,LI1,IIOB,NIOB,JJOB,NJOB,OCSTR(:,K12),NELK,NSTRK,REO(:,I12),Z(:,I12),NOCOB,KMAX,KMIN,IEND,SCLFAC)

end subroutine ADADST_GAS
