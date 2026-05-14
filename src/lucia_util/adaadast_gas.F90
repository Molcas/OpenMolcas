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
! I1(KSTR) = ISTR if a+/a IORB a+/a JORB !KSTR> = +/-!ISTR>, ISTR is in
! ISPGP,ISM,IGRP.
! (numbering relative to TS start)
! Only excitations IOB >= JOB are included
! The orbitals are in GROUP-SYM IOBTP,IOBSM, JOBTP,JOBSM respectively,
! and IOB (JOB) is the first orbital to be used, and the number of orbitals
! to be checked is NIOB (NJOB).
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
! Jeppe Olsen, August of 95   (adadst)
!              November 1997 : annihilation added

use Symmetry_Info, only: Mul
use lucia_data, only: IBGPSTR, IBSPGPFTP, IOBPTS, ISPGPFTP, MXPNGAS, NELFGP, NELFSPGP, NELFTP, NELIS, NGAS, NGPSTR, NOBPT, NOCOB, &
                      NSTRKS, OCSTR, REO, Z, ZSCR
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IOB, IOBSM, IOBTP, NIOB, IAC, JOB, JOBSM, JOBTP, NJOB, JAC, ISPGP, ISM, ITP, KMIN, KMAX, LI1, &
                                 IFRST, KFRST, I12, K12
integer(kind=iwp), intent(out) :: I1(LI1*NIOB*NJOB), NK, IEND
real(kind=wp), intent(out) :: XI1S(LI1*NIOB*NJOB)
real(kind=wp), intent(in) :: SCLFAC
integer(kind=iwp) :: IDELTA, IDUM(1), IEL, IGRP, IIGRP, IIOB, ISPGPABS, ITRIVIAL, JDELTA, JEL, JGRP, JJGRP, JJOB, K1SM, &
                     KGRP(MXPNGAS), KSM, NELI, NELK, NSTRI, NSTRK
integer(kind=iwp), save :: NSTRI_

! Some dummy initializations
IIGRP = 0 ! jwk-cleanup
JJGRP = 0 ! jwk-cleanup

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ======================'
write(u6,*) ' ADAADST_GAS in service'
write(u6,*) ' ======================'
write(u6,*)
write(u6,*) ' IOB,IOBSM,IOBTP,IAC ',IOB,IOBSM,IOBTP,IAC
write(u6,*) ' JOB,JOBSM,JOBTP,JAC ',JOB,JOBSM,JOBTP,JAC
write(u6,*) ' I12, K12 ',I12,K12
write(u6,*) ' IFRST,KFRST',IFRST,KFRST
#endif

! Internal affairs

if ((I12 > size(Z,2)) .or. (K12 > size(OCSTR,2))) then
  write(u6,*) ' ADST_GAS : Illegal value of I12 or K12 ',I12,K12
  !stop ' ADST_GAS : Illegal value of I12 or K12'
  call SYSABENDMSG('lucia_util/adst_gas','Internal error','')
  return
end if

! Supergroup and symmetry of K strings

K1SM = Mul(IOBSM,ISM)
KSM = Mul(JOBSM,K1SM)
#ifdef _DEBUGPRINT_
write(u6,*) ' K1SM,KSM : ',K1SM,KSM
#endif
ISPGPABS = IBSPGPFTP(ITP)-1+ISPGP
IDELTA = -1
if (IAC == 2) IDELTA = 1
JDELTA = -1
if (JAC == 2) JDELTA = 1
#ifdef _DEBUGPRINT_
write(u6,*) ' IDELTA, JDELTA',IDELTA,JDELTA
#endif
! Occupation of K-strings
if (IOBTP == JOBTP) then
  IEL = NELFSPGP(IOBTP,ISPGPABS)-IDELTA-JDELTA
  JEL = IEL
else
  IEL = NELFSPGP(IOBTP,ISPGPABS)-IDELTA
  JEL = NELFSPGP(JOBTP,ISPGPABS)-JDELTA
end if
#ifdef _DEBUGPRINT_
write(u6,*) ' IEL, JEL',IEL,JEL
#endif
! Trivial zero? (Nice, then mission is complete)
ITRIVIAL = 0
if ((IEL < 0) .or. (JEL < 0) .or. (IEL > NOBPT(IOBTP)) .or. (JEL > NOBPT(JOBTP))) then
  ! No strings with this number of elecs - be happy : No work
  NK = 0
# ifdef _DEBUGPRINT_
  write(u6,*) ' Trivial zero excitations'
# endif
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
  KGRP(1:NGAS) = ISPGPFTP(1:NGAS,ISPGPABS)
  KGRP(IOBTP) = IIGRP
  KGRP(JOBTP) = JJGRP
# ifdef _DEBUGPRINT_
  write(u6,*) ' Groups in KGRP'
  call IWRTMA(KGRP,1,NGAS,1,NGAS)
# endif
end if

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
# ifdef _DEBUGPRINT_
  write(u6,*) ' Info on I strings generated'
  write(u6,*) ' NSTRI = ',NSTRI
  write(u6,*) ' REORDER array'
  call IWRTMA(REO(:,I12),1,NSTRI,1,NSTRI)
# endif
  NSTRI_ = NSTRI

end if
#ifdef _DEBUGPRINT_
write(u6,*) ' REORDER array for I STRINGS'
call IWRTMA(REO(:,I12),1,NSTRI,1,NSTRI)
#endif

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
#ifdef _DEBUGPRINT_
write(u6,*) ' NELK = ',NELK
#endif
if (KFRST /= 0) then
  ! Generate occupation of K STRINGS
  IDUM(1) = 0
  call GETSTR2_TOTSM_SPGP(KGRP,NGAS,KSM,NELK,NSTRK,OCSTR(:,K12),NOCOB,0,IDUM,IDUM)
  !    GETSTR2_TOTSM_SPGP(IGRP,NIGRP,ISPGRPSM,NEL,NSTR,ISTR,NORBT,IDOREO,IZ,IREO)
  NSTRKS(K12) = NSTRK
# ifdef _DEBUGPRINT_
  write(u6,*) ' K strings generated'
  write(u6,*) ' Reorder array after generation of K strings'
  call IWRTMA(REO(:,I12),1,NSTRI,1,NSTRI)
# endif
end if

NSTRK = NSTRKS(K12)

IIOB = IOBPTS(IOBTP,IOBSM)+IOB-1
JJOB = IOBPTS(JOBTP,JOBSM)+JOB-1

I1(:) = 0
!OLD XI1S(:) = Zero

call ADAADAS1_GAS(NK,I1,XI1S,LI1,IIOB,NIOB,IAC,JJOB,NJOB,JAC,OCSTR(:,K12),NELK,NSTRK,REO(:,I12),Z(:,I12),NOCOB,KMAX,KMIN,IEND, &
                  SCLFAC,NSTRI_)

#ifdef _DEBUGPRINT_
write(u6,*) ' Reorder array after ADAADAS1'
call IWRTMA(REO(:,I12),1,NSTRI,1,NSTRI)
#endif

end subroutine ADAADAST_GAS
