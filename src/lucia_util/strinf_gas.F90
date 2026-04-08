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
! Copyright (C) 1994, Jeppe Olsen                                      *
!               2024, Giovanni Li Manni                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine STRINF_GAS()
! Obtain string information for GAS expansion
!
! =====
! Input
! =====
!
! /LUCINP/,/ORBINP/,/CSM/, /CGAS/, /GASSTR/
!
! =====
! Output
! =====
!
! /STRINP/,/STINF/,STRBAS and string information in STIN
!
! modification Jeppe + Giovanni + Dongxia.
! G. Li Manni, June 2024: Scale-up capability for single SD ROHF type calculations

use lucia_data, only: IBSPGPFTP, IGSFGP, IGSOCC, INGRP_VAL, IOCLS, IPHGAS, ISMDFGP, ISMSCR, ISPGPFTP, ISTAC, ISTSGP, ISTSO, &
                      MAX_STR_OC_BLK, MAX_STR_SPGP, MINMAX_SM_GP, MS2, MXNSTR, MXPNGAS, MXPNSMST, NACOB, NACTEL, NACTSYM, NELFGP, &
                      NELFSPGP, NELFTP, NGAS, NGRP, NHLFSPGP, NIRREP, NMXOCCLS, NOBPT, NORB1, NORB2, NORB3, NSPGPFTP, NSTFGP, &
                      NSTFSMGP, NSTFSMSPGP, NSTSGP, NSTSO, NSTTYP, NTSPGP, OCCSTR, SPGPAN, SPGPCR, STREO, STSTM, ZMAT
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: i, IEC, IEL, IGAS, IGP, IGRP, IGRPABS, IGSOB, IIEL, IISPGP, IOCTYPX, ISM, ISPGP, ISTTYP, ISTTYPA, ISTTYPC, &
                     ITP, JGRP, JSTTYP, LAC, LROW, MAXSCR, MN, MNRS1X, MNRS3X, MX, MXNSTRFSG, MXRS1X, MXRS3X, NACOB_EFFECTIVE, &
                     NEL, NGSOBP, NHOLE, NSMCLS, NSMCLSE, NSMCLSE1, NSTINI, NSTR, ZERO_ARR(1)
integer(kind=iwp), allocatable :: FREEL(:)

! Some dummy initializations
LAC = 0 ! jwk-cleanup

! 2 : Number of classes per string type and mappings between
!     string types (/STINF/)

if ((NActEl /= MS2) .or. (NActEl /= NACOB)) call ZSTINF_GAS()

! 3 : Static memory for string information

call MEMSTR_GAS()

!  4 : Info about group of strings

! First free address

! Find maximum needed length of scratch.
!   MAXSCR = 2*NACOB+(IEL+1)(NACOB+1)
!   with IEL = MAX(NELFGP(IGRP=1,NGRP)

IEL = max(0,maxval(NELFGP(1:NGRP)))
NACOB_EFFECTIVE = NACOB
if (NACOB == 0) NACOB_EFFECTIVE = 1
MAXSCR = 2*NACOB_EFFECTIVE+(IEL+1)*(NACOB_EFFECTIVE+1)+NIRREP
call mma_allocate(FREEL,MAXSCR,Label='FREEL')
do IGRP=1,NGRP
  ! A gas group can be considered as a RAS group with 0 electrons in RAS1, RAS3!
  IGAS = IGSFGP(IGRP)
  if (IGAS == 1) then
    NORB1 = 0
  else
    NORB1 = sum(NOBPT(1:IGAS-1))
  end if
  NORB2 = NOBPT(IGAS)
  NORB3 = NACOB-NORB1-NORB2
  MNRS1X = 0
  MXRS1X = 0
  MNRS3X = 0
  MXRS3X = 0
  IEL = NELFGP(IGRP)
  IOCTYPX = 1
  ! Reverse lexical adresing schemes for each group of string
  call WEIGHT_LUCIA(Zmat(IGRP)%A,IEL,NORB1,NORB2,NORB3,MNRS1X,MXRS1X,MNRS3X,MXRS3X,FREEL)
  ! Number of strings per symmetry in a given group
  call NSTRSO_GAS(IEL,NORB1,NORB2,NORB3,MNRS1X,MXRS1X,MNRS3X,MXRS3X,FREEL,NSTSGP,ISTSGP,NIRREP,IGRP)
  ! Construct the strings ordered by symmetry
  call GENSTR_GAS(IEL,MNRS1X,MXRS1X,MNRS3X,MXRS3X,ISTSGP,IGRP,IOCTYPX,NIRREP,Zmat(IGRP)%A,FREEL,STREO(IGRP)%A,OCCSTR(IGRP)%A, &
                  FREEL(1+IOCTYPX*NIRREP))

  NSTFSMGP(1:NIRREP,IGRP) = NSTSGP((IGRP-1)*NIRREP+1:IGRP*NIRREP)
  !ISTFSMGP(1:NIRREP,IGRP) = ISTSGP((IGRP-1)*NIRREP+1:IGRP*NIRREP)
end do
call mma_deallocate(FREEL)

INGRP_VAL = NGRP
call mma_allocate(ISMDFGP,NIRREP*NGRP,Label='ISMDFGP')
call mma_allocate(NACTSYM,NGRP,Label='NACTSYM')
call mma_allocate(ISMSCR,NGRP,Label='ISMSCR')
call SMDFGP_GEN(NGRP,NIRREP,MXPNSMST,NSTFSMGP,NACTSYM,ISMDFGP)

#ifdef _DEBUGPRINT_
write(u6,*) 'NGRP',NGRP
write(u6,*) 'NIRREP*NGRP',NIRREP*NGRP
write(u6,*) ' Number of strings per group and symmetry'
call IWRTMA10(NSTSGP,NIRREP,NGRP,NIRREP,NGRP)
write(u6,*) ' Number of strings per group and symmetry(2)'
call IWRTMA10(NSTFSMGP,NIRREP,NGRP,MXPNSMST,NGRP)
#endif

! Min and max of sym for each group

do IGP=1,NGRP
  MX = 1
  do ISM=1,NIRREP
    if (NSTFSMGP(ISM,IGP) > 0) MX = ISM
  end do

  MN = NIRREP
  do ISM=NIRREP,1,-1
    if (NSTFSMGP(ISM,IGP) > 0) MN = ISM
  end do

  MINMAX_SM_GP(1,IGP) = MN
  MINMAX_SM_GP(2,IGP) = MX

end do
#ifdef _DEBUGPRINT_
write(u6,*) ' MINMAX array for sym of groups'
write(u6,*) ' =============================='
call IWRTMA(MINMAX_SM_GP,2,NGRP,2,NGRP)
#endif

! 4.5 : Creation/Annihilation mappings between different
!       types of strings

do IGRP=1,NGRP

  IGAS = IGSFGP(IGRP)
  NGSOBP = NOBPT(IGAS)
  ! First orbital in GAS spacce
  IGSOB = sum(NOBPT(1:IGAS-1))+1
  IEL = NELFGP(IGRP)
  NSTINI = NSTFGP(IGRP)

  ! Type of mapping : Only creation                  (LAC = 1)
  !                   Only annihilation              (LAC = 2)
  !                   Both annihilation and creation (LAC = 3)
  ! If only annihilation is present the string mapping arrays
  ! will only be over electronns
  if ((ISTAC(IGRP,1) /= 0) .and. (ISTAC(IGRP,2) /= 0)) then
    LAC = 3
    IEC = 1
    LROW = NGSOBP
  else if ((ISTAC(IGRP,1) /= 0) .and. (ISTAC(IGRP,2) == 0)) then
    LAC = 1
    IEC = 2
    LROW = IEL
  else if ((ISTAC(IGRP,1) == 0) .and. (ISTAC(IGRP,2) /= 0)) then
    LAC = 2
    IEC = 0
    LROW = NGSOBP
  else if ((ISTAC(IGRP,1) == 0) .and. (ISTAC(IGRP,2) == 0)) then
    LAC = 0
    IEC = 0
    LROW = 0
  end if
  ! Zero
  if (LAC /= 0) then
    STSTM(IGRP,1)%A(1:LROW*NSTINI) = 0
    STSTM(IGRP,2)%A(1:LROW*NSTINI) = 0
  end if

  if (ISTAC(IGRP,2) /= 0) then
    JGRP = ISTAC(IGRP,2)
    call CRESTR_GAS(OCCSTR(IGRP)%A,NSTFGP(IGRP),NSTFGP(JGRP),IEL,NGSOBP,IGSOB,Zmat(JGRP)%A,STREO(JGRP)%A,STSTM(IGRP,1)%A, &
                    STSTM(IGRP,2)%A,NACOB)

  end if
  if (ISTAC(IGRP,1) /= 0) then
    JGRP = ISTAC(IGRP,1)
    call ANNSTR_GAS(OCCSTR(IGRP)%A,NSTFGP(IGRP),NSTFGP(JGRP),IEL,NGSOBP,IGSOB,Zmat(JGRP)%A,STREO(JGRP)%A,STSTM(IGRP,1)%A, &
                    STSTM(IGRP,2)%A,NACOB,IEC,LROW)

  end if
end do

! Now to supergroups, i.e. strings of with given number of elecs in each GAspace

NSTFSMSPGP(:,1:NTSPGP) = 0
MXNSTR = -1
do ITP=1,NSTTYP
  ! Loop over supergroups of given type . i.e. strings
  ! with given occupation in each GAS space
  do IGRP=1,NSPGPFTP(ITP)
    IGRPABS = IGRP-1+IBSPGPFTP(ITP)
    call NSTPTP_GAS(NGAS,ISPGPFTP(1,IGRPABS),NSTSGP,NIRREP,NSTSO(ITP)%A,IGRP,MXNSTRFSG,NSMCLS,NSMCLSE,NSMCLSE1)

    !MXSMCLS = max(MXSMCLS,NSMCLS)
    !MXSMCLSE = max(MXSMCLSE,NSMCLSE)
    !MXSMCLSE1 = max(MXSMCLSE1,NSMCLSE1)

    MXNSTR = max(MXNSTR,MXNSTRFSG)
  end do

  do i=1,NSPGPFTP(ITP)
    NSTFSMSPGP(1:NIRREP,IBSPGPFTP(ITP)+i-1) = NSTSO(ITP)%A((i-1)*NIRREP+1:i*NIRREP)
  end do
  ! Corresponding offset array : Each supergroup is generated individually
  ! so each supergroup starts with offset 1 !
  call ZSPGPIB(NSTSO(ITP)%A,ISTSO(ITP)%A,NSPGPFTP(ITP),NIRREP)

# ifdef _DEBUGPRINT_
  write(u6,*) ' Number of strings per sym (row) and supergroup(column) for type = ',ITP
  call IWRTMA(NSTSO(ITP)%A,NIRREP,NSPGPFTP(ITP),NIRREP,NSPGPFTP(ITP))
  write(u6,'(A,3I6)') ' NSMCLS,NSMCLSE,NSMCLSE1=',NSMCLS,NSMCLSE,NSMCLSE1
  write(u6,*)
# endif

end do
! Number of electron in each AS for each supergroup
call ZNELFSPGP()

! Number of holes per supergroup
do IISPGP=1,NTSPGP
  NHOLE = 0
  do IGAS=1,NGAS
    if (IPHGAS(IGAS) == 2) NHOLE = NHOLE+NELFSPGP(IGAS,IISPGP)
  end do
  NHLFSPGP(IISPGP) = NHOLE
end do
#ifdef _DEBUGPRINT_
write(u6,*) ' Number of electrons in hole spaces per supergroup'
call IWRTMA(NHLFSPGP,1,NTSPGP,1,NTSPGP)
#endif
! Largest number of strings belonging to given supergroup
! Largest Occupation block for given supergroup and sym
MAX_STR_OC_BLK = -1
MAX_STR_SPGP = 0
do ISPGP=1,NTSPGP
  NSTR = sum(NSTFSMSPGP(1:NIRREP,ISPGP))
  MAX_STR_SPGP = max(MAX_STR_SPGP,NSTR)
  NEL = sum(NELFSPGP(1:NGAS,ISPGP))
  MAX_STR_OC_BLK = max(MAX_STR_OC_BLK,(NEL+4)*maxval(NSTFSMSPGP(1:NIRREP,ISPGP)))
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Largest number of strings of given supergroup        ',MAX_STR_SPGP
write(u6,*) ' Largest block of string occupations ',MAX_STR_OC_BLK
write(u6,*) ' Largest number of strings of given supergroup and sym',MXNSTR
#endif
!write(u6,'(A,3I6)') ' MXSMCLS,MXSMCLSE,MXSMCLSE1 = ',MXSMCLS,MXSMCLSE,MXSMCLSE1

! Possible occupation classes

ZERO_ARR(1) = 0
call OCCLS(2,NMXOCCLS,IOCLS,NACTEL,NGAS,IGSOCC(1,1),IGSOCC(1,2),0,ZERO_ARR,NOBPT)

! Maps creation/annihilation of given gas orb from given supergroup
! gives new supergroup.

SPGPCR(:) = 0
SPGPAN(:) = 0

do ISTTYP=1,NSTTYP
  if (NSPGPFTP(ISTTYP) > 0) then
    ! Creation map from this type
    IIEL = NELFTP(ISTTYP)
    ! Type of string with one elec more
    ISTTYPC = 0
    do JSTTYP=1,NSTTYP
      if ((mod(ISTTYP,2) == mod(JSTTYP,2)) .and. (NELFTP(JSTTYP) == IIEL+1)) ISTTYPC = JSTTYP
    end do
    !write(u6,*) ' ISTTYP and ISTTYPC ',ISTTYP,ISTTYPC
    if (ISTTYPC >= 1) then
      if (NSPGPFTP(ISTTYPC) > 0) &
        call SPGP_AC(NELFSPGP(1,1),NSPGPFTP(ISTTYP),NELFSPGP(1,1),NSPGPFTP(ISTTYPC),NGAS,MXPNGAS,2,SPGPCR,IBSPGPFTP(ISTTYP), &
                     IBSPGPFTP(ISTTYPC))
    end if

    ! Annihilation maps
    ISTTYPA = 0
    do JSTTYP=1,NSTTYP
      if ((mod(ISTTYP,2) == mod(JSTTYP,2)) .and. (NELFTP(JSTTYP) == IIEL-1)) ISTTYPA = JSTTYP
    end do
    !write(g6,*) 'ISTTYP, ISTTYPA', ISTTYP,ISTTYPA
    if (ISTTYPA >= 1) then
      if (NSPGPFTP(ISTTYPA) > 0) &
        call SPGP_AC(NELFSPGP(1,1),NSPGPFTP(ISTTYP),NELFSPGP(1,1),NSPGPFTP(ISTTYPA),NGAS,MXPNGAS,1,SPGPAN,IBSPGPFTP(ISTTYP), &
                     IBSPGPFTP(ISTTYPA))
    end if
  end if
end do

end subroutine STRINF_GAS
