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

subroutine STRINF_GAS(IPRNT)
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

use stdalloc, only: mma_allocate, mma_deallocate
use strbas, only: ZMAT, NSTSGP, ISTSGP, STREO, OCSTR, STSTM, NSTSO, ISTSO, IOCLS, SPGPAN, SPGPCR
use lucia_data, only: NGAS, IGSOCC, IPHGAS, NMXOCCLS
! modification Jeppe + Giovanni + Dongxia.
! G. Li Manni, June 2024: Scale-up capability for single SD ROHF type calculations
use distsym, only: INGRP_VAL, ISMDFGP, ISMSCR, NACTSYM
use lucia_data, only: MS2
use lucia_data, only: NGRP, NTSPGP, MXNSTR, MXSMCLS, MXSMCLSE, MXSMCLSE1, MAX_STR_OC_BLK, MAX_STR_SPGP, MINMAX_SM_GP, IBSPGPFTP, &
                      IGSFGP, ISPGPFTP, ISTFSMGP, NELFGP, NELFSPGP, NELFTP, NHLFSPGP, NSPGPFTP, NSTFGP, NSTFSMGP, NSTFSMSPGP
use lucia_data, only: NACTEL
use lucia_data, only: NACOB, NORB1, NORB2, NORB3, NOBPT
use lucia_data, only: ISTAC
use lucia_data, only: NSTTYP
use lucia_data, only: MXPNSMST, MXPNGAS
use csm_data, only: NSMST
use Definitions, only: u6

implicit none
integer IPRNT
integer ZERO_ARR(1), IDUM(1)
integer, allocatable :: FREEL(:)
!  A bit of scratch
integer, external :: IELSUM
integer LAC, NTEST, IEL, IGRP, NACOB_EFFECTIVE, MAXSCR, IGAS, MNRS1X, MXRS1X, MNRS3X, MXRS3X, IOCTYPX, IGP, MX, ISM, MN, NGSOBP, &
        IGSOB, NSTINI, IEC, LROW, IZERO, JGRP, ITP, IGRPABS, NSMCLS, NSMCLSE, NSMCLSE1, IISPGP, NHOLE, ISPGP, NSTR, NEL, ISTSM, &
        ISTTYP, IIEL, ISTTYPC, JSTTYP, ISTTYPA, MXNSTRFSG

! Some dummy initializtions
LAC = 0 ! jwk-cleanup

NTEST = 0
NTEST = max(NTEST,IPRNT)

! 2 : Number of classes per string type and mappings between
!     string types (/STINF/)

if ((NActEl /= MS2) .or. (NActEl /= NACOB)) call ZSTINF_GAS(IPRNT)

! 3 : Static memory for string information

call MEMSTR_GAS()

!  4 : Info about group of strings

! First free address

! Find maximum needed length of scratch.
!   MAXSCR = 2*NACOB+(IEL+1)(NACOB+1)
!   with IEL = MAX(NELFGP(IGRP=1,NGRP)

IEL = 0
do IGRP=1,NGRP
  IEL = max(IEL,NELFGP(IGRP))
end do
NACOB_EFFECTIVE = NACOB
if (NACOB == 0) NACOB_EFFECTIVE = 1
MAXSCR = 2*NACOB_EFFECTIVE+(IEL+1)*(NACOB_EFFECTIVE+1)+NSMST
call mma_allocate(FREEL,MAXSCR,Label='FREEL')
do IGRP=1,NGRP
  ! A gas group can be considered as a RAS group with 0 electrons in RAS1, RAS3!
  IGAS = IGSFGP(IGRP)
  if (IGAS == 1) then
    NORB1 = 0
  else
    NORB1 = IELSUM(NOBPT,IGAS-1)
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
  call WEIGHT_LUCIA(Zmat(IGRP)%I,IEL,NORB1,NORB2,NORB3,MNRS1X,MXRS1X,MNRS3X,MXRS3X,FREEL,IPRNT)
  ! Number of strings per symmetry in a given group
  call NSTRSO_GAS(IEL,NORB1,NORB2,NORB3,MNRS1X,MXRS1X,MNRS3X,MXRS3X,FREEL,NACOB,NSTSGP(1)%I,ISTSGP(1)%I,IOCTYPX,NSMST,IGRP,IPRNT)
  ! Construct the strings ordered by symmetry
  call GENSTR_GAS(IEL,MNRS1X,MXRS1X,MNRS3X,MXRS3X,ISTSGP(1)%I,IGRP,IOCTYPX,NSMST,Zmat(IGRP)%I,FREEL,STREO(IGRP)%I,OCSTR(IGRP)%I, &
                  FREEL(1+IOCTYPX*NSMST),IGRP,IPRNT)

  call ICOPVE2(NSTSGP(1)%I,1+(IGRP-1)*NSMST,NSMST,NSTFSMGP(1,IGRP))
  call ICOPVE2(ISTSGP(1)%I,1+(IGRP-1)*NSMST,NSMST,ISTFSMGP(1,IGRP))
end do
call mma_deallocate(FREEL)

INGRP_VAL = NGRP
call mma_allocate(ISMDFGP,NSMST*NGRP,Label='ISMDFGP')
call mma_allocate(NACTSYM,NGRP,Label='NACTSYM')
call mma_allocate(ISMSCR,NGRP,Label='ISMSCR')
call SMDFGP_GEN(NGRP,NSMST,MXPNSMST,NSTFSMGP,NACTSYM,ISMDFGP)

if (NTEST >= 10) then
  write(u6,*) 'NGRP',NGRP
  write(u6,*) 'NSMST*NGRP',NSMST*NGRP
  write(u6,*) ' Number of strings per group and symmetry'
  call IWRTMA10(NSTSGP(1)%I,NSMST,NGRP,NSMST,NGRP)
  write(u6,*) ' Number of strings per group and symmetry(2)'
  call IWRTMA10(NSTFSMGP,NSMST,NGRP,MXPNSMST,NGRP)
end if

! Min and max of sym for each group

do IGP=1,NGRP
  MX = 1
  do ISM=1,NSMST
    if (NSTFSMGP(ISM,IGP) > 0) MX = ISM
  end do

  MN = NSMST
  do ISM=NSMST,1,-1
    if (NSTFSMGP(ISM,IGP) > 0) MN = ISM
  end do

  MINMAX_SM_GP(1,IGP) = MN
  MINMAX_SM_GP(2,IGP) = MX

end do
if (NTEST > 5) then
  write(u6,*) ' MINMAX array for sym of groups'
  write(u6,*) ' =============================='
  call IWRTMA(MINMAX_SM_GP,2,NGRP,2,NGRP)
end if

! 4.5 : Creation/Annihilation mappings between different
!       types of strings

do IGRP=1,NGRP

  IGAS = IGSFGP(IGRP)
  NGSOBP = NOBPT(IGAS)
  ! First orbital in GAS spacce
  IGSOB = IELSUM(NOBPT,IGAS-1)+1
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
    IZERO = 0
    call ISETVC(STSTM(IGRP,1)%I,IZERO,LROW*NSTINI)
    call ISETVC(STSTM(IGRP,2)%I,IZERO,LROW*NSTINI)
  end if

  if (ISTAC(IGRP,2) /= 0) then
    JGRP = ISTAC(IGRP,2)
    call CRESTR_GAS(OCSTR(IGRP)%I,NSTFGP(IGRP),NSTFGP(JGRP),IEL,NGSOBP,IGSOB,Zmat(JGRP)%I,STREO(JGRP)%I,0,IDUM,IDUM, &
                    STSTM(IGRP,1)%I,STSTM(IGRP,2)%I,NACOB,IPRNT)

  end if
  if (ISTAC(IGRP,1) /= 0) then
    JGRP = ISTAC(IGRP,1)
    call ANNSTR_GAS(OCSTR(IGRP)%I,NSTFGP(IGRP),NSTFGP(JGRP),IEL,NGSOBP,IGSOB,Zmat(JGRP)%I,STREO(JGRP)%I,0,IDUM,IDUM, &
                    STSTM(IGRP,1)%I,STSTM(IGRP,2)%I,NACOB,IEC,LROW,IPRNT)

  end if
end do

! Now to supergroups, i.e. strings of with given number of elecs in each GAspace

call ISETVC(NSTFSMSPGP,0,MXPNSMST*NTSPGP)
MXNSTR = -1
do ITP=1,NSTTYP
  ! Loop over supergroups of given type . i.e. strings
  ! with given occupation in each GAS space
  do IGRP=1,NSPGPFTP(ITP)
    IGRPABS = IGRP-1+IBSPGPFTP(ITP)
    call NSTPTP_GAS(NGAS,ISPGPFTP(1,IGRPABS),NSTSGP(1)%I,NSMST,NSTSO(ITP)%I,IGRP,MXNSTRFSG,NSMCLS,NSMCLSE,NSMCLSE1)

    MXSMCLS = max(MXSMCLS,NSMCLS)
    MXSMCLSE = max(MXSMCLSE,NSMCLSE)
    MXSMCLSE1 = max(MXSMCLSE1,NSMCLSE1)

    MXNSTR = max(MXNSTR,MXNSTRFSG)
  end do

  call ICOPMT(NSTSO(ITP)%I,NSMST,NSPGPFTP(ITP),NSTFSMSPGP(1,IBSPGPFTP(ITP)),MXPNSMST,NSPGPFTP(ITP))
  ! Corresponding offset array : Each supergroup is generated individually
  ! so each supergroup starts with offset 1 !
  call ZSPGPIB(NSTSO(ITP)%I,ISTSO(ITP)%I,NSPGPFTP(ITP),NSMST)

  if (NTEST >= 5) then
    write(u6,*) ' Number of strings per sym (row) and supergroup(column) for type = ',ITP
    call IWRTMA(NSTSO(ITP)%I,NSMST,NSPGPFTP(ITP),NSMST,NSPGPFTP(ITP))
    write(u6,'(A,3I6)') ' NSMCLS,NSMCLSE,NSMCLSE1=',NSMCLS,NSMCLSE,NSMCLSE1
    write(u6,*)
  end if

end do
! Number of electron in each AS for each supergroup
call ZNELFSPGP(IPRNT)

! Number of holes per supergroup
do IISPGP=1,NTSPGP
  NHOLE = 0
  do IGAS=1,NGAS
    if (IPHGAS(IGAS) == 2) NHOLE = NHOLE+NELFSPGP(IGAS,IISPGP)
  end do
  NHLFSPGP(IISPGP) = NHOLE
end do
if (NTEST >= 10) then
  write(u6,*) ' Number of electrons in hole spaces per supergroup'
  call IWRTMA(NHLFSPGP,1,NTSPGP,1,NTSPGP)
end if
! Largest number of strings belonging to given supergroup
! Largest Occupation block for given supergroup and sym
MAX_STR_OC_BLK = -1
MAX_STR_SPGP = 0
do ISPGP=1,NTSPGP
  NSTR = IELSUM(NSTFSMSPGP(1,ISPGP),NSMST)
  MAX_STR_SPGP = max(MAX_STR_SPGP,NSTR)
  NEL = IELSUM(NELFSPGP(1,ISPGP),NGAS)
  do ISTSM=1,NSMST
    MAX_STR_OC_BLK = max(MAX_STR_OC_BLK,(NEL+4)*NSTFSMSPGP(ISTSM,ISPGP))
    !MOD = max(MAX_STR_OC_BLK,NEL*NSTFSMSPGP(ISTSM,ISPGP))
  end do
end do

if (NTEST >= 2) then
  write(u6,*) ' Largest number of strings of given supergroup        ',MAX_STR_SPGP
  write(u6,*) ' Largest block of string occupations ',MAX_STR_OC_BLK

  write(u6,*) ' Largest number of strings of given supergroup and sym',MXNSTR
end if
!write(u6,'(A,3I6)') ' MXSMCLS,MXSMCLSE,MXSMCLSE1 = ',MXSMCLS,MXSMCLSE,MXSMCLSE1

! Possible occupation classes

ZERO_ARR(1) = 0
call OCCLS(2,NMXOCCLS,IOCLS,NACTEL,NGAS,IGSOCC(1,1),IGSOCC(1,2),0,ZERO_ARR,NOBPT)

! Maps creation/annihilation of given gas orb from given supergroup
! gives new supergroup.

IZERO = 0
call ISETVC(SPGPCR,IZERO,NGAS*NTSPGP)
call ISETVC(SPGPAN,IZERO,NGAS*NTSPGP)

do ISTTYP=1,NSTTYP
  ! Creation map from this type
  IIEL = NELFTP(ISTTYP)
  ! Type of string with one elec more
  ISTTYPC = 0
  do JSTTYP=1,NSTTYP
    if ((mod(ISTTYP,2) == mod(JSTTYP,2)) .and. (NELFTP(JSTTYP) == IIEL+1)) ISTTYPC = JSTTYP
  end do
  !write(u6,*) ' ISTTYP and ISTTYPC ',ISTTYP,ISTTYPC
  if (NSPGPFTP(ISTTYP) > 0) then

    if ((ISTTYPC >= 1) .and. (NSPGPFTP(ISTTYPC) > 0)) &
      call SPGP_AC(NELFSPGP(1,1),NSPGPFTP(ISTTYP),NELFSPGP(1,1),NSPGPFTP(ISTTYPC),NGAS,MXPNGAS,2,SPGPCR,IBSPGPFTP(ISTTYP), &
                   IBSPGPFTP(ISTTYPC))
    ! Annihilation maps
    ISTTYPA = 0
    do JSTTYP=1,NSTTYP
      if ((mod(ISTTYP,2) == mod(JSTTYP,2)) .and. (NELFTP(JSTTYP) == IIEL-1)) ISTTYPA = JSTTYP
    end do
    !write(g6,*) 'ISTTYP, ISTTYPA', ISTTYP,ISTTYPA
    if ((ISTTYPA >= 1) .and. (NSPGPFTP(ISTTYPA) > 0)) &
      call SPGP_AC(NELFSPGP(1,1),NSPGPFTP(ISTTYP),NELFSPGP(1,1),NSPGPFTP(ISTTYPA),NGAS,MXPNGAS,1,SPGPAN,IBSPGPFTP(ISTTYP), &
                   IBSPGPFTP(ISTTYPA))
  end if
end do

end subroutine STRINF_GAS
