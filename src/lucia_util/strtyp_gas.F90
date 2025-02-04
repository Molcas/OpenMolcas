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

subroutine STRTYP_GAS(IPRNT)
! Find groups of strings in each GA space
!
! Output : /GASSTR/
!
! Jeppe Olsen,  Oct 1994
! G. Li Manni, June 2024: Scale-up capability for single SD ROHF type calculations

use lucia_data, only: NGAS, NCISPC, NMXOCCLS, IGSOCC, IGSOCCX, IPHGAS, IPHGAS1
use lucia_data, only: MS2
use lucia_data, only: NGRP, NSTTP, NTSPGP, NSPGPFTP, IBSPGPFTP, NELFTP, NGPSTR, MXGSOC, MNGSOC, IBGPSTR, IGSFGP, ISPGPFTP, &
                      MNELFGP, MNGSOC, MXELFGP, MXGSOC, NELFGP, NELFSPGP, NSTFGP
use lucia_data, only: NACTEL
use lucia_data, only: MXTSOB_P, MXTSOB_H, NOBPT
use lucia_data, only: NOCTYP, ISTAC
use lucia_data, only: NSTTYP, NELEC
use lucia_data, only: MXPSTT
use Definitions, only: u6

implicit none
integer IPRNT
! Local scratch
integer IOCTYP(MXPSTT), IREOSPGP(MXPSTT), ISCR(MXPSTT)
integer IOCCLS(1), IBASSPC(1)
integer NTESTL, NTEST, IGAS, MINI, MAXI, ICISPC, IPHGASL, NPHGAS, NOCCLS, NAEL, NBEL, MAXSUB, IGRP, MXAL, MNAL, MXBL, MNBL, MXA1, &
        MXB1, MNA1, MNB1, MXA, MXB, MNA, MNB, MNAB, MXAB, IADD, JGRP, IEL, IITYPE, ITP, IOFF, NABEL, NSPGP_TOT, ITYP, NELEC_REF, &
        IDEL, NSPGP, IONE, NEL, IRED, NONEW, I_AM_OKAY, IOELMX, IBTYP, ISPGP, IDIM, JGAS, ISPGP_N, ISPGP_O
integer, external :: IBION_LUCIA

NTESTL = 0
NTEST = max(IPRNT,NTESTL)
! As input NCISPC GAS spaces IGSOCCX are given.
! Obtain space that cantains all these as special cases

!write(u6,*) ' NCISPC ',NCISPC
do IGAS=1,NGAS
  MINI = IGSOCCX(IGAS,1,1)
  MAXI = IGSOCCX(IGAS,2,1)
  !write(u6,*) ' MINI and MAXI for ISPC = 1 ',MINI,MAXI
  do ICISPC=2,NCISPC
    MINI = min(MINI,IGSOCCX(IGAS,1,ICISPC))
    MAXI = max(MAXI,IGSOCCX(IGAS,2,ICISPC))
    !write(u6,*) ' MINI and MAXI for ISPC =  ',ICISPC,MINI,MAXI
  end do
  IGSOCC(IGAS,1) = MINI
  IGSOCC(IGAS,2) = MAXI
end do

if (NTEST >= 5) then
  write(u6,*) ' Compound GAS space :'
  write(u6,*) ' ===================='
  write(u6,'(A)')
  write(u6,'(A)') '         Min. occ    Max. occ'
  write(u6,'(A)') '         ========    ========'
  do IGAS=1,NGAS
    write(u6,'(A,I2,3X,I3,9X,I3)') '   GAS',IGAS,IGSOCC(IGAS,1),IGSOCC(IGAS,2)
  end do
end if

! Find min and max number of elecs in each subspace

do IGAS=1,NGAS
  if (IGAS == 1) then
    MNGSOC(IGAS) = IGSOCC(IGAS,1)
    MXGSOC(IGAS) = IGSOCC(IGAS,2)
  else
    MXGSOC(IGAS) = IGSOCC(IGAS,2)-IGSOCC(IGAS-1,1)
    MNGSOC(IGAS) = max(0,IGSOCC(IGAS,1)-IGSOCC(IGAS-1,2))
  end if
end do

! Particle and hole spaces  :
! Hole spaces are always more than half occupied

IPHGASL = 0
NPHGAS = 0
do IGAS=1,NGAS
  !if (IUSE_PH == 1) then
  !  ! P/H separation for compound space
  !  if (MNGSOC(IGAS) > NOBPT(IGAS)) then
  !    IPHGAS(IGAS) = 2
  !    IPHGASL = IGAS
  !    NPHGAS = NPHGAS+1
  !  else
  !    IPHGAS(IGAS) = 1
  !  end if
  !  ! P/H separation for initial space
  !  if (IGAS == 1) then
  !    MIN_OC1 = IGSOCCX(IGAS,1,1)
  !  else
  !    MIN_OC1 = IGSOCCX(IGAS,1,1)-IGSOCCX(IGAS-1,2,1)
  !  end if
  !  if (MIN_OC1 > NOBPT(IGAS)) then
  !    IPHGAS1(IGAS) = 2
  !  else
  !    IPHGAS1(IGAS) = 1
  !  end if
  !else if (IUSE_PH == 0) then
  IPHGAS(IGAS) = 1
  IPHGAS1(IGAS) = 1
  !end if
end do
! Large number of particle and hole orbitals
MXTSOB_P = 0
MXTSOB_H = 0
do IGAS=1,NGAS
  if (IPHGAS1(IGAS) == 1) then
    MXTSOB_P = max(MXTSOB_P,NOBPT(IGAS))
  else
    MXTSOB_H = max(MXTSOB_H,NOBPT(IGAS))
  end if
end do
if (NTEST > 0) write(u6,*) ' MXTSOB_H, MXTSOB_P = ',MXTSOB_H,MXTSOB_P

!if (IUSE_PH == 1) then
!  IPHGAS(1) = 2
!  do I=1,100
!    write(u6,*) ' First space enforced to hole spaces'
!  end do
!end if

! In the following I assume that the hole spaces are the first NPGAS spaces
! (only used when calculating min number of electrons, so it can be modified easily)
! Min number of electrons in hole spaces
!if (NTEST > 0) write(u6,*) ' IPHGASL, NPHGAS ',IPHGASL,NPHGAS
!_Jesper if (IPHGASL /= NPHGAS) then
!_Jesper   write(u6,*) ' The hole spaces are not the first orbital spaces'
!_Jesper   stop ' The hole spaces are not the first orbital spaces'
!_Jesper end if
!MNHL = IGSOCC(IPHGASL,1)
!MNHL = 0
!do IGAS = 1, NGAS
!  if (IPHGAS(IGAS) == 2) MNHL = MNHL+MNGSOC(IGAS)
!end do

if (NTEST > 0) then
  write(u6,*) ' IPHGASL,NPHGAS:',IPHGASL,NPHGAS
  if (IPHGASL > 0) write(u6,*) 'IGSOCC(IPHGASL,1):',IGSOCC(IPHGASL,1)
end if

if (NTEST >= 5) then
  write(u6,*)
  write(u6,'(A)') ' Min and Max occupation in each GAS space:'
  write(u6,'(A)') ' ========================================='
  write(u6,*)
  do IGAS=1,NGAS
    write(u6,'(A,I2,4X,2I3)') '  GAS',IGAS,MNGSOC(IGAS),MXGSOC(IGAS)
  end do

  write(u6,*) ' Particle(1) or hole(2) spaces (for compound space)'
  call IWRTMA(IPHGAS,1,NGAS,1,NGAS)
  write(u6,*) ' Particle(1) or hole(2) spaces (for initial space)'
  call IWRTMA(IPHGAS1,1,NGAS,1,NGAS)
end if

! Occupation classes corresponding to largest CI space

IBASSPC(1) = 0
call OCCLS(1,NOCCLS,IOCCLS,NACTEL,NGAS,IGSOCC(1,1),IGSOCC(1,2),0,IBASSPC,NOBPT)
NMXOCCLS = NOCCLS

! Split into alpha and beta parts

! Number of alpha. and beta electrons

NAEL = (MS2+NACTEL)/2
NBEL = (NACTEL-MS2)/2

if (NAEL+NBEL /= NACTEL) then
  write(u6,*) '  MS2 NACTEL NAEL NBEL'
  write(u6,'(5I4)') MS2,NACTEL,NAEL,NBEL
  write(u6,*) ' STOP : NUMBER OF ELECTRONS AND MULTIPLICITY INCONSISTENT'
  !stop ' NUMBER OF ELECTRONS INCONSISTENT WITH MULTIPLICITY'
  call SYSABENDMSG('lucia_util/strtyp_gas','Internal error','')
end if

if (NTEST >= 5) then
  write(u6,*) '  MS2 NACTEL NAEL NBEL'
  write(u6,'(5I6)') MS2,NACTEL,NAEL,NBEL
end if

if (NAEL+NBEL /= NACTEL) then
  write(u6,*) '  MS2 NACTEL NAEL NBEL'
  write(u6,'(5I4)') MS2,NACTEL,NAEL,NBEL
  write(u6,*) ' STOP : NUMBER OF ELECTRONS AND MULTIPLICITY INCONSISTENT'
  !stop ' NUMBER OF ELECTRONS INCONSISTENT WITH MULTIPLICITY'
  call SYSABENDMSG('lucia_util/strtyp_gas','Internal error','')
end if

! Number of electrons to be subtracted or added

MAXSUB = 2
! electrons are only added for systems that atleast have halffilled shells
IGRP = 0
MXAL = NAEL
MNAL = NAEL
MXBL = NBEL
MNBL = NBEL
do IGAS=1,NGAS
  ! occupation constraints 1
  MXA1 = min(MXGSOC(IGAS),NOBPT(IGAS),MXAL)
  MXB1 = min(MXGSOC(IGAS),NOBPT(IGAS),MXBL)
  MNA1 = max(0,MNGSOC(IGAS)-MXA1)
  MNB1 = max(0,MNGSOC(IGAS)-MXB1)

  ! Additional checks can be made here
  MXA = MXA1
  MXB = MXB1
  MNA = MNA1
  MNB = MNB1

  MXAL = MXAL-MNA
  MNAL = max(0,MNAL-MXA)
  MXBL = MXBL-MNB
  MNBL = max(0,MNBL-MXB)

  if (NTEST >= 100) then
    write(u6,*) ' Occupation numbers for IGAS = ',IGAS
    write(u6,*) ' MXAL MNAL MXBL MNBL ',MXAL,MNAL,MXBL,MNBL
    write(u6,*) ' MXA MNA MXB MNB ',MXA,MNA,MXB,MNB
  end if

  MNAB = min(MNA,MNB)
  MXAB = max(MXA,MXB)

  ! Additional holes only allowed in particle spaces
  if (IPHGAS(IGAS) == 1) then
    MNAB = max(0,MNAB-MAXSUB)
  else if (IPHGAS(IGAS) == 2) then
    MNAB = MNAB
  end if
  ! For coupled cluster- could be refined ...
  MNAB = 0
  if (IPHGAS(IGAS) == 2) MXAB = min(MXAB+2,NOBPT(IGAS))

  if (NTEST >= 100) write(u6,*) ' MNAB,MXAB',MNAB,MXAB
  NGPSTR(IGAS) = MXAB-MNAB+1
  if ((Nactel == MS2) .and. (Nactel > 2) .and. (Nactel == NOBPT(2)) .and. (IGAS == 2)) then
    NGPSTR(IGAS) = 4 ! Either EMPTY, (FULL-2), (FULL -1), FULL
    MNAB = NAEL-2
  end if
  IBGPSTR(IGAS) = IGRP+1
  MNELFGP(IGAS) = MNAB
  MXELFGP(IGAS) = MXAB

  IADD = 0
  do JGRP=IGRP+1,IGRP+NGPSTR(IGAS)
    if (JGRP > MXPSTT) then
      write(u6,*) ' Too many string groups'
      write(u6,*) ' Current limit ',MXPSTT
      write(u6,*) ' STOP : GASSTR, Too many string groups'
      !stop' GASSTR, Too many string groups'
      call SYSABENDMSG('lucia_util/gasstr','Internal error','')
    end if

    if ((Nactel == MS2) .and. (Nactel > 2) .and. (Nactel == NOBPT(2)) .and. (IGAS == 2) .and. (JGRP == 2)) then
      IEL = 0
    else
      IADD = IADD+1
      IEL = MNAB-1+IADD
    end if
    NELFGP(JGRP) = IEL
    IGSFGP(JGRP) = IGAS
    NSTFGP(JGRP) = IBION_LUCIA(NOBPT(IGAS),IEL)
  end do
  IGRP = IGRP+NGPSTR(IGAS)
end do
NGRP = IGRP

if (NTEST >= 5) then
  write(u6,*)
  write(u6,'(A)') ' Information about Groups of strings'
  write(u6,'(A)') ' ==================================='
  write(u6,*)
  write(u6,*) '     GAS  MNEL  MXEL IBGRP  NGRP'
  write(u6,*) '    ============================'
  do IGAS=1,NGAS
    write(u6,'(5(2X,I4))') IGAS,MNELFGP(IGAS),MXELFGP(IGAS),IBGPSTR(IGAS),NGPSTR(IGAS)
  end do
  write(u6,'(A,I3)') ' Total number of groups generated ',NGRP

  write(u6,'(A)') ' Information about each string group'
  write(u6,'(A)') ' ==================================='
  write(u6,*)
  IITYPE = 0
  write(u6,'(A)') ' GROUP  GAS   NEL      NSTR'
  write(u6,'(A)') ' =========================='
  do IGRP=1,NGRP
    IITYPE = IITYPE+1
    write(u6,'(3(2X,I4),2X,I8)') IITYPE,IGSFGP(IGRP),NELFGP(IGRP),NSTFGP(IGRP)
  end do
end if

! Creation-annihilation connections between groups

do IGRP=1,NGRP
  ISTAC(IGRP,1) = 0
  ISTAC(IGRP,2) = 0
  do JGRP=1,NGRP
    if ((IGSFGP(IGRP) == IGSFGP(JGRP)) .and. (NELFGP(IGRP) == NELFGP(JGRP)-1)) ISTAC(IGRP,2) = JGRP
    if ((IGSFGP(IGRP) == IGSFGP(JGRP)) .and. (NELFGP(IGRP) == NELFGP(JGRP)+1)) ISTAC(IGRP,1) = JGRP
  end do
end do

if (NTEST >= 5) then
  write(u6,*)
  write(u6,*) ' ======================================'
  write(u6,*) ' Annihilation / Creation connections'
  write(u6,*) ' ======================================'
  write(u6,*)
  call IWRTMA(ISTAC,NGRP,2,MXPSTT,2)
end if

! Construct number of type ( combinations of groups ) with nael and nbel strings

! Type 1 : NAEL electrons
!      2 : NBEL ELECTRONS
!      3 : NAEL -1 ELECTRONS
!      4 : NBEL -1 ELECTRONS
!      5 : NAEL -2 ELECTRONS
!      6 : NBEL -2 ELECTRONS

NSTTYP = 6
NSTTP = 6

!if (IUSE_PH == 1) then
!  ! allow N+1,N+2 resolution string
!  NSTTYP = 10
!  NSTTP = 10
!end if
! alpha
NELEC(1) = NAEL
NELFTP(1) = NAEL
! beta
NELEC(2) = NBEL
NELFTP(2) = NBEL
! alpha -1
NELEC(3) = NAEL-1
NELFTP(3) = NAEL-1
! beta  -1
NELEC(4) = NBEL-1
NELFTP(4) = NBEL-1
! alpha -2
NELEC(5) = NAEL-2
NELFTP(5) = NAEL-2
! beta  -2
NELEC(6) = NBEL-2
NELFTP(6) = NBEL-2

!if (IUSE_PH == 1) then
!  ! Alpha + 1
!  NELEC(7) = NAEL+1
!  NELFTP(7) = NAEL+1
!  ! beta  + 1
!  NELEC(8) = NBEL+1
!  NELFTP(8) = NBEL+1
!  ! Alpha + 2
!  NELEC(9) = NAEL+2
!  NELFTP(9) = NAEL+2
!  ! beta  + 2
!  NELEC(10) = NBEL+2
!  NELFTP(10) = NBEL+2
!end if
! Can easily be extended to relativistic case !!
do ITP=1,NSTTYP
  NOCTYP(ITP) = 0
  NSPGPFTP(ITP) = 0
end do

! Loop over types, i.e.  given number of electrons

IOFF = 1
NABEL = NAEL+NBEL
NSPGP_TOT = 0
do ITYP=1,NSTTYP
  ! Number of electrons in reference space ( alpha or beta )
  if (mod(ITYP,2) == 1) then
    ! alpha type
    NELEC_REF = NELEC(1)
  else
    ! beta type
    NELEC_REF = NELEC(2)
  end if
  ! If we are studying beta type, and number of alpha and beta
  ! electrons are identical, just refer to alpha
  if ((NAEL == NBEL) .and. (mod(ITYP,2) == 0)) then
    IBSPGPFTP(ITYP) = IBSPGPFTP(ITYP-1)
    NOCTYP(ITYP) = NOCTYP(ITYP-1)
    NSPGPFTP(ITYP) = NSPGPFTP(ITYP-1)
  else
    ! Number of electrons removed compared to reference
    IDEL = NELEC(ITYP)-NELEC_REF
    !write(u6,*) '  GASSPC : ITYP IDEL ',ITYP,IDEL
    ! Initial type of strings, relative to offset for given group
    do IGAS=1,NGAS
      IOCTYP(IGAS) = 1
    end do
    NSPGP = 0
    IBSPGPFTP(ITYP) = IOFF
    if (NELEC(ITYP) < 0) then
      NOCTYP(ITYP) = 0
      NSPGPFTP(ITYP) = 0
      goto 2000
    end if
    ! Number of electrons in present type
    ! Loop over  SUPER GROUPS with current nomenclature!
    IONE = 1
1000 continue
    ! Number of electrons in present supergroup
    NEL = 0
    do IGAS=1,NGAS
      NEL = NEL+NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
    end do

    if (NEL > NELEC(ITYP)) then
      ! If the number of electrons is to large find next number that
      ! can be correct.
      ! The following uses that within a given GAS space
      ! the number of elecs increases as the type number increases

      ! First integer  that can be reduced
      IRED = 0
      do IGAS=1,NGAS
        if (IOCTYP(IGAS) /= 1) then
          IRED = IGAS
          goto 888
        end if
      end do
888   continue
      if (IRED == NGAS) then
        NONEW = 1
      else if (IRED < NGAS) then
        IOCTYP(IRED) = 1
        ! Increase remanining part
        call NXTNUM2(IOCTYP(IRED+1),NGAS-IRED,IONE,NGPSTR(IRED+1),NONEW)
      end if
      goto 2803
    end if

    if (NEL == NELEC(ITYP)) then
      ! test 1 has been passed, check additional occupation constraints

      I_AM_OKAY = 1
      ! Number of extra holes in hole spaces
      !E if (IUSE_PH == 1) then
      !E   IDELP = 0
      !E   IDELM = 0
      !E   do IGAS=1,NGAS
      !E     if (IPHGAS(IGAS) == 2) then
      !E       NELH = NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
      !E       if (NELH < MNGSOC(IGAS)) IDELM = IDELM+MNGSOC(IGAS)-NELH
      !E       if (NELH > MXGSOC(IGAS)) IDELP = IDELP+NELH-MXGSOC(IGAS)
      !E     end if
      !E   end do
      !E   if ((IDELM > 0) .or. (IDELP > MAX(0,IDEL))) then
      !E     I_AM_OKAY = 0
      !E     write(u6,*) ' P/H rejected supergroup'
      !E     call IWRTMA(IOCTYP,1,NGAS,1,NGAS)
      !E     write(u6,*) ' IDELM, IDELP ',IDELM,IDELP
      !E   end if
      !E end if

      ! Check from above

      do IGAS=NGAS,1,-1
        ! Number of electrons when all electrons of AS IGAS have been added
        if (IGAS == NGAS) then
          IEL = max(NABEL,NABEL+IDEL)
        else
          IEL = IEL-NELFGP(IOCTYP(IGAS+1)+IBGPSTR(IGAS+1)-1)
          if (IEL < max(IGSOCC(IGAS,1),IGSOCC(IGAS,1)+IDEL)) I_AM_OKAY = 0
        end if
      end do

      ! Check from below

      IEL = 0
      IOELMX = 0
      do IGAS=1,NGAS
        IEL = IEL+NELFGP(IOCTYP(IGAS)+IBGPSTR(IGAS)-1)
        IOELMX = IOELMX+NOBPT(IGAS)
        if (IEL+IOELMX < min(IGSOCC(IGAS,1),IGSOCC(IGAS,1)+IDEL)) I_AM_OKAY = 0
      end do

      if (I_AM_OKAY == 1) then
        ! passed !!!
        NSPGP = NSPGP+1
        ! Copy supergroup to ISPGPFTP with absolute group numbers
        do IGAS=1,NGAS
          ISPGPFTP(IGAS,IOFF-1+NSPGP) = IOCTYP(IGAS)+IBGPSTR(IGAS)-1
        end do
      end if

    end if
    ! Next type of strings
    IONE = 1
    call NXTNUM2(IOCTYP,NGAS,IONE,NGPSTR,NONEW)
2803 continue
    if (NONEW == 0) goto 1000
    ! End of loop over possible supergroups, save information about current type
    IOFF = IOFF+NSPGP
    NOCTYP(ITYP) = NSPGP
    NSPGPFTP(ITYP) = NSPGP
    NSPGP_TOT = NSPGP_TOT+NSPGP
  end if
2000 continue
end do
NTSPGP = NSPGP_TOT

if (NSPGP_TOT > MXPSTT) then
  write(u6,*) ' Too many super groups = ',NSPGP_TOT
  write(u6,*) ' Increase MXPSTT to this value'
  write(u6,*) ' See you later'
  write(u6,*)
  write(u6,*) ' STOP Increase MXPSTT'
  !stop ' Increase MXPSTT'
  call SYSABENDMSG('lucia_util/strtyp_gas','Internal error','')
end if

! Reorder supergroups according to dimension

do ITYP=1,NSTTP
  IBTYP = IBSPGPFTP(ITYP)
  NSPGP = NSPGPFTP(ITYP)
  ! Dimension of supergroups
  do ISPGP=1,NSPGP
    IDIM = 1
    do JGAS=1,NGAS
      IDIM = IDIM*NSTFGP(ISPGPFTP(JGAS,ISPGP+IBTYP-1))
    end do
    IOCTYP(ISPGP) = IDIM
  end do
  ! Reorder
  !    ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
  call ORDINT(IOCTYP,ISCR,NSPGP,IREOSPGP,NTEST)
  !write(u6,*) ' IREO array'
  !call IWRTMA(IREOSPGP,1,NSPGP,1,NSPGP)
  ! And reorder the definition of supergroups
  do ISPGP=1,NSPGP
    call ICOPVE(ISPGPFTP(1,ISPGP+IBTYP-1),NELFSPGP(1,ISPGP),NGAS)
  end do
  do ISPGP_N=1,NSPGP
    ISPGP_O = IREOSPGP(ISPGP_N)
    call ICOPVE(NELFSPGP(1,ISPGP_O),ISPGPFTP(1,ISPGP_N+IBTYP-1),NGAS)
  end do

end do
! End of loop over types

if (NTEST >= 2) then
  write(u6,*) ' Total number of super groups ',NTSPGP
  write(u6,*) ' Number of alpha supergroups  ',NSPGPFTP(1)
  write(u6,*) ' Number of beta  supergroups  ',NSPGPFTP(2)
  write(u6,*)
  write(u6,*)
end if

if (NTEST >= 5) then
  write(u6,*) ' Information about types of strings'
  write(u6,*) ' =================================='
  write(u6,*)
  do ITYP=1,NSTTYP
    write(u6,*)
    write(u6,*) '      Type : ',ITYP
    write(u6,*) '      ==============='
    write(u6,*) '      Number of electrons  ',NELFTP(ITYP)
    write(u6,*) '      Number of super groups ',NSPGPFTP(ITYP)
    write(u6,*) '      Supergroups'
    do ISPGP=1,NSPGPFTP(ITYP)
      IOFF = IBSPGPFTP(ITYP)
      call IWRTMA(ISPGPFTP(1,IOFF-1+ISPGP),1,NGAS,1,NGAS)
    end do
  end do
end if

end subroutine STRTYP_GAS
