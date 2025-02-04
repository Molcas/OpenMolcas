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
! Copyright (C) 1991,1994-1996, Jeppe Olsen                            *
!***********************************************************************

subroutine ADSTN_GAS(OFFI,IOBSM,IOBTP,ISPGP,ISPGPSM,ISPGPTP,I1,XI1S,NKSTR,SCLFAC)
! Obtain mappings
! a+IORB !KSTR> = +/-!ISTR> for orbitals of symmetry IOBSM and type IOBTP
! and I strings belonging to supergroup ISPGP wih symmetry ISPGPSM
! and type ISPGPTP(=1=>alpha,=2=>beta)
!
! The results are given in the form
! I1(KSTR,IORB) =  ISTR if A+IORB !KSTR> = +/-!ISTR>
! (numbering relative to TS start)
! Above +/- is stored in XI1S
!
! Jeppe Olsen, Winter of 1991
!              January 1994 : modified to allow for several orbitals
!              August 95    : GAS version
!              October 96   : Improved version

use strbas, only: NSTSGP, ISTSGP, STSTM
use Constants, only: Zero
use lucia_data, only: NGAS
use lucia_data, only: IBSPGPFTP, ISPGPFTP, NELFGP, NSTFSMSPGP
use lucia_data, only: IOBPTS, NOBPT, NOBPTS
use lucia_data, only: MXPNGAS, MXPNSMST
use csm_data, only: NSMST
use Definitions, only: wp, u6

implicit none
real*8 :: OFFI(*)
integer IOBSM, IOBTP, ISPGP, ISPGPSM, ISPGPTP, NKSTR
real*8 SCLFAC
! ======
! Output
! ======
integer I1(*)
real*8 XI1S(*)
! Local scratch
integer NELFGS(MXPNGAS), ISMFGS(MXPNGAS), ITPFGS(MXPNGAS)
integer maxval(MXPNGAS), minval(MXPNGAS)
integer NNSTSGP(MXPNSMST,MXPNGAS)
integer IISTSGP(MXPNSMST,MXPNGAS)
integer IACIST(MXPNSMST), NACIST(MXPNSMST)
integer, parameter :: MXLNGAS = 20
integer, external :: IELSUM
integer NTEST, ISPGRPABS, KSM, KSPGRPABS, NORBTS, IZERO, IBORBSP, IBORBSPS, NGASL, IGAS, NELB, NACGSOB, ISMST, IFIRST, NSTRINT, &
        NONEW, ISTSMM1, JSTSMM1, ISMGSN, NSTRII, IOFF, MULT, KACGRP, KFIRST, KSTRBS, NSTRIK, ISAVE, IACSM, IBSTRINI, NSTB, NSTA, &
        NIAC, IIAC, NKAC, IKAC, NKSD, IORB, IORBR

! Will be stored as an matrix of dimension
! (NKSTR,*), Where NKSTR is the number of K-strings of
! correct symmetry . Nk is provided by this routine.

if (NGAS > MXLNGAS) then
  write(u6,*) ' Ad hoc programming in ADSTN (IOFFI)'
  write(u6,*) ' Must be changed - or redimensioned'
  !stop 'ADST : IOFFI problem'
  call SYSABENDMSG('lucia_util/adstn_gas','Internal error','')
end if

NTEST = 0
if (NTEST >= 100) then
  write(u6,*)
  write(u6,*) ' ===================='
  write(u6,*) ' ADSTN_GAS in service'
  write(u6,*) ' ===================='
  write(u6,*)
  write(u6,*) '  IOBTP IOBSM : ',IOBTP,IOBSM
  write(u6,*) '  ISPGP ISPGPSM ISPGPTP :  ',ISPGP,ISPGPSM,ISPGPTP
end if

!if (SCLFAC /= One) then
!  write(u6,*) ' Problemo : ADSTN_GAS'
!  write(u6,*) ' SCLFAC /= 1'
!end if

! Supergroup and symmetry of K strings

ISPGRPABS = IBSPGPFTP(ISPGPTP)-1+ISPGP
call NEWTYP(ISPGRPABS,1,IOBTP,KSPGRPABS)
call SYMCOM(2,IOBSM,KSM,ISPGPSM)
NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
if (NTEST >= 200) write(u6,*) ' KSM, KSPGPRABS, NKSTR : ',KSM,KSPGRPABS,NKSTR
if (NKSTR == 0) goto 9999

NORBTS = NOBPTS(IOBTP,IOBSM)
call SETVEC(XI1S,Zero,NORBTS*NKSTR)
IZERO = 0
call ISETVC(I1,IZERO,NORBTS*NKSTR)

! First orbital of given GASSpace
IBORBSP = IELSUM(NOBPT,IOBTP-1)+1
! First orbital of fiven GASSPace and Symmetry
IBORBSPS = IOBPTS(IOBTP,IOBSM)

! Information about I strings
! ===========================

! structure of group of strings defining I strings
NGASL = 1
do IGAS=1,NGAS
  ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPABS)
  NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
  if (NELFGS(IGAS) > 0) NGASL = IGAS
end do
! Number of electrons before active type
NELB = 0
do IGAS=1,IOBTP-1
  NELB = NELB+NELFGS(IGAS)
end do
! Number of electrons in active space
NACGSOB = NOBPT(IOBTP)

! Number of strings per symmetry for each symmetry
do IGAS=1,NGAS
  call ICOPVE2(NSTSGP(1)%I,(ITPFGS(IGAS)-1)*NSMST+1,NSMST,NNSTSGP(1,IGAS))
end do
! Offset and dimension for active group in I strings
call ICOPVE2(ISTSGP(1)%I,(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,IACIST)
call ICOPVE2(NSTSGP(1)%I,(ITPFGS(IOBTP)-1)*NSMST+1,NSMST,NACIST)
!write(u6,*) ' IACIST and NACIST arrays'
!call IWRTMA(IACIST,1,NSMST,1,NSMST)
!call IWRTMA(NACIST,1,NSMST,1,NSMST)

! Generate offsets for I strings with given symmetry in each space

do IGAS=1,NGAS
  do ISMST=1,NSMST
    if (NNSTSGP(ISMST,IGAS) > 0) maxval(IGAS) = ISMST
  end do
  do ISMST=NSMST,1,-1
    if (NNSTSGP(ISMST,IGAS) > 0) minval(IGAS) = ISMST
  end do
end do
IFIRST = 1
NSTRINT = 0
2000 continue
if (IFIRST == 1) then
  do IGAS=1,NGASL-1
    ISMFGS(IGAS) = minval(IGAS)
  end do
else
  ! Next distribution of symmetries in NGAS -1
  call NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
  if (NONEW /= 0) goto 2001
end if
IFIRST = 0
! Symmetry of NGASL -1 spaces given, symmetry of full space
ISTSMM1 = 1
do IGAS=1,NGASL-1
  call SYMCOM(3,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
  ISTSMM1 = JSTSMM1
end do
! sym of SPACE NGASL
call SYMCOM(2,ISTSMM1,ISMGSN,ISPGPSM)
ISMFGS(NGASL) = ISMGSN
if (NTEST >= 200) then
  write(u6,*) ' next symmetry of NGASL spaces'
  call IWRTMA(ISMFGS,1,NGASL,1,NGASL)
end if
! Number of strings with this symmetry combination
NSTRII = 1
do IGAS=1,NGASL
  NSTRII = NSTRII*NNSTSGP(ISMFGS(IGAS),IGAS)
end do
! Offset for this symmetry distribution in IOFFI
IOFF = 1
MULT = 1
do IGAS=1,NGASL
  IOFF = IOFF+(ISMFGS(IGAS)-1)*MULT
  MULT = MULT*NSMST
end do

if (NTEST >= 1) then !SJS
  write(u6,*)
  write(u6,*) ' ============================'
  write(u6,*) ' If program is crashing here,'
  write(u6,*) ' LOFFI needs to be increased.'
  write(u6,*) ' ============================'
  write(u6,*)
end if
OFFI(IOFF) = real(NSTRINT,kind=wp)+1.001_wp
NSTRINT = NSTRINT+NSTRII
if (NTEST >= 200) write(u6,*) ' IOFF, OFFI(IOFF) NSTRII ',IOFF,OFFI(IOFF),NSTRII

if (NGASL-1 > 0) goto 2000
2001 continue

! Supergroup and symmetry of K strings

!M call NEWTYP(ISPGRPABS,1,IOBTP,KSPGRPABS)
!M call SYMCOM(2,IOBSM,KSM,ISPGPSM)
!M NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
!M if (NTEST >= 200) write(u6,*) ' KSM, KSPGPRABS, NKSTR : ',KSM,KSPGRPABS,NKSTR

! Gas structure of K strings

NGASL = 1
do IGAS=1,NGAS
  ITPFGS(IGAS) = ISPGPFTP(IGAS,KSPGRPABS)
  NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
  if (NELFGS(IGAS) > 0) NGASL = IGAS
end do
! Active group of K-strings
KACGRP = ITPFGS(IOBTP)
! Number of strings per symmetry distribution
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

! Loop over symmetry distribtions of K strings

KFIRST = 1
KSTRBS = 1
1000 continue
if (KFIRST == 1) then
  do IGAS=1,NGASL-1
    ISMFGS(IGAS) = minval(IGAS)
  end do
else
  ! Next distribution of symmetries in NGAS -1
  call NXTNUM3(ISMFGS,NGASL-1,MINVAL,MAXVAL,NONEW)
  if (NONEW /= 0) goto 1001
end if
KFIRST = 0
if (NTEST >= 200) then
  write(u6,*) ' next symmetry of NGASL-1 spaces'
  call IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
end if
! Symmetry of NGASL -1 spaces given, symmetry of total space
ISTSMM1 = 1
do IGAS=1,NGASL-1
  call SYMCOM(3,ISTSMM1,ISMFGS(IGAS),JSTSMM1)
  ISTSMM1 = JSTSMM1
end do
! required sym of SPACE NGASL
call SYMCOM(2,ISTSMM1,ISMGSN,KSM)
!write(u6,*) ' after  SYMCOM'
!write(u6,*) ' ngasl istsmm1 ksm',ngasl,istsmm1,ksm
ISMFGS(NGASL) = ISMGSN

do IGAS=NGASL+1,NGAS
  ISMFGS(IGAS) = 1
end do
if (NTEST >= 200) then
  write(u6,*) ' Next symmetry distribution'
  call IWRTMA(ISMFGS,1,NGAS,1,NGAS)
end if
! Number of strings of this symmetry distribution
NSTRIK = 1
do IGAS=1,NGASL
  NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
end do
! Offset for corresponding I strings
ISAVE = ISMFGS(IOBTP)
call SYMCOM(3,IOBSM,ISMFGS(IOBTP),IACSM)
ISMFGS(IOBTP) = IACSM
IOFF = 1
MULT = 1
do IGAS=1,NGAS
  IOFF = IOFF+(ISMFGS(IGAS)-1)*MULT
  MULT = MULT*NSMST
end do
ISMFGS(IOBTP) = ISAVE
IBSTRINI = int(OFFI(IOFF))
!write(u6,*) ' IOFF IBSTRINI ',IOFF,IBSTRINI
! Number of strings before active GAS space
NSTB = 1
do IGAS=1,IOBTP-1
  NSTB = NSTB*NNSTSGP(ISMFGS(IGAS),IGAS)
end do
! Number of strings before active GAS space
NSTA = 1
do IGAS=IOBTP+1,NGAS
  NSTA = NSTA*NNSTSGP(ISMFGS(IGAS),IGAS)
end do
! Number and offset for active group
!write(u6,*) ' IACSM = ',IACSM
NIAC = NACIST(IACSM)
IIAC = IACIST(IACSM)

NKAC = NNSTSGP(ISMFGS(IOBTP),IOBTP)
IKAC = IISTSGP(ISMFGS(IOBTP),IOBTP)
! I and K strings of given symmetry distribution
NKSD = NSTB*NKAC*NSTA
!write(u6,*) ' nstb nsta niac nkac ',nstb,nsta,niac,nkac
! Obtain annihilation n mapping for all strings of this type

NORBTS = NOBPTS(IOBTP,IOBSM)

!write(u6,*) ' KACGRP ',KACGRP
call ADSTN_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,STSTM(KACGRP,1)%I,STSTM(KACGRP,2)%I,IBORBSPS,IBORBSP,NORBTS,NKAC,NIAC,NKSTR, &
                 NELB,NACGSOB,I1,XI1S,SCLFAC)
KSTRBS = KSTRBS+NKSD
if (NGASL-1 > 0) goto 1000
1001 continue

9999 continue

if (NTEST >= 100) then
  write(u6,*) ' Output from ADSTN_GAS'
  write(u6,*) ' ====================='
  write(u6,*) ' Total number of K strings ',NKSTR
  if (NKSTR /= 0) then
    do IORB=IBORBSPS,IBORBSPS+NORBTS-1
      IORBR = IORB-IBORBSPS+1
      write(u6,*) ' Info for orbital ',IORB
      write(u6,*) ' Excited strings and sign'
      call IWRTMA(I1((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
      call WRTMAT(XI1S((IORBR-1)*NKSTR+1),1,NKSTR,1,NKSTR)
    end do
  end if
end if

! PAM Mars-2006: This flush moved outside of this subroutine
!call MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'ADSTN ')

end subroutine ADSTN_GAS
