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

!#define _DEBUGPRINT_
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

use Symmetry_Info, only: Mul
use lucia_data, only: IBSPGPFTP, IOBPTS, ISPGPFTP, ISTSGP, MXPNGAS, MXPNSMST, NELFGP, NGAS, NIRREP, NOBPT, NOBPTS, NSTFSMSPGP, &
                      NSTSGP, STSTM
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: OFFI(*), XI1S(*)
integer(kind=iwp), intent(in) :: IOBSM, IOBTP, ISPGP, ISPGPSM, ISPGPTP
integer(kind=iwp), intent(_OUT_) :: I1(*)
integer(kind=iwp), intent(out) :: NKSTR
real(kind=wp), intent(in) :: SCLFAC
integer(kind=iwp) :: IACIST(MXPNSMST), IACSM, IBORBSP, IBORBSPS, IBSTRINI, IFIRST, IGAS, IIAC, IISTSGP(MXPNSMST,MXPNGAS), IKAC, &
                     IOFF, ISAVE, ISMFGS(MXPNGAS), ISMST, ISPGRPABS, ISTSMM1, ITPFGS(MXPNGAS), KACGRP, KFIRST, KSM, KSPGRPABS, &
                     KSTRBS, MNVAL(MXPNGAS), MULT, MXVAL(MXPNGAS), NACGSOB, NACIST(MXPNSMST), NELB, NELFGS(MXPNGAS), NGASL, NIAC, &
                     NKAC, NKSD, NNSTSGP(MXPNSMST,MXPNGAS), NONEW, NORBTS, NSTA, NSTB, NSTRII, NSTRIK, NSTRINT
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IORB, IORBR
#endif
integer(kind=iwp), parameter :: MXLNGAS = 20

! Will be stored as an matrix of dimension
! (NKSTR,*), Where NKSTR is the number of K-strings of
! correct symmetry . Nk is provided by this routine.

if (NGAS > MXLNGAS) then
  write(u6,*) ' Ad hoc programming in ADSTN (IOFFI)'
  write(u6,*) ' Must be changed - or redimensioned'
  !stop 'ADST : IOFFI problem'
  call SYSABENDMSG('lucia_util/adstn_gas','Internal error','')
end if

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ===================='
write(u6,*) ' ADSTN_GAS in service'
write(u6,*) ' ===================='
write(u6,*)
write(u6,*) '  IOBTP IOBSM : ',IOBTP,IOBSM
write(u6,*) '  ISPGP ISPGPSM ISPGPTP :  ',ISPGP,ISPGPSM,ISPGPTP
#endif

!if (SCLFAC /= One) then
!  write(u6,*) ' Problemo : ADSTN_GAS'
!  write(u6,*) ' SCLFAC /= 1'
!end if

! Supergroup and symmetry of K strings

ISPGRPABS = IBSPGPFTP(ISPGPTP)-1+ISPGP
call NEWTYP(ISPGRPABS,1,IOBTP,KSPGRPABS)
KSM = Mul(IOBSM,ISPGPSM)
NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
#ifdef _DEBUGPRINT_
write(u6,*) ' KSM, KSPGPRABS, NKSTR : ',KSM,KSPGRPABS,NKSTR
#endif
if (NKSTR /= 0) then

  NORBTS = NOBPTS(IOBTP,IOBSM)
  XI1S(1:NORBTS*NKSTR) = Zero
  I1(1:NORBTS*NKSTR) = 0

  ! First orbital of given GASSpace
  IBORBSP = sum(NOBPT(1:IOBTP-1))+1
  ! First orbital of fiven GASSPace and Symmetry
  IBORBSPS = IOBPTS(IOBTP,IOBSM)

  ! Information about I strings
  ! ===========================

  ! structure of group of strings defining I strings
  NGASL = 1
  ITPFGS(:) = 0
  do IGAS=1,NGAS
    ITPFGS(IGAS) = ISPGPFTP(IGAS,ISPGRPABS)
    NELFGS(IGAS) = NELFGP(ITPFGS(IGAS))
    if (NELFGS(IGAS) > 0) NGASL = IGAS
  end do
  ! Number of electrons before active type
  NELB = sum(NELFGS(1:IOBTP-1))
  ! Number of electrons in active space
  NACGSOB = NOBPT(IOBTP)

  ! Number of strings per symmetry for each symmetry
  do IGAS=1,NGAS
    NNSTSGP(1:NIRREP,IGAS) = NSTSGP((ITPFGS(IGAS)-1)*NIRREP+1:ITPFGS(IGAS)*NIRREP)
  end do
  ! Offset and dimension for active group in I strings
  IACIST(1:NIRREP) = ISTSGP((ITPFGS(IOBTP)-1)*NIRREP+1:ITPFGS(IOBTP)*NIRREP)
  NACIST(1:NIRREP) = NSTSGP((ITPFGS(IOBTP)-1)*NIRREP+1:ITPFGS(IOBTP)*NIRREP)
  !write(u6,*) ' IACIST and NACIST arrays'
  !call IWRTMA(IACIST,1,NIRREP,1,NIRREP)
  !call IWRTMA(NACIST,1,NIRREP,1,NIRREP)

  ! Generate offsets for I strings with given symmetry in each space

  do IGAS=1,NGAS
    do ISMST=1,NIRREP
      if (NNSTSGP(ISMST,IGAS) > 0) MXVAL(IGAS) = ISMST
    end do
    do ISMST=NIRREP,1,-1
      if (NNSTSGP(ISMST,IGAS) > 0) MNVAL(IGAS) = ISMST
    end do
  end do
  IFIRST = 1
  NSTRINT = 0
  do
    if (IFIRST == 1) then
      ISMFGS(1:NGASL-1) = MNVAL(1:NGASL-1)
    else
      ! Next distribution of symmetries in NGAS -1
      call NXTNUM3(ISMFGS,NGASL-1,MNVAL,MXVAL,NONEW)
      if (NONEW /= 0) exit
    end if
    IFIRST = 0
    ! Symmetry of NGASL -1 spaces given, symmetry of full space
    ISTSMM1 = 1
    do IGAS=1,NGASL-1
      ISTSMM1 = Mul(ISTSMM1,ISMFGS(IGAS))
    end do
    ! sym of SPACE NGASL
    ISMFGS(NGASL) = Mul(ISTSMM1,ISPGPSM)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' next symmetry of NGASL spaces'
    call IWRTMA(ISMFGS,1,NGASL,1,NGASL)
#   endif
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
      MULT = MULT*NIRREP
    end do

#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) ' ============================'
    write(u6,*) ' If program is crashing here,'
    write(u6,*) ' LOFFI needs to be increased.'
    write(u6,*) ' ============================'
    write(u6,*)
#   endif
    OFFI(IOFF) = real(NSTRINT,kind=wp)+1.001_wp
    NSTRINT = NSTRINT+NSTRII
#   ifdef _DEBUGPRINT_
    write(u6,*) ' IOFF, OFFI(IOFF) NSTRII ',IOFF,OFFI(IOFF),NSTRII
#   endif

    if (NGASL <= 1) exit
  end do

  ! Supergroup and symmetry of K strings

  !M call NEWTYP(ISPGRPABS,1,IOBTP,KSPGRPABS)
  !M KSM = Mul(IOBSM,ISPGPSM)
  !M NKSTR = NSTFSMSPGP(KSM,KSPGRPABS)
# ifdef _DEBUGPRINT_
  !M write(u6,*) ' KSM, KSPGPRABS, NKSTR : ',KSM,KSPGRPABS,NKSTR
# endif

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

  ! Loop over symmetry distribtions of K strings

  KFIRST = 1
  KSTRBS = 1
  do
    if (KFIRST == 1) then
      ISMFGS(1:NGASL-1) = MNVAL(1:NGASL-1)
    else
      ! Next distribution of symmetries in NGAS -1
      call NXTNUM3(ISMFGS,NGASL-1,MNVAL,MXVAL,NONEW)
      if (NONEW /= 0) exit
    end if
    KFIRST = 0
#   ifdef _DEBUGPRINT_
    write(u6,*) ' next symmetry of NGASL-1 spaces'
    call IWRTMA(ISMFGS,NGASL-1,1,NGASL-1,1)
#   endif
    ! Symmetry of NGASL -1 spaces given, symmetry of total space
    ISTSMM1 = 1
    do IGAS=1,NGASL-1
      ISTSMM1 = Mul(ISTSMM1,ISMFGS(IGAS))
    end do
    ! required sym of SPACE NGASL
    !write(u6,*) ' after  SYMCOM'
    !write(u6,*) ' ngasl istsmm1 ksm',ngasl,istsmm1,ksm
    ISMFGS(NGASL) = Mul(ISTSMM1,KSM)

    ISMFGS(NGASL+1:NGAS) = 1
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Next symmetry distribution'
    call IWRTMA(ISMFGS,1,NGAS,1,NGAS)
#   endif
    ! Number of strings of this symmetry distribution
    NSTRIK = 1
    do IGAS=1,NGASL
      NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
    end do
    ! Offset for corresponding I strings
    ISAVE = ISMFGS(IOBTP)
    IACSM = Mul(IOBSM,ISMFGS(IOBTP))
    ISMFGS(IOBTP) = IACSM
    IOFF = 1
    MULT = 1
    do IGAS=1,NGAS
      IOFF = IOFF+(ISMFGS(IGAS)-1)*MULT
      MULT = MULT*NIRREP
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
    call ADSTN_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,STSTM(KACGRP,1)%A,STSTM(KACGRP,2)%A,IBORBSPS,IBORBSP,NORBTS,NKAC,NIAC, &
                     NKSTR,NELB,NACGSOB,I1,XI1S,SCLFAC)
    KSTRBS = KSTRBS+NKSD
    if (NGASL <= 1) exit
  end do

end if

#ifdef _DEBUGPRINT_
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
#endif

! PAM Mars-2006: This flush moved outside of this subroutine
!call MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'ADSTN ')

end subroutine ADSTN_GAS
