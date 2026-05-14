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
! Copyright (C) 1991,1994-1997, Jeppe Olsen                            *
!               2012, Giovanni Li Manni                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ADAST_GAS(IOBSM,IOBTP,NIGRP,IGRP,ISPGPSM,I1,XI1S,NKSTR,KACT,SCLFAC,IAC)
! Obtain creation or annihilation mapping
!
! IAC = 2 : Creation map
! a+IORB !KSTR> = +/-!ISTR>
!
! IAC = 1 : Annihilation map
! a IORB !KSTR> = +/-!ISTR>
!
! for orbitals of symmetry IOBSM and type IOBTP
! and Istrings defined by the NIGRP groups IGRP and symmetry ISPGPSM
!
! The results are given in the form
! I1(KSTR,IORB) =  ISTR if A+IORB !KSTR> = +/-!ISTR>
! (numbering relative to TS start)
! Above +/- is stored in XI1S
!
! if some nonvanishing excitations were found, KACT is set to 1,
! else it is zero
!
!
! Jeppe Olsen, Winter of 1991
!              January 1994 : modified to allow for several orbitals
!              August 95    : GAS version
!              October 96   : Improved version
!              September 97 : annihilation mappings added
!                             I groups defined by IGRP
!
! Giovanni Li Manni, February 2012
! Smart Loop over symmetry distributions
! in order to make the code faster

use Symmetry_Info, only: Mul
use lucia_data, only: IBGPSTR, IGSFGP, IOBPTS, ISMDFGP, ISMSCR, ISTAC, ISTSGP, LOFFI, MXPNGAS, MXPNSMST, NACTSYM, NELFGP, NGPSTR, &
                      NGRP, NIRREP, NOBPT, NOBPTS, NSTSGP, STSTM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IOBSM, IOBTP, NIGRP, IGRP(NIGRP), ISPGPSM, IAC
integer(kind=iwp), intent(inout) :: I1(*)
real(kind=wp), intent(inout) :: XI1S(*)
integer(kind=iwp), intent(out) :: NKSTR, KACT
real(kind=wp), intent(in) :: SCLFAC
integer(kind=iwp) :: I, IACGAS, IACGRP, IACIST(MXPNSMST), IACSM, IBORBSP, IBORBSPS, IBSTRINI, IDELTA, IEC, IGAS, IIAC, &
                     IISTSGP(MXPNSMST,MXPNGAS), IKAC, ISAVE, ISMFGS(MXPNGAS), JGRP, KACGRP, KFIRST, KGRP(MXPNGAS), KSM, KSTRBS, &
                     LROW_IN, MNVLI(MXPNGAS), MXVLI(MXPNGAS), NACIST(MXPNSMST), NELB, NGASL, NIAC, NIEL, NIGASL, NKAC, NKDIST, &
                     NKEL, NKSD, NNSTSGP(MXPNSMST,MXPNGAS), NONEW, NORBT, NORBTS, NSTA, NSTB, NSTRIK
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IORB, IORBR
#endif
integer(kind=iwp), allocatable :: IOFFI(:)
integer(kind=iwp), external :: IOFF_SYM_DIST

! Will be stored as a matrix of dimension
! (NKSTR,*), Where NKSTR is the number of K-strings of
! correct symmetry . Nk is provided by this routine.

#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' ===================='
write(u6,*) ' ADAST_GAS in service'
write(u6,*) ' ===================='
write(u6,*)
write(u6,*) 'GAS space (IOBTP), Symm of it (IOBSM) :',IOBTP,IOBSM
write(u6,*) ' Supergroup in action :'
write(u6,'(A,I3  )') ' Number of active spaces ',NIGRP
write(u6,'(A,20I3)') ' The active groups       ',(IGRP(I),I=1,NIGRP)
write(u6,*) '  Symmetry of supergroup : ',ISPGPSM
write(u6,*) ' SCLFAC = ',SCLFAC
if (IAC == 1) then
  write(u6,*) ' Annihilation mapping'
else if (IAC == 2) then
  write(u6,*) ' Creation mapping'
else
  write(u6,*) ' Unknown IAC parameter in ADAST ',IAC
  call SYSABENDMSG('lucia_util/adast_gas','Internal error','')
end if
#endif
! A few preparations
NORBTS = NOBPTS(IOBTP,IOBSM)
NORBT = NOBPT(IOBTP)
IACGAS = IOBTP
! First orbital of given GASpace
IBORBSP = sum(NOBPT(1:IOBTP-1))+1
! First orbital of given GASpace and Symmetry
IBORBSPS = IOBPTS(IOBTP,IOBSM)
#ifdef _DEBUGPRINT_
write(u6,*) ' NORBTS per GAS and sym          :',NORBTS
write(u6,*) ' NORBT per GAS                   :',NORBT
write(u6,*) ' IACGAS GAS involved             :',IACGAS
write(u6,*) ' IBORBSP 1st orb per GAS         :',IBORBSP
write(u6,*) ' IBORBSPS 1st orb per GAS and sym:',IBORBSPS
#endif

!===================================================
! K strings : Supergroup, symmetry and distributions
!===================================================
if (IAC == 1) then
  IDELTA = +1
else
  IDELTA = -1
end if
! Is required mapping contained within current set of maps?
! a:) Is active GASpace included in IGRP - must be
IACGRP = 0
do JGRP=1,NIGRP
  if (IGSFGP(IGRP(JGRP)) == IACGAS) IACGRP = JGRP
end do
! Note : IACGRP is not the actual active group, it is the address of the
!        active group in IGRP
if (IACGRP == 0) then
  write(u6,*) ' ADAST in problems'
  write(u6,*) ' Active GASpace not included in IGRP'
  write(u6,*) ' Active GASpace : ',IACGAS
  write(u6,'(A,20I3)') ' The active groups       ',(IGRP(I),I=1,NIGRP)
  call SYSABENDMSG('lucia_util/adast_gas','Internal error','')
end if
! b:) active group in K strings
NIEL = NELFGP(IGRP(IACGRP))
NKEL = NIEL+IDELTA
#ifdef _DEBUGPRINT_
write(u6,*) ' NIEL and NKEL ',NIEL,NKEL
#endif
if ((NKEL == -1) .or. (NKEL == NOBPT(IACGAS)+1)) then
  ! No strings with this number of elecs - be happy : No work
  NKSTR = 0
  KACT = 0
  KACGRP = 0
else
  ! Find group with NKEL electrons in IACGAS
  KACGRP = 0
  do JGRP=IBGPSTR(IACGAS),IBGPSTR(IACGAS)+NGPSTR(IACGAS)-1
    if (NELFGP(JGRP) == NKEL) KACGRP = JGRP
  end do
# ifdef _DEBUGPRINT_
  write(u6,*) ' KACGRP = ',KACGRP
# endif
  ! KACGRP is the Active group itself
  if (KACGRP == 0) then
    write(u6,*) ' ADAST : cul de sac, active K group not found'
    write(u6,*) ' GAS space and number of electrons ',IACGAS,NKEL
    call SYSABENDMSG('lucia_util/adast_gas','Internal error','')
  end if
  ! Okay active K group was found and is nontrivial
  KSM = Mul(IOBSM,ISPGPSM)
  ! The K supergroup
  KGRP(1:NIGRP) = IGRP(:)
  KGRP(IACGRP) = KACGRP
  ! Number of strings and symmetry distributions of K strings
  call NST_SPGRP(NIGRP,KGRP,KSM,NSTSGP,NIRREP,NKSTR,NKDIST)
# ifdef _DEBUGPRINT_
  write(u6,*) 'KSM,NKSTR,NKDIST:',KSM,NKSTR,NKDIST
# endif
  if (NKSTR /= 0) then
    ! Last active space in K strings and number of strings per group and sym
    NGASL = 1
    do JGRP=1,NIGRP
      if (NELFGP(KGRP(JGRP)) > 0) NGASL = JGRP
      NNSTSGP(1:NIRREP,JGRP) = NSTSGP((KGRP(JGRP)-1)*NIRREP+1:KGRP(JGRP)*NIRREP)
      IISTSGP(1:NIRREP,JGRP) = ISTSGP((KGRP(JGRP)-1)*NIRREP+1:KGRP(JGRP)*NIRREP)
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) 'NGASL',NGASL
#   endif
    ! MIN/MAX for Kstrings
    !call MINMAX_FOR_SYM_DIST(NIGRP,KGRP,MNVLK,MXVLK,NKDIST_TOT)
#   ifdef _DEBUGPRINT_
    !write(u6,*) 'MNVLK and MXVLK'
    !call IWRTMA(MNVLK,1,NIGRP,1,NIGRP)
    !call IWRTMA(MXVLK,1,NIGRP,1,NIGRP)
#   endif
    ! (NKDIST_TOT is number of distributions, all symmetries)
    ! =========
    ! I Strings
    ! =========
    call mma_allocate(IOFFI,LOFFI,label='IOFFI')
    ! Generate symmetry distributions of I strings with given symmetry
    call TS_SYM_PNT2(IGRP,NIGRP,MXVLI,MNVLI,ISPGPSM,IOFFI,LOFFI)
    ! Offset and dimension for active group in I strings
    IACIST(1:NIRREP) = ISTSGP((IGRP(IACGRP)-1)*NIRREP+1:IGRP(IACGRP)*NIRREP)
    NACIST(1:NIRREP) = NSTSGP((IGRP(IACGRP)-1)*NIRREP+1:IGRP(IACGRP)*NIRREP)
    ! Last entry in IGRP with a nonvanisking number of strings
    NIGASL = 1
    do JGRP=1,NIGRP
      if (NELFGP(IGRP(JGRP)) > 0) NIGASL = JGRP
    end do
    ! Number of electrons before active space
    NELB = 0
    do JGRP=1,IACGRP-1
      NELB = NELB+NELFGP(IGRP(JGRP))
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) ' NELB = ',NELB
#   endif
    I1(1:NORBTS*NKSTR) = 0
    ! Loop over symmetry distribtions of K strings
    KFIRST = 1
    KSTRBS = 1
    ISMFGS(1:NIGRP) = 1
    do
      !GLM if (KFIRST == 1) ISMFGS(1:NIGRP) = MNVLI(IGAS)
      ! Next distribution
      call NEXT_SYM_DISTR_NEW(NIRREP,NGRP,KGRP,NIGRP,ISMFGS,KSM,KFIRST,NONEW,ISMDFGP,NACTSYM,ISMSCR)
      !GLM   if (NONEW == 1) exit
      !GLM end if
#     ifdef _DEBUGPRINT_
      write(u6,*) ' Symmetry distribution'
      call iwrtma(ISMFGS,1,NIGRP,1,NIGRP)
#     endif
      if (NONEW == 1) exit
      KFIRST = 0
      ! Number of strings of this symmetry distribution
      NSTRIK = 1
      do IGAS=1,NGASL
        NSTRIK = NSTRIK*NNSTSGP(ISMFGS(IGAS),IGAS)
      end do
      ! Offset for corresponding I strings
      ISAVE = ISMFGS(IACGRP)
      IACSM = Mul(IOBSM,ISMFGS(IACGRP))
      ISMFGS(IACGRP) = IACSM
      IBSTRINI = IOFF_SYM_DIST(ISMFGS,NIGASL,IOFFI,MXVLI,MNVLI)
      !write(u6,*) 'IBSTRINI :',IBSTRINI
      ISMFGS(IACGRP) = ISAVE
      ! Number of strings before active GAS space
      NSTB = 1
      !do IGAS=1,IOBTP-1
      do IGAS=1,IACGRP-1
        NSTB = NSTB*NNSTSGP(ISMFGS(IGAS),IGAS)
      end do
      ! Number of strings After active GAS space
      NSTA = 1
      !do IGAS=IOBTP+1,NIGRP
      do IGAS=IACGRP+1,NIGRP
        NSTA = NSTA*NNSTSGP(ISMFGS(IGAS),IGAS)
      end do
      ! Number and offset for active group
      NIAC = NACIST(IACSM)
      IIAC = IACIST(IACSM)
      NKAC = NNSTSGP(ISMFGS(IACGRP),IACGRP)
      IKAC = IISTSGP(ISMFGS(IACGRP),IACGRP)
      ! I and K strings of given symmetry distribution
      NKSD = NSTB*NKAC*NSTA
#     ifdef _DEBUGPRINT_
      write(u6,*) ' nstb nsta niac nkac ',nstb,nsta,niac,nkac
#     endif
      ! Obtain annihilation/creation mapping for all strings of this type
      ! Are group mappings in expanded or compact form
      if ((IAC == 1) .and. (ISTAC(KACGRP,2) == 0)) then
        IEC = 2
        LROW_IN = NKEL
      else
        IEC = 1
        LROW_IN = NORBT
      end if

      if (NSTA*NSTB*NIAC*NKAC /= 0) &
        call ADAST_GASSM(NSTB,NSTA,IKAC,IIAC,IBSTRINI,KSTRBS,STSTM(KACGRP,1)%A,STSTM(KACGRP,2)%A,IBORBSPS,IBORBSP,NORBTS,NKAC, &
                         NIAC,NKSTR,NELB,I1,XI1S,SCLFAC,IAC,LROW_IN,IEC)
      KSTRBS = KSTRBS+NKSD
    end do

    call mma_deallocate(IOFFI)
  end if
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from ADAST_GAS'
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

end subroutine ADAST_GAS
