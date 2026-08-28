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
! Copyright (C) 1996, Jeppe Olsen                                      *
!***********************************************************************

subroutine GSBBD2B_MCLR(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,SB,CB, &
                        NOBPTS,IOBPTS,MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,NSM,IUSEAB,CJRES,SIRES,NORB,ieaw)
! alpha-beta contribution to two-particle density matrix
! from given c-block and s-block.
!
! =====
! Input
! =====
! IASM,IATP : Symmetry and type of alpha  strings in sigma
! IBSM,IBTP : Symmetry and type of beta   strings in sigma
! JASM,JATP : Symmetry and type of alpha  strings in C
! JBSM,JBTP : Symmetry and type of beta   strings in C
! NIA,NIB   : Number of alpha-(beta-) strings in sigma
! NJA,NJB   : Number of alpha-(beta-) strings in C
! IAGRP     : String group of alpha strings
! IBGRP     : String group of beta strings
! IAEL1(3)  : Number of electrons in RAS1(3) for alpha strings in sigma
! IBEL1(3)  : Number of electrons in RAS1(3) for beta  strings in sigma
! JAEL1(3)  : Number of electrons in RAS1(3) for alpha strings in C
! JBEL1(3)  : Number of electrons in RAS1(3) for beta  strings in C
! CB        : Input C block
! NTSOB     : Number of orbitals per type and symmetry
! IBTSOB    : base for orbitals of given type and symmetry
! IBORB     : Orbitals of given type and symmetry
! NSM       : Number of symmetries of orbitals
! MAXK      : Largest number of inner resolution strings treated at simult.
!
! ======
! Output
! ======
! SB : updated sigma block
!
! =======
! Scratch
! =======
! I1, XI1S : at least MXSTSO : Largest number of strings of given
!            type and symmetry
! I2, XI2S : at least MXSTSO : Largest number of strings of given
!            type and symmetry
! X        : Space for block of two-electron integrals
!
! Jeppe Olsen, Fall of 1996

use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: RHO2(*)
integer(kind=iwp), intent(in) :: IASM, IATP, IBSM, IBTP, NIA, NIB, JASM, JATP, JBSM, JBTP, NJA, NJB, IAGRP, IBGRP, IAOC(*), &
                                 IBOC(*), JAOC(*), JBOC(*), NOBPTS(3,*), IOBPTS(3,*), MAXK, NSM, IUSEAB, NORB, ieaw
integer(kind=iwp), intent(out) :: NGAS
real(kind=wp), intent(in) :: SB(*), CB(*)
real(kind=wp), intent(_OUT_) :: XI1S(*), XI2S(*), XI3S(*), XI4S(*), X(*), CJRES(*), SIRES(*)
integer(kind=iwp), intent(_OUT_) :: I1(*), I2(*), I3(*), I4(*)
integer(kind=iwp) :: IDOCOMP, II, IJSM, IJTYP, IKABTC, IKORD, IOFF, ISM, ITP(3), ITYP, itype, JJ, JOFF, JSM, JTP(3), JTYP, KABOT, &
                     KAEND, KATOP, KBBOT, KBEND, KBTOP, KLSM, KLTYP, KOFF, KSM, KTP(3), KTYP, LKABTC, LOFF, LSM, LTP(3), LTYP, NI, &
                     NIJTYP, NJ, NK, NKABTC, NKAEFF, NKASTR, NKBSTR, NKLTYP, NL

NGAS = 3

! Symmetry of allowed excitations
IJSM = Mul(IASM,JASM)
KLSM = Mul(IBSM,JBSM)
itype = 2
if (ieaw == 1) itype = 3
if ((IJSM == 0) .or. (KLSM == 0)) return
! Types of SX that connects the two strings
call SXTYP_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC)
call SXTYP_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC)
if ((NIJTYP == 0) .or. (NKLTYP == 0)) return
do IJTYP=1,NIJTYP
  ITYP = ITP(IJTYP)
  JTYP = JTP(IJTYP)
  do ISM=1,NSM
    JSM = Mul(ISM,IJSM)
    if (JSM == 0) cycle
    IOFF = IOBPTS(ITYP,ISM)
    JOFF = IOBPTS(JTYP,JSM)
    NI = NOBPTS(ITYP,ISM)
    NJ = NOBPTS(JTYP,JSM)
    if ((NI == 0) .or. (NJ == 0)) cycle
    !EAW
    ! Find Ka strings that connect with Ja strings for given group of Jorbs
    KABOT = 1
    ! Obtain all strings
    KATOP = -1
    call ADST(JOFF,NJ,JATP,JASM,IAGRP,KABOT,KATOP,I1,XI1S,MAXK,NKASTR,KAEND)
    call ADST(IOFF,NI,IATP,IASM,IAGRP,KABOT,KATOP,I3,XI3S,MAXK,NKASTR,KAEND)
    IDOCOMP = 1
    if (IDOCOMP == 1) then
      call COMPRS2LST(I1,XI1S,NJ,I3,XI3S,NI,NKASTR,NKAEFF)
    else
      NKAEFF = NKASTR
    end if

    ! Loop over batches of KA strings
    NKABTC = NKAEFF/MAXK
    if (NKABTC*MAXK < NKAEFF) NKABTC = NKABTC+1
    do IKABTC=1,NKABTC
      KABOT = (IKABTC-1)*MAXK+1
      KATOP = min(KABOT+MAXK-1,NKAEFF)
      LKABTC = KATOP-KABOT+1
      ! Obtain C(ka,J,JB) for Ka in batch
      do JJ=1,NJ
        call GET_CKAJJB(CB,NJ,NJA,CJRES,LKABTC,NJB,JJ,I1(KABOT+(JJ-1)*NKASTR),XI1S(KABOT+(JJ-1)*NKASTR))
      end do
      ! Obtain S(ka,i,Ib) for Ka in batch
      do II=1,NI
        call GET_CKAJJB(SB,NI,NIA,SIRES,LKABTC,NIB,II,I3(KABOT+(II-1)*NKASTR),XI3S(KABOT+(II-1)*NKASTR))
      end do

      do KLTYP=1,NKLTYP
        KTYP = KTP(KLTYP)
        LTYP = LTP(KLTYP)

        do KSM=1,NSM
          LSM = Mul(KSM,KLSM)
          if (LSM == 0) cycle
          KOFF = IOBPTS(KTYP,KSM)
          LOFF = IOBPTS(LTYP,LSM)
          NK = NOBPTS(KTYP,KSM)
          NL = NOBPTS(LTYP,LSM)
          ! If IUSEAB is used, only terms with i >= k will be generated so
          IKORD = 0
          if ((IUSEAB == 1) .and. (ISM > KSM)) cycle
          if ((IUSEAB == 1) .and. (ISM == KSM) .and. (ITYP < KTYP)) cycle
          if ((IUSEAB == 1) .and. (ISM == KSM) .and. (ITYP == KTYP)) IKORD = 1

          if ((NK == 0) .or. (NL == 0)) cycle
          !EAW
          ! Obtain all connections a+l!Kb> = +/-/0!Jb>
          ! NKBSTR must be given as input
          ! obtain cb(KA,KB,jl) =  sum(JA,JB)<KA!a la!JA><KB!a jb !JB>C(JA,JB)

          KBBOT = 1
          ! Obtain all strings
          KBTOP = -1
          call ADST(LOFF,NL,JBTP,JBSM,IBGRP,KBBOT,KBTOP,I2,XI2S,MAXK,NKBSTR,KBEND)
          call ADST(KOFF,NK,IBTP,IBSM,IBGRP,KBBOT,KBTOP,I4,XI4S,MAXK,NKBSTR,KBEND)

          !if (NKBSTR == 0) cycle
          X(1:NI*NJ*NK*NL) = Zero

          call ABTOR2(SIRES,CJRES,LKABTC,NKBSTR,X,NI,NJ,NK,NL,NKBSTR,I4,XI4S,I2,XI2S,IKORD)
          ! contributions to Rho2(ij,kl) has been obtained, scatter out
          call ADTOR2_MCLR(RHO2,X,itype,NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)

        end do
      end do
    end do
    ! End of loop over partitioning of alpha strings
  end do
end do

end subroutine GSBBD2B_MCLR
