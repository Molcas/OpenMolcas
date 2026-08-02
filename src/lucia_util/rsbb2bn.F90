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
! Copyright (C) 1991-1994,1996,1997,2000, Jeppe Olsen                  *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine RSBB2BN(IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,SB,CB,NOBPTS,MAXK, &
                   I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,NSMOB,IUSEAB,CJRES,SIRES,SCLFAC,IPHGAS)
! SUBROUTINE RSBB2BN --> 52
!
! Combined alpha-beta double excitation
! contribution from given C block to given S block
! If IUSAB only half the terms are constructed
!
! =====
! Input
! =====
!
! IASM,IATP : Symmetry and type of alpha  strings in sigma
! IBSM,IBTP : Symmetry and type of beta   strings in sigma
! JASM,JATP : Symmetry and type of alpha  strings in C
! JBSM,JBTP : Symmetry and type of beta   strings in C
! NIA,NIB : Number of alpha-(beta-) strings in sigma
! NJA,NJB : Number of alpha-(beta-) strings in C
! IAGRP : String group of alpha strings
! IBGRP : String group of beta strings
! IAEL1(3) : Number of electrons in RAS1(3) for alpha strings in sigma
! IBEL1(3) : Number of electrons in RAS1(3) for beta  strings in sigma
! JAEL1(3) : Number of electrons in RAS1(3) for alpha strings in C
! JBEL1(3) : Number of electrons in RAS1(3) for beta  strings in C
! CB   : Input C block
! NTSOB  : Number of orbitals per type and symmetry
! IBTSOB : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! NSMOB  : Number of symmetries of orbitals
! MAXK   : Largest number of inner resolution strings treated at simult.
!
! ======
! Output
! ======
! SB : updated sigma block
!
! =======
! Scratch
! =======
!
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! I2, XI2S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! XINT  : Space for two electron integrals
!
! Jeppe Olsen, Winter of 1991
!
! Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
! October 1993 : IUSEAB added
! January 1994 : Loop restructured + CJKAIB introduced
! February 1994 : Fetching and adding to transposed blocks
! October 96 : New routines for accessing annihilation information
!             Cleaned and shaved, only IROUTE = 3 option active
! October 97 : allowing for N-1/N+1 switch
!
! Last change : Aug 2000

use Symmetry_Info, only: Mul
use Para_Info, only: MyRank, nProcs
use lucia_data, only: MXPNGAS, TSIGMA
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IASM, IATP, IBSM, IBTP, NIA, NIB, JASM, JATP, JBSM, JBTP, NJA, NJB, IAGRP, IBGRP, NGAS, IAOC(*), &
                                 IBOC(*), JAOC(*), JBOC(*), NOBPTS(MXPNGAS,*), MAXK, NSMOB, IUSEAB, IPHGAS(*)
real(kind=wp), intent(inout) :: SB(*), XI1S(*), XI2S(*), XI3S(*), XI4S(*)
real(kind=wp), intent(in) :: CB(*), SCLFAC
integer(kind=iwp), intent(inout) :: I1(*), I2(*), I3(*), I4(*)
real(kind=wp), intent(_OUT_) :: XINT(*), CJRES(*), SIRES(*)
integer(kind=iwp) :: IASPGP(20), IBSPGP(20), ICOUL, IDOCOMP, II, IJ_DIM(2), IJ_REO(2), IJ_SYM(2), IJ_TYP(2), IJAC, IJSM, IJTYP, &
                     IKABTC, IKORD, IROUTE, ISM, ITP(20), ITYP, IXCHNG, JASPGP(20), JBSPGP(20), JJ, JSM, JTP(20), JTYP, KABOT, &
                     KACT, KATOP, KL_DIM(2), KL_REO(2), KL_SYM(2), KL_TYP(2), KLAC, KLSM, KLTYP, KSM, KTP(20), KTYP, LKABTC, LSM, &
                     LTP(20), LTYP, NI, NIJTYP, NJ, NK, NKABTC, NKABTCSZ, NKAEFF, NKASTR, NKBSTR, NKLTYP, NL
real(kind=wp) :: CPU, CPU0, CPU1, FACS, SIGNIJ2, SIGNKL, WALL, WALL0, WALL1

#ifdef _DEBUGPRINT_
write(u6,*) ' ================'
write(u6,*) ' RSBB2BN speaking'
write(u6,*) ' ================'
#endif
!  Groups defining each supergroup
call GET_SPGP_INF(IATP,IAGRP,IASPGP)
call GET_SPGP_INF(JATP,IAGRP,JASPGP)
call GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
call GET_SPGP_INF(JBTP,IBGRP,JBSPGP)

! Symmetry of allowed excitations
IJSM = Mul(IASM,JASM)
KLSM = Mul(IBSM,JBSM)
if ((IJSM == 0) .or. (KLSM == 0)) return
#ifdef _DEBUGPRINT_
write(u6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
write(u6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
#endif
! Types of SX that connects the two strings
call SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC,IPHGAS)
call SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC,IPHGAS)
if ((NIJTYP == 0) .or. (NKLTYP == 0)) return

do IJTYP=1,NIJTYP

  ITYP = ITP(IJTYP)
  JTYP = JTP(IJTYP)
  do ISM=1,NSMOB
    JSM = Mul(ISM,IJSM)
    if (JSM == 0) cycle
    NI = NOBPTS(ITYP,ISM)
    NJ = NOBPTS(JTYP,JSM)
    if ((NI == 0) .or. (NJ == 0)) cycle
    ! Should N-1 or N+1 projection be used for alpha strings
    IJ_TYP(1) = ITYP
    IJ_TYP(2) = JTYP
    !IJ_AC(1) = 2
    !IJ_AC(2) = 1
    !NOP = 2
    !if (IUSE_PH == 1) then
    !  call ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TYP,IJ_AC,IJ_REO,SIGNIJ)
    !else
    ! Enforced a+ a
    IJ_REO(1) = 1
    IJ_REO(2) = 2
    !end if
    ! Two choices here :
    !  1 : <Ia!a+ ia!Ka><Ja!a+ ja!Ka> (good old creation mapping)
    !  2 :-<Ia!a  ja!Ka><Ja!a  ia!Ka> + delta(i,j)
    !write(u6,*) ' RSBB2BN : IOP_REO : ',(IOP_REO(II),II=1,2)
    if ((IJ_REO(1) == 1) .and. (IJ_REO(2) == 2)) then
      ! Business as usual i.e. creation map
      IJAC = 2
      SIGNIJ2 = SCLFAC

      IJ_DIM(1) = NI
      IJ_DIM(2) = NJ
      IJ_SYM(1) = ISM
      IJ_SYM(2) = JSM
      IJ_TYP(1) = ITYP
      IJ_TYP(2) = JTYP
    else
      ! Terra Nova, annihilation map
      IJAC = 1
      SIGNIJ2 = -SCLFAC

      IJ_DIM(1) = NJ
      IJ_DIM(2) = NI
      IJ_SYM(1) = JSM
      IJ_SYM(2) = ISM
      IJ_TYP(1) = JTYP
      IJ_TYP(2) = ITYP
    end if

    ! Generate creation- or annihilation- mappings for all Ka strings

    ! For operator connecting to |Ka> and |Ja> i.e. operator 2
    call ADAST_GAS(IJ_SYM(2),IJ_TYP(2),NGAS,JASPGP,JASM,I1,XI1S,NKASTR,KACT,SIGNIJ2,IJAC)
    !call ADAST_GAS(JSM,JTYP,JATP,JASM,IAGRP,I1,XI1S,NKASTR,KACT,SCLFACS,IJ_AC)
    ! For operator connecting |Ka> and |Ia>, i.e. operator 1
    call ADAST_GAS(IJ_SYM(1),IJ_TYP(1),NGAS,IASPGP,IASM,I3,XI3S,NKASTR,KACT,One,IJAC)
    !call ADAST_GAS(ISM,ITYP,NGAS,IASPGP,IASM,I3,XI3S,NKASTR,KACT,One,IJ_AC)
    ! Compress list to common nonvanishing elements
    IDOCOMP = 0
    if (IDOCOMP == 1) then
      call COMPRS2LST(I1,XI1S,IJ_DIM(2),I3,XI3S,IJ_DIM(1),NKASTR,NKAEFF)
    else
      NKAEFF = NKASTR
    end if

    ! Loop over batches of KA strings
    NKABTC = 0
    do
      NKABTC = NKABTC+NPROCS
      NKABTCSZ = max(NKAEFF-1,0)/NKABTC+1
      if (NKABTCSZ <= MAXK) exit
    end do

    do IKABTC=1+MYRANK,NKABTC,NPROCS
      KABOT = (IKABTC-1)*NKABTCSZ+1
      KATOP = min(KABOT+NKABTCSZ-1,NKAEFF)
      LKABTC = KATOP-KABOT+1
      if (LKABTC <= 0) exit
      ! Obtain C(ka,J,JB) for Ka in batch
      call TIMING(CPU0,CPU,WALL0,WALL)
      do JJ=1,IJ_DIM(2)
        call GET_CKAJJB(CB,IJ_DIM(2),NJA,CJRES,LKABTC,NJB,JJ,I1(KABOT+(JJ-1)*NKASTR),XI1S(KABOT+(JJ-1)*NKASTR))
      end do
      call TIMING(CPU1,CPU,WALL1,WALL)
      TSIGMA(4) = TSIGMA(4)+(WALL1-WALL0)
#     ifdef _DEBUGPRINT_
      write(u6,*) ' Updated CJRES as C(Kaj,Jb)'
      call WRTMAT(CJRES,NKASTR*NJ,NJB,NKASTR*NJ,NJB)
#     endif

      !MXACJ = MAX(MXACJ,NIB*LKABTC*IJ_DIM(1),NJB*LKABTC*IJ_DIM(2))
      SIRES(1:NIB*LKABTC*IJ_DIM(1)) = Zero
      FACS = One

      do KLTYP=1,NKLTYP
        KTYP = KTP(KLTYP)
        LTYP = LTP(KLTYP)
        ! Allowed double excitation ?
        !IJKL_ACT = 1
        !if (IJKL_ACT == 0) cycle
#       ifdef _DEBUGPRINT_
        write(u6,*) ' KTYP, LTYP',KTYP,LTYP
#       endif
        ! Should this group of excitations be included
        !if (NSEL2E /= 0) then
        !  IAMOKAY = 0
        !  if ((ITYP == JTYP) .and. (ITYP == KTYP) .and. (ITYP == LTYP)) then
        !    do JSEL2E=1,NSEL2E
        !      if (ISEL2E(JSEL2E) == ITYP) IAMOKAY = 1
        !    end do
        !  end if
        !  if (IAMOKAY == 0) cycle
        !end if

        KL_TYP(1) = KTYP
        KL_TYP(2) = LTYP
        !KL_AC(1) = 2
        !KL_AC(2) = 1
        !NOP = 2
        !if (IUSE_PH == 1) then
        !  call ALG_ROUTERX(IBOC,JBOC,NOP,KL_TYP,KL_AC,KL_REO,SIGNKL)
        !else
        ! Enforced a+ a
        KL_REO(1) = 1
        KL_REO(2) = 2
        SIGNKL = One
        !end if

        do KSM=1,NSMOB
          LSM = Mul(KSM,KLSM)
#         ifdef _DEBUGPRINT_
          write(u6,*) ' KSM, LSM',KSM,LSM
#         endif
          if (LSM == 0) cycle
          NK = NOBPTS(KTYP,KSM)
          NL = NOBPTS(LTYP,LSM)

          if ((KL_REO(1) == 1) .and. (KL_REO(2) == 2)) then
            ! Business as usual i.e. creation map
            KLAC = 2
            KL_DIM(1) = NK
            KL_DIM(2) = NL
            KL_SYM(1) = KSM
            KL_SYM(2) = LSM
            KL_TYP(1) = KTYP
            KL_TYP(2) = LTYP
          else
            ! Terra Nova, annihilation map
            KLAC = 1
            KL_DIM(1) = NL
            KL_DIM(2) = NK
            KL_SYM(1) = LSM
            KL_SYM(2) = KSM
            KL_TYP(1) = LTYP
            KL_TYP(2) = KTYP
          end if
          ! If IUSEAB is used, only terms with i >= k will be generated so
          IKORD = 0
          if ((IUSEAB == 1) .and. (ISM > KSM)) cycle
          if ((IUSEAB == 1) .and. (ISM == KSM) .and. (ITYP < KTYP)) cycle
          if ((IUSEAB == 1) .and. (ISM == KSM) .and. (ITYP == KTYP)) IKORD = 1

          if ((NK == 0) .or. (NL == 0)) cycle
          ! Obtain all connections a+l!Kb> = +/-/0!Jb>
          ! currently we are using creation mappings for kl
          call ADAST_GAS(KL_SYM(2),KL_TYP(2),NGAS,JBSPGP,JBSM,I2,XI2S,NKBSTR,KACT,SIGNKL,KLAC)
          if (NKBSTR == 0) cycle
          ! Obtain all connections a+k!Kb> = +/-/0!Ib>
          call ADAST_GAS(KL_SYM(1),KL_TYP(1),NGAS,IBSPGP,IBSM,I4,XI4S,NKBSTR,KACT,One,KLAC)
          if (NKBSTR == 0) cycle

          ! Fetch Integrals as (iop2 iop1 | k l)

          IXCHNG = 0
          ICOUL = 1
          ! Normal integrals with conjugation symmetry
          call GETINT(XINT,IJ_TYP(2),IJ_SYM(2),IJ_TYP(1),IJ_SYM(1),KL_TYP(1),KL_SYM(1),KL_TYP(2),KL_SYM(2),IXCHNG,0,0,ICOUL)

          ! S(Ka,i,Ib) = sum(j,k,l,Jb)<Ib!a+kba lb!Jb>C(Ka,j,Jb)*(ji!kl)

          IROUTE = 3
          call TIMING(CPU0,CPU,WALL0,WALL)
          call SKICKJ(SIRES,CJRES,LKABTC,NKBSTR,XINT,IJ_DIM(1),IJ_DIM(2),KL_DIM(1),KL_DIM(2),NKBSTR,I4,XI4S,I2,XI2S,IKORD,FACS, &
                      IROUTE)
          call TIMING(CPU1,CPU,WALL1,WALL)
          TSIGMA(5) = TSIGMA(5)+(WALL1-WALL0)

#         ifdef _DEBUGPRINT_
          write(u6,*) ' Updated Sires as S(Kai,Ib)'
          call WRTMAT(SIRES,LKABTC*NI,NIB,LKABTC*NI,NIB)
#         endif

        end do
        ! End of loop over KSM
      end do
      ! End of loop over KLTYP

      ! Scatter out from s(Ka,Ib,i)

#     ifdef _DEBUGPRINT_
      write(u6,*) ' S(Ka,Ib,i) as S(Ka,Ibi)'
      call WRTMAT(SIRES,LKABTC,NIB*IJ_DIM(1),LKABTC,IJ_DIM(1))
#     endif

      call TIMING(CPU0,CPU,WALL0,WALL)
      do II=1,IJ_DIM(1)
        call ADD_SKAIIB(SB,IJ_DIM(1),NIA,SIRES,LKABTC,NIB,II,I3(KABOT+(II-1)*NKASTR),XI3S(KABOT+(II-1)*NKASTR))
      end do
      call TIMING(CPU1,CPU,WALL1,WALL)
      TSIGMA(6) = TSIGMA(6)+(WALL1-WALL0)
    end do
    ! End of loop over partitioning of alpha strings
  end do
  ! End of loop over ISM
end do
! End of loop over IJTYP

end subroutine RSBB2BN
