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
! Copyright (C) 1991, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine RSSBCB2(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,NGAS,IAOC,IBOC,JAOC,JBOC,NAEL,NBEL,IJAGRP,IJBGRP,SB,CB,JDOH2,NOBPTS, &
                   MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,NSMOB,NSMST,NIA,NIB,NJA,NJB,IDC,CJRES,SIRES,I3,XI3S,I4,XI4S,MOCAA, &
                   SCLFAC,IPHGAS,I_RES_AB)
! SUBROUTINE RSSBCB2 --> 82
!
! Contributions to sigma block (iasm iatp, ibsm ibtp) from
! C block (jasm jatp, jbsm, jbtp)
!
! =====
! Input
! =====
!
! IASM,IATP : Symmetry and type of alpha strings in sigma
! IBSM,IBTP : Symmetry and type of beta  strings in sigma
! JASM,JATP : Symmetry and type of alpha strings in C
! JBSM,JBTP : Symmetry and type of beta  strings in C
! NGAS      : Number of active spaces in calculation
! IAOC,IBOC : Number of electrons in each AS for sigma supergroups
! JAOC,JBOC : Number of electrons in each AS for C     supergroups
! NAEL : Number of alpha electrons
! NBEL : Number of  beta electrons
! IJAGRP    : IA and JA belongs to this group of strings
! IJBGRP    : IB and JB belongs to this group of strings
! CB : Input c block
! IDOH2 : = 0 => no two electron operator
! IDOH2 : = 1 =>    two electron operator
! NTSOB  : Number of orbitals per type and symmetry
! IBTSOB : base for orbitals of given type and symmetry
! IBORB  : Orbitals of given type and symmetry
! MAXI   : Largest Number of "spectator strings" treated simultaneously
! MAXK   : Largest number of inner resolution strings treated at simult.
!
! IHAPR : /= 0 implies thatt the exact Hamiltonian shoulf not be uses
! In the case IPTSPC and JPTSPC defined the perturbation spaces
! a nonvanishing perturbation is allowed inside each subspace.
! The actual type of approximate Hamiltonian in each subspace is defined by
! IHFORM
! NNSEL2E : Only selected 2e terms will be included
! ISEL2E : orbital spaces in which 2e terms are included
!          (Currently : all indices identical)
!
! ======
! Output
! ======
! SB : fresh sigma block
!
! =======
! Scratch
! =======
! SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
!              largest number of orbital pairs of given symmetries and
!              types.
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! I1, XI1S   : at least MXSTSO : Largest number of strings of given
!              type and symmetry
! C2 : Must hold largest STT block of sigma or C
!
! XINT : Scratch space for integrals.
!
! Jeppe Olsen, Winter of 1991

use lucia_data, only: TSIGMA
use Constants, only: Zero, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: IASM, IATP, IBSM, IBTP, JASM, JATP, JBSM, JBTP, NGAS, IAOC(NGAS), IBOC(NGAS), JAOC(NGAS), &
                                 JBOC(NGAS), NAEL, NBEL, IJAGRP, IJBGRP, JDOH2, NOBPTS(*), MAXI, MAXK, NSMOB, NSMST, NIA, NIB, &
                                 NJA, NJB, IDC, MOCAA, IPHGAS(*), I_RES_AB
real(kind=wp), intent(inout) :: SB(NIA*NIB), CB(NJA*NJB), XI1S(*), XI2S(*), XI3S(*), XI4S(*)
real(kind=wp), intent(_OUT_) :: SSCR(*), CSCR(*), C2(*), XINT(*), CJRES(*), SIRES(*)
integer(kind=iwp), intent(inout) :: I1(*), I2(*), I3(*), I4(*)
real(kind=wp), intent(in) :: SCLFAC
integer(kind=iwp) :: I12, IBLOCK(8), IDIAG, IDOH2, IIDC, IIITRNS, ITASK, IUSEAB, JJJTRNS, LADVICE
real(kind=wp) :: CPU, CPU0, CPU1, FACTOR, WALL, WALL0, WALL1

! For H(apr)

#ifdef _DEBUGPRINT_
write(u6,*) ' ==============================='
write(u6,*) ' RSSBCB2 :  C block (transposed)'
write(u6,*) ' ================================'
call WRTMAT(CB,NJB,NJA,NJB,NJA)
write(u6,*) ' ======================================'
write(u6,*) ' RSSBCB2 : Initial  S block(transposed)'
write(u6,*) ' ======================================'
call WRTMAT(SB,NIA,NIB,NIA,NIB)
write(u6,*) ' Overall scalefactor ',SCLFAC
write(u6,*) ' JDOH2 = ',JDOH2
write(u6,*) ' I_RES_AB = ',I_RES_AB

write(u6,*) ' IAOC and IBOC'
call IWRTMA(IAOC,1,NGAS,1,NGAS)
call IWRTMA(IBOC,1,NGAS,1,NGAS)
write(u6,*) ' JAOC and JBOC  :'
call IWRTMA(JAOC,1,NGAS,1,NGAS)
call IWRTMA(JBOC,1,NGAS,1,NGAS)
write(u6,*) ' IASM IATP JASM JATP ',IASM,IATP,JASM,JATP
write(u6,*) ' IBSM IBTP JBSM JBTP ',IBSM,IBTP,JBSM,JBTP
write(u6,*) ' NAEL NBEL ',NAEL,NBEL
#endif
! Should the corresponding Hamiltonian matrix block be
! calculated exactly or approximately
!if (IHAPR /= 0) then
!  call HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,IMZERO,IDIAG)
!# ifdef _DEBUGPRINT_
!  write(u6,*) ' RSSBCBN : ',NNSEL2E,ISEL2E(1)
!# endif
!  NSEL2E = NNSEL2E
!  if (IMZERO /= 0) goto 9999
!else
! Operator specified by input
IDOH2 = JDOH2
IDIAG = 0
!end if
#ifdef _DEBUGPRINT_
write(u6,*) ' IDIAG IDOH2 ',IDIAG,IDOH2
#endif

if ((IDC == 2) .and. (IATP == IBTP) .and. (IASM == IBSM) .and. (I_RES_AB == 0) .and. (JASM == JBSM) .and. (JATP == JBTP)) then
  ! Diagonal sigma block, use alpha-beta symmetry to reduce computations.
  IUSEAB = 1
else
  IUSEAB = 0
end if

if (IDIAG == 0) then

  ! Calculate block exactly

  if ((I_RES_AB /= 1) .and. (IUSEAB == 0) .and. (IATP == JATP) .and. (JASM == IASM)) then

    ! =============================
    ! Sigma beta beta contribution
    ! =============================

    ! Sigma aa(IA,IB) = sum(i > k,j > l)<IB!Eb(ij)Eb(kl)!JB>
    !                 * ((ij!kl)-(il!kj)) C(IA,JB)
    !                 + sum(ij) <IB!Eb(ij)!JB> H(ij) C(IA,JB)
    ! One electron part
    call TRPMT3(SB,NIB,NIA,C2)
    SB(1:NIA*NIB) = C2(1:NIA*NIB)
    call TRPMT3(CB,NJB,NJA,C2)
    CB(1:NJA*NJB) = C2(1:NJA*NJB)

    if (NBEL >= 0) then
#     ifdef _DEBUGPRINT_
      write(u6,*) ' SB before RSBB1E'
      call wrtmat(sb,nia,nib,nia,nib)
      write(u6,*) ' I am going to call RSBB1E'
#     endif
      call TIMING(CPU0,CPU,WALL0,WALL)
      call RSBB1E(IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NGAS,IBOC,JBOC,SB,CB,NOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,NSMOB,MOCAA, &
                  SCLFAC,IPHGAS)
      call TIMING(CPU1,CPU,WALL1,WALL)
      TSIGMA(1) = TSIGMA(1)+(WALL1-WALL0)

      ! CALL RSBB1E --> 33

#     ifdef _DEBUGPRINT_
      write(u6,*) ' SB after RSBB1E'
      call wrtmat(sb,nib,nia,nib,nia)
      write(u6,*) ' first element of SB after RSBB1E',SB(1)
#     endif

    end if
    if ((IDOH2 /= 0) .and. (NBEL >= 0)) then
      ! Two electron part
#     ifdef _DEBUGPRINT_
      write(u6,*) ' I am going to call RSBB2A'
#     endif
      call TIMING(CPU0,CPU,WALL0,WALL)
      call RSBB2A(IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NGAS,IBOC,JBOC,SB,CB,NOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,XINT,NSMOB,NSMST,SCLFAC, &
                  IPHGAS)
      call TIMING(CPU1,CPU,WALL1,WALL)
      TSIGMA(2) = TSIGMA(2)+(WALL1-WALL0)

      ! CALL RSBB2A --> 46

#     ifdef _DEBUGPRINT_
      write(u6,*) ' SB after RSBB2a'
      call wrtmat(sb,nib,nia,nib,nia)
      write(u6,*) ' first element of SB after RSBB1E',SB(1)
#     endif
    end if
    call TRPMT3(SB,NIA,NIB,C2)
    SB(1:NIA*NIB) = C2(1:NIA*NIB)
    call TRPMT3(CB,NJA,NJB,C2)
    CB(1:NJA*NJB) = C2(1:NJA*NJB)
  end if

  ! =============================
  ! Sigma alpha beta contribution
  ! =============================

  if ((IDOH2 /= 0) .and. (NAEL >= 0) .and. (NBEL >= 0)) then
#   ifdef _DEBUGPRINT_
    write(u6,*) ' I am going to call RSBB2B'
#   endif
    IIITRNS = 1
    if (IIITRNS == 1) then
      ! Call advice routine
      !    ADVICE_SIGMA(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,LADVICE)
      call ADVICE_SIGMA(IAOC,IBOC,JAOC,JBOC,LADVICE)
      ! LADVICE = 2 => implies transpose
      if (LADVICE == 2) then
        JJJTRNS = 1
      else
        JJJTRNS = 0
      end if
    end if

    !write(u6,*) ' IUSE_PA = ',IUSE_PA

    if (JJJTRNS == 0) then
      !if (IUSE_PA == 0) then
      call TIMING(CPU0,CPU,WALL0,WALL)
      call RSBB2BN(IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,SB,CB,NOBPTS, &
                   MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,NSMOB,IUSEAB,CJRES,SIRES,SCLFAC,IPHGAS)
      call TIMING(CPU1,CPU,WALL1,WALL)
      TSIGMA(3) = TSIGMA(3)+(WALL1-WALL0)

      ! CALL RSBB2BN --> 52

      !else if (IUSE_PA == 1) then
      !  !write(u6,*) ' RSBB2BVN will be called'
      !  call RSBB2BVN(IASM,IATP,IBSM,IBTP,NIA,NIB,JASM,JATP,JBSM,JBTP,NJA,NJB,IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,SB,CB, &
      !                MXPNGAS,NOBPTS,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,NSMOB,IUSEAB,CJRES,SIRES,SCLFAC,0,0,0, &
      !                IUSE_PH,IPHGAS,CJPA,SIPA)
      !end if

    else if (JJJTRNS == 1) then
      ! well lets give the transpose routine some more practice : Transpose back
      call TRPMT3(SB,NIB,NIA,C2)
      SB(1:NIA*NIB) = C2(1:NIA*NIB)
      call TRPMT3(CB,NJB,NJA,C2)
      CB(1:NJA*NJB) = C2(1:NJA*NJB)
      !write(u6,*) ' RSSBCB2 : Transpose path choosen'

      !if (IUSE_PA == 0) then
      ! No division into active/passive
      call TIMING(CPU0,CPU,WALL0,WALL)
      call RSBB2BN(IBSM,IBTP,IASM,IATP,NIB,NIA,JBSM,JBTP,JASM,JATP,NJB,NJA,IJBGRP,IJAGRP,NGAS,IBOC,IAOC,JBOC,JAOC,SB,CB,NOBPTS, &
                   MAXK,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,NSMOB,IUSEAB,CJRES,SIRES,SCLFAC,IPHGAS)
      call TIMING(CPU1,CPU,WALL1,WALL)
      TSIGMA(3) = TSIGMA(3)+(WALL1-WALL0)

      ! CALL RSBB2BN --> 52

      !else
      !  ! Divide into active/passive
      !  call RSBB2BVN(IBSM,IBTP,IASM,IATP,NIB,NIA,JBSM,JBTP,JASM,JATP,NJB,NJA,IJBGRP,IJAGRP,NGAS,IBOC,IAOC,JBOC,JAOC,SB,CB, &
      !                MXPNGAS,NOBPTS,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,XINT,NSMOB,IUSEAB,CJRES,SIRES,SCLFAC,0,0,0, &
      !                IUSE_PH,IPHGAS,CJPA,SIPA)
      !end if

      ! Transpose (To compensate later transposition)
      call TRPMT3(SB,NIA,NIB,C2)
      SB(1:NIA*NIB) = C2(1:NIA*NIB)
      call TRPMT3(CB,NJA,NJB,C2)
      CB(1:NJA*NJB) = C2(1:NJA*NJB)
    end if
#   ifdef _DEBUGPRINT_
    write(u6,*) ' SB after RSBB2B, first element'
    call wrtmat(sb,1,1,nia,nib)
    write(u6,*) ' SB after RSBB2b'
    call wrtmat(sb,nia,nib,nia,nib)
#   endif
  end if

  ! =============================
  ! Sigma alpha alpha contribution
  ! =============================

  if ((I_RES_AB /= -1) .and. (NAEL >= 0) .and. (IBTP == JBTP) .and. (IBSM == JBSM)) then

    ! alpha single excitation

#   ifdef _DEBUGPRINT_
    write(u6,*) ' I am going to call RSBB1E (last time)'
#   endif
    call TIMING(CPU0,CPU,WALL0,WALL)
    call RSBB1E(IASM,IATP,JASM,JATP,IJAGRP,NIB,NGAS,IAOC,JAOC,SB,CB,NOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,NSMOB,MOCAA, &
                SCLFAC,IPHGAS)
    call TIMING(CPU1,CPU,WALL1,WALL)
    TSIGMA(1) = TSIGMA(1)+(WALL1-WALL0)

    ! CALL RSBB1E --> 33

#   ifdef _DEBUGPRINT_
    write(u6,*) ' SB transposed after RSBB1, first element'
    call wrtmat(sb,1,1,nia,nib)
    write(u6,*) ' SB transposed  after RSBB1E'
    call wrtmat(SB,nib,nia,nib,nia)
#   endif

    ! alpha double excitation

    if ((IDOH2 /= 0) .and. (NAEL >= 0)) then
#     ifdef _DEBUGPRINT_
      write(u6,*) ' I am going to call RSBB2A (last time)'
#     endif
      call TIMING(CPU0,CPU,WALL0,WALL)
      call RSBB2A(IASM,IATP,JASM,JATP,IJAGRP,NIB,NGAS,IAOC,JAOC,SB,CB,NOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,XINT,NSMOB,NSMST,SCLFAC, &
                  IPHGAS)
      call TIMING(CPU1,CPU,WALL1,WALL)
      TSIGMA(2) = TSIGMA(2)+(WALL1-WALL0)

      ! CALL RSBB2A --> 46

    end if

#   ifdef _DEBUGPRINT_
    write(u6,*) ' SB transposed after RSBB2A, first element'
    call wrtmat(sb,1,1,nia,nib)
    write(u6,*) ' SB after RSBB2A'
    call wrtmat(sb,nia,nib,nia,nib)
#   endif
  end if

else if (IDIAG == 1) then

  ! Diagonal approxiation (IDIAG = 1)
  ! or complete orbital conserving part of Ham (IH_OCC_CONS = 1)

  IBLOCK(1) = IATP
  IBLOCK(2) = IBTP
  IBLOCK(3) = IASM
  IBLOCK(4) = IBSM
  IBLOCK(5) = 1
  IBLOCK(6) = 1
  if (IDOH2 == 0) then
    I12 = 1
  else
    I12 = 2
  end if
  !write(u6,*) ' IDOH2, I12 ', IDOH2,I12
  ITASK = 2
  FACTOR = Zero
  ! Well, we are not using transposed matrices here so
  call TRPMT3(CB,NJB,NJA,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)

  if ((IATP == JATP) .and. (IBTP == JBTP) .and. (IASM == JASM) .and. (IBSM == JBSM)) then
    !write(u6,*) ' DIATERM2_GAS will be called'
    C2(1:NJA*NJB) = CB(1:NJA*NJB)
    ! Input is in det basis
    IIDC = 1
    call DIATERM2_GAS(FACTOR,ITASK,C2,1,IBLOCK,1,I12,IIDC)
  else
    C2(1:NIA*NIB) = Zero
  end if
  ! Remaining occupation conserving operator
  !if (IH_OCC_CONS == 1) &
  !  call HCONFDIA_BBM(NAEL,NBEL,IJAGRP,IJBGRP,IASM,IATP,IAOC,NIA,IBSM,IBTP,IBOC,NIB,JASM,JATP,JAOC,NJA,JBSM,JBTP,JBOC,NJB,XINT, &
  !                    CB,C2)

  if (IUSEAB == 0) then
    FACTOR = SCLFAC
  else
    FACTOR = Half*SCLFAC
  end if
  !    MAT_P_MATT(A,B,NR,NC,COEF)
  call MAT_P_MATT(SB,C2,NIB,NIA,FACTOR)
  call TRPMT3(CB,NJA,NJB,C2)
  CB(1:NJA*NJB) = C2(1:NJA*NJB)
end if

!9999 continue
! Clean up
!if (IHAPR /= 0) call HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,IMZERO,IDIAG)

#ifdef _DEBUGPRINT_
write(u6,*) ' ==================================='
write(u6,*) ' RSSBCB : Final S block (transposed)'
write(u6,*) ' ==================================='
call WRTMAT(SB,NIB,NIA,NIB,NIA)
#endif

end subroutine RSSBCB2
