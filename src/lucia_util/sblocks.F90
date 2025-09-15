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
! Copyright (C) 1991,1999, Jeppe Olsen                                 *
!***********************************************************************

subroutine SBLOCKS(NSBLOCK,ISBLOCK,CB,SB,C2,ICOCOC,ICSM,NSSOA,NSSOB,NAEL,IAGRP,NBEL,IBGRP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,NSMST,NSMOB, &
                   NOBPTS,MXPNGAS,MAXK,MAXI,XINT,CSCR,SSCR,NGAS,NELFSPGP,IDC,I1,XI1S,I2,XI2S,IDOH2,ISTRFL,PS,LUC,ICJKAIB,CJRES, &
                   SIRES,I3,XI3S,I4,XI4S,MOCAA,LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,IRESTRICT,ICONSPA,ICONSPB,SCLFAC,IH0SPC, &
                   ICBAT_RES,ICBAT_INI,ICBAT_END,IPHGAS,I_RES_AB)
! SUBROUTINE SBLOCKS --> 91
!
! Direct RAS routine employing combined MOC/n-1 resolution method
!
! Jeppe Olsen, Winter of 1991
!              Last modification : April 99
!
! =====
! Input
! =====
!
! NSBLOCK : Number of BLOCKS included
! ISBLOCK : Blocks included
!   ISBLOCK(1,*) : alpha type of block
!   ISBLOCK(2,*) : beta type of block
!   ISBLOCK(3,*) : sym of alpha in block
!   ISBLOCK(4,*) : Offset of block
!
! ICOCOC : Allowed type combinations for C
! ICSM   : Symmetry for C
! NACOB  : Number of active orbitals
! NSSOA  : Number of strings per type and symmetry for alpha strings
! NAEL   : Number of active alpha electrons
! NSSOB  : Number of strings per type and symmetry for beta strings
! NBEL   : Number of active beta electrons
! NTSOB  : Number of orbitals per type and symmetry
! NOBPTS : Orbitals of given type and symmetry
!
! MAXIJ  : Largest allowed number of orbital pairs treated simultaneously
! MAXK   : Largest number of N-2,N-1 strings treated simultaneously
! MAXI   : Max number of N strings treated simultaneously
!
! LI     : Length of scratch array for integrals
! LS     : Length of scratch array for S
! XINT   : Scratch array for integrals
! CSCR   : Scratch array for C vector
! SSCR   : Scratch array for S vector
!
! ICJKAIB = 1 => construct C(Ka,Jb,j) and S(Ka,IB,i) as intermediate terms
!         = 0 => do not construct the above montioned matrices
! CJRES,SIRES : Space for above matrices
! The C and S vectors are accessed through routines that
! either fetches/disposes symmetry blocks or
! Symmetry-occupation-occupation blocks
!
! If IRESTRICT /= 0 THEN we are after :
! sigma(iblk) = summa(jblk <= iblk) (2-delta(iblk,jblk))/2
!                                                 * <Iblk!H!Jblk>C(Jblk)

use lucia_data, only: IDISK
use spinfo, only: DOBKAP
#ifdef _DEBUGPRINT_
use spinfo, only: NGASBK, IOCCPSPC
#endif
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NSBLOCK, ISBLOCK(8,*), NOCTPA, NOCTPB, ICOCOC(NOCTPA,NOCTPB), NSMST, ICSM, NSSOA(NSMST,*), &
                                 NSSOB(NSMST,*), NAEL, IAGRP, NBEL, IBGRP, NSMOB, MXPNGAS, NOBPTS(MXPNGAS,*), MAXK, MAXI, NGAS, &
                                 NELFSPGP(MXPNGAS,*), IDC, IDOH2, ISTRFL(*), LUC, ICJKAIB, MOCAA, IRESTRICT, &
                                 ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB), IH0SPC(NOCTPA,NOCTPB), ICBAT_RES, ICBAT_INI, &
                                 ICBAT_END, IPHGAS(*), I_RES_AB
real(kind=wp), intent(inout) :: CB(*), SB(*), XI1S(*), XI2S(*), XI3S(*), XI4S(*)
real(kind=wp), intent(_OUT_) :: C2(*), XINT(*), CSCR(*), SSCR(*), CJRES(*), SIRES(*), SCLFAC(*)
integer(kind=iwp), intent(out) :: IOCTPA, IOCTPB
integer(kind=iwp), intent(inout) :: I1(*), I2(*), I3(*), I4(*)
real(kind=wp), intent(in) :: PS
integer(kind=iwp), intent(_OUT_) :: LCBLOCK(*), LECBLOCK(*), I1CBLOCK(*), ICBLOCK(8,*)
integer(kind=iwp) :: I_DO_EXACT_BLK, IASM, IATP, IBSM, IBTP, ICBLK, ICOFF, ICOOSC(1), iDUMMY(1), INTERACT, IOFF, IPERM, IPTSPC, &
                     ISBLK, ISCALE, ISOFF, JASM, JATP, JBLOCK, JBSM, JBTP, JCBAT_END, JCBAT_INI, JCBATCH, JJCBLOCK, JOFF, JPTSPC, &
                     JSBLOCK, LASM(4), LATP(4), LBL, LBSM(4), LBTP(4), LLASM, LLATP, LLBSM, LLBTP, LSGN(5), LTRP(5), MXEXC, NASTR, &
                     NBSTR, NCBATCH, NIA, NIB, NJA, NJB, NJBLOCK, NLLA, NLLB, NPERM
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IBLOCK, IGAS, II
#endif
real(kind=wp) :: C(1), FACTOR, PL, XFAC
! IH_OCC_CONS = 1 implies that we should employ occupation conserving part of Hamiltonian
integer(kind=iwp), parameter :: IH_OCC_CONS = 0

#ifdef _DEBUGPRINT_
write(u6,*) ' ================='
write(u6,*) ' SBLOCKS speaking :'
write(u6,*) ' ================='
write(u6,*)
write(u6,*) ' Number of sigma blocks to be calculated ',NSBLOCK
write(u6,*) ' TTSS for each ACTIVE sigma block'
do IBLOCK=1,NSBLOCK
  if (ISBLOCK(1,IBLOCK) > 0) write(u6,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
end do
write(u6,*) ' IDC PS',IDC,PS
write(u6,*) ' IDOH2 = ',IDOH2
write(u6,*) ' I_RES_AB=',I_RES_AB
if (DoBKAP) then
  write(u6,*) ' I am doing BK-type of approximation'
  write(u6,*) ' It is based on orbital splitting'
  write(u6,*) ' Min and Max for subspace with exact Hamiltonian'
  write(u6,*) ' ==============================================='
  write(u6,*) 'NGASBK : ',NGASBK
  write(u6,*) '              Min. Occ.      Max. Occ.'
  do IGAS=1,NGASBK
    write(u6,'(A,I2,10X,I3,9X,I3)') '   GAS',IGAS,IOCCPSPC(IGAS,1),IOCCPSPC(IGAS,2)
  end do
end if

write(u6,*) ' Initial C vector'
call WRTVCD(CB,LUC,1,-1)
#endif
! ==========================
! 1 : Arrays for accessing C
! ==========================
! Find batches of C - strings
call PART_CIV2(IDC,NSSOA,NSSOB,NOCTPA,NOCTPB,NSMST,ICOCOC,ICSM,NCBATCH,LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,0)
! Find the active blocks on LUC, store info in SCLFAC
call FIND_ACTIVE_BLOCKS(LUC,-1,SCLFAC,CB)

! Initialize sigma blocks
do JSBLOCK=1,NSBLOCK
  IATP = ISBLOCK(1,JSBLOCK)
  IBTP = ISBLOCK(2,JSBLOCK)
  IASM = ISBLOCK(3,JSBLOCK)
  IBSM = ISBLOCK(4,JSBLOCK)
  IOFF = ISBLOCK(5,JSBLOCK)
  NASTR = NSSOA(IASM,IATP)
  NBSTR = NSSOB(IBSM,IBTP)
  if (ISBLOCK(1,JSBLOCK) > 0) SB(IOFF:IOFF+NASTR*NBSTR-1) = Zero
end do
! Loop over batches over C blocks
if (IDOH2 == 1) then
  MXEXC = 2
else
  MXEXC = 1
end if
IDISK(LUC) = 0
if (ICBAT_RES == 1) then
  write(u6,*) ' Restricted set of C batches'
  write(u6,*) ' ICBAT_INI ICBAT_END',ICBAT_INI,ICBAT_END
  JCBAT_INI = ICBAT_INI
  JCBAT_END = ICBAT_END
else
  JCBAT_INI = 1
  JCBAT_END = NCBATCH
end if

JOFF = 0 ! jwk-cleanup
do JCBATCH=JCBAT_INI,JCBAT_END

  ! Read C blocks into core

  ICOFF = 1
  NJBLOCK = LCBLOCK(JCBATCH)
  do JJCBLOCK=1,NJBLOCK
    JBLOCK = I1CBLOCK(JCBATCH)-1+JJCBLOCK
    ! Will this block be needed ??
    INTERACT = 0
    if (SCLFAC(JBLOCK) == One) then
      JATP = ICBLOCK(1,JBLOCK)
      JBTP = ICBLOCK(2,JBLOCK)
      JASM = ICBLOCK(3,JBLOCK)
      JBSM = ICBLOCK(4,JBLOCK)
      JOFF = ICBLOCK(5,JBLOCK)
      PL = One
      call PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
      do IPERM=1,NPERM
        LLASM = LASM(IPERM)
        LLBSM = LBSM(IPERM)
        LLATP = LATP(IPERM)
        LLBTP = LBTP(IPERM)
        ! Loop over Sigma blocks in batch
        do JSBLOCK=1,NSBLOCK
          if (ISBLOCK(1,JSBLOCK) > 0) then
            IATP = ISBLOCK(1,JSBLOCK)
            IBTP = ISBLOCK(2,JSBLOCK)
            IASM = ISBLOCK(3,JSBLOCK)
            IBSM = ISBLOCK(4,JSBLOCK)
            ! Are the two blocks connected by allowed excitation
            call CON_BLOCKS(IATP,IBTP,LLATP,LLBTP,IASM,IBSM,LLASM,LLBSM,ICONSPA,ICONSPB,NOCTPA,NOCTPB,MXEXC,IH_OCC_CONS,INTERACT)

          end if
        end do
      end do
      ! End of checking whether C-block is needed
    end if
    ! Checking was only done for nonvanishing blocks

    ISCALE = 0
    if (INTERACT == 1) then
      call GSTTBL(C,CB(JOFF),JATP,JASM,JBTP,JBSM,NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,PL,LUC,C2,NSMST,ISCALE,SCLFAC(JBLOCK))
      ! Note in GSTTBL : ICOOSC only used for CI vectors in core,
    else
      ! not relevant
      call IDAFILE(LUC,2,iDUMMY,1,IDISK(LUC))
      LBL = iDUMMY(1)
      call IDAFILE(LUC,2,iDUMMY,1,IDISK(LUC))
      call SKPRCD2(LBL,-1,LUC)
      SCLFAC(JBLOCK) = Zero
    end if

#   ifdef _DEBUGPRINT_
    if (INTERACT == 1) then
      write(u6,*) ' TTSS for C block read in'
      call IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
    else
      write(u6,*) ' TTSS for C block skipped'
      call IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
    end if
#   endif

  end do
  ! End of loop over Blocks

  ! Loop over blocks of sigma and C in core and obtain  contribution from
  ! given C block to given S block
  ! Loop over C blocks
  do ICBLK=I1CBLOCK(JCBATCH),I1CBLOCK(JCBATCH)-1+NJBLOCK
    JATP = ICBLOCK(1,ICBLK)
    JBTP = ICBLOCK(2,ICBLK)
    JASM = ICBLOCK(3,ICBLK)
    JBSM = ICBLOCK(4,ICBLK)
    ICOFF = ICBLOCK(5,ICBLK)
    NJA = NSSOA(JASM,JATP)
    NJB = NSSOB(JBSM,JBTP)

    if (SCLFAC(ICBLK) /= Zero) then
      ! Other symmetry blocks that can be obtained from this block
      call PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
      ! Start with transposed block
      do IPERM=NPERM,1,-1
        LLASM = LASM(IPERM)
        LLBSM = LBSM(IPERM)
        LLATP = LATP(IPERM)
        LLBTP = LBTP(IPERM)
        NLLA = NSSOA(LLASM,LLATP)
        NLLB = NSSOB(LLBSM,LLBTP)
        ! The routines assumes on input that the blocks are transposed, so,
        ! Initial tour, IPERM = 1 corresponds always to no transpose, so transpose!
        if (IPERM == 1) then
          if ((IDC == 2) .and. (JATP == JBTP) .and. (JASM == JBSM)) then
            ! Diagonal blocks, Transposing corresponds to scaling
            if (PS == -One) CB(ICOFF:ICOFF+NJA*NJB-1) = -CB(ICOFF:ICOFF+NJA*NJB-1)
          else
            ! off-diagonal blocks, explicit transposing
            call TRPMT3(CB(ICOFF),NJA,NJB,C2)
            CB(ICOFF:ICOFF+NJA*NJB-1) = C2(1:NJA*NJB)
          end if
        end if

        do ISBLK=1,NSBLOCK
          if (ISBLOCK(1,ISBLK) <= 0) cycle
          IATP = ISBLOCK(1,ISBLK)
          IBTP = ISBLOCK(2,ISBLK)
          IASM = ISBLOCK(3,ISBLK)
          IBSM = ISBLOCK(4,ISBLK)
          ISOFF = ISBLOCK(5,ISBLK)
          NIA = NSSOA(IASM,IATP)
          NIB = NSSOB(IBSM,IBTP)

          if (NIA*NIB == 0) cycle
          if ((IRESTRICT == 1) .and. &
              ((JASM > IASM) .or. (JASM == IASM) .and. (JATP > IATP) .or. (JASM == IASM) .and. (JATP == IATP) .and. &
               (JBTP > IBTP))) cycle
          ! Are the two blocks connected by allowed excitation
          call CON_BLOCKS(IATP,IBTP,LLATP,LLBTP,IASM,IBSM,LLASM,LLBSM,ICONSPA,ICONSPB,NOCTPA,NOCTPB,MXEXC,IH_OCC_CONS,INTERACT)

          ! IF BK approximation is active, check whether block should
          ! be calculated exactly (1), by diagonal (-1) or is set to zero (0).
          I_DO_EXACT_BLK = 1
          if (DoBKAP) call CHECK_BLOCKS_FOR_BK_APPROX(IATP,IBTP,LLATP,LLBTP,IASM,IBSM,LLASM,LLBSM,IOCTPA,IOCTPB,I_DO_EXACT_BLK)
          ! BK-like approximation stuff
          if ((INTERACT == 0) .or. (I_DO_EXACT_BLK == 0)) cycle

#         ifdef _DEBUGPRINT_
          write(u6,*) ' Next s block in batch :'
          write(u6,*) ' ISBLK IASM IBSM IATP IBTP'
          write(u6,'(5I5)') ISBLK,IASM,IBSM,IATP,IBTP
#         endif

          if ((IDC == 2) .and. (IASM == IBSM) .and. (IATP == IBTP) .and. &
              ((LLBSM > LLASM) .or. (LLASM == LLBSM) .and. (LLBTP > LLATP))) cycle

#         ifdef _DEBUGPRINT_
          write(u6,*) ' RSSBCB will be called for'
          write(u6,*) ' Sigma block :'
          write(u6,*) ' ISOFF ',ISOFF
          write(u6,*) ' ISBLK IASM IBSM IATP IBTP'
          write(u6,'(5I5)') ISBLK,IASM,IBSM,IATP,IBTP
          write(u6,*) ' C     block :'
          write(u6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
          write(u6,'(5I5)') ICBLK,LLASM,LLBSM,LLATP,LLBTP
          write(u6,*) ' ICOFF ',ICOFF
          write(u6,*) ' Overall scale',SCLFAC(ICBLK)
#         endif

          if ((IRESTRICT == 1) .and. &
              (((IASM == LLASM) .and. (IBSM == LLBSM) .and. (IATP == LLATP) .and. (IBTP == LLBTP)) .or. &
               ((IDC == 2) .and. (IASM == LLBSM) .and. (IBSM == LLASM) .and. (IATP == LLBTP) .and. (IBTP == LLATP)))) then
            XFAC = Half*SCLFAC(ICBLK)
          else
            XFAC = SCLFAC(ICBLK)
          end if
          ! Form of operator in action
          !if (IPERTOP /= 0) then
          ! Not exact Hamiltonian in use
          IPTSPC = IH0SPC(IATP,IBTP)
          JPTSPC = IH0SPC(JATP,JBTP)

          if (IPTSPC /= JPTSPC) cycle
          ! BK-like approximation stuff
          if (I_DO_EXACT_BLK == 1) then
            call RSSBCB2(IASM,IATP,IBSM,IBTP,LLASM,LLATP,LLBSM,LLBTP,NGAS,NELFSPGP(:,IATP+IOCTPA-1),NELFSPGP(:,IBTP+IOCTPB-1), &
                         NELFSPGP(:,LLATP+IOCTPA-1),NELFSPGP(:,LLBTP+IOCTPB-1),NAEL,NBEL,IAGRP,IBGRP,SB(ISOFF),CB(ICOFF),IDOH2, &
                         NOBPTS,MAXI,MAXK,SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,NSMOB,NSMST,NIA,NIB,NLLA,NLLB,IDC,CJRES,SIRES,I3,XI3S, &
                         I4,XI4S,MOCAA,XFAC,IPHGAS,I_RES_AB)
          else if (I_DO_EXACT_BLK == -1) then
            ! Giovanni.... transposing sigma and CI vectors:
            call TRPMT3(SB(ISOFF),NIB,NIA,C2)
            SB(ISOFF:ISOFF+NIA*NIB-1) = C2(1:NIA*NIB)
            call TRPMT3(CB(ICOFF),NLLB,NLLA,C2)
            CB(ICOFF:ICOFF+NLLA*NLLB-1) = C2(1:NLLA*NLLB)
            FACTOR = Zero
            call ADDDIA_TERM(FACTOR,CB(ICOFF),SB(ISOFF),IATP,IBTP,IASM,IBSM)
            ! Giovanni.... transposing back sigma and CI vectors:
            call TRPMT3(SB(ISOFF),NIA,NIB,C2)
            SB(ISOFF:ISOFF+NIA*NIB-1) = C2(1:NIA*NIB)
            call TRPMT3(CB(ICOFF),NLLA,NLLB,C2)
            CB(ICOFF:ICOFF+NLLA*NLLB-1) = C2(1:NLLA*NLLB)
          end if ! End BK stuff
          ! CALL RSSBCB2 --> 82
        end do
        ! End of loop over sigma blocks
      end do
    end if
    ! End of C-block is nonvanishing
  end do
  ! End of loop over C blocks in Batch
end do
! End of loop over batches of C blocks

! Order
do ISBLK=1,NSBLOCK
  if (ISBLOCK(1,ISBLK) > 0) then
    IATP = ISBLOCK(1,ISBLK)
    IBTP = ISBLOCK(2,ISBLK)
    IASM = ISBLOCK(3,ISBLK)
    IBSM = ISBLOCK(4,ISBLK)
    ISOFF = ISBLOCK(5,ISBLK)
    NIA = NSSOA(IASM,IATP)
    NIB = NSSOB(IBSM,IBTP)
    if (ICJKAIB /= 0) then
      ! Tranpose sigma block was obtained, transpose to obtain correct block
      call TRPMT3(SB(ISOFF),NSSOB(IBSM,IBTP),NSSOA(IASM,IATP),C2)
      SB(ISOFF:ISOFF+NSSOA(IASM,IATP)*NSSOB(IBSM,IBTP)-1) = C2(1:NSSOA(IASM,IATP)*NSSOB(IBSM,IBTP))
    end if
    if ((IDC == 2) .and. (IASM == IBSM) .and. (IATP == IBTP)) call TRPAD3(SB(ISOFF),PS,NSSOA(IASM,IATP))

  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' output blocks from SBLOCKS'
call WRTTTS(SB,ISBLOCK,NSBLOCK,NSMST,NSSOA,NSSOB,1)
#endif

end subroutine SBLOCKS
