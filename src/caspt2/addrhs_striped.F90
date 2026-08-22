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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************

#include "compiler_features.h"
#ifdef _MOLCAS_MPP_

module ADDRHS_STRIPED

! Striped counterparts of the ADDRHS routines, selected with PRHS = 4 or STRIPED (iParRHS = 4).
!
! The public routines take the same leading arguments as their replicated counterparts,
! with the work array and its size in place of the case block and the two scatter buffers,
! so that PROCESS_RHS_BLOCK can dispatch to either from the same variables.

#include "macros.fh"

use Constants, only: Zero, One, Two, Three, Half, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
private

integer(kind=iwp) :: iOffRHSLoc(8,13) = 0         ! offset of a block in it, 0 if not in memory
integer(kind=iwp) :: iRHSLocGrp = 0               ! case group currently in memory, 0 if none was ever loaded
integer(kind=iwp) :: iRHSLocSym = 0               ! symmetry of that group, only meaningful at tier 2
integer(kind=iwp) :: iRHSLocTier = 0              ! how much stays in memory: 0 all blocks, 1 one case group, 2 one symmetry of one
integer(kind=iwp) :: ISYMTGT = 0                  ! the only symmetry the kernels accumulate at tier 2, 0 when none is skipped
integer(kind=iwp) :: jLoRHSLoc(8,13) = 0          ! lower bound of the block owned by this process
integer(kind=iwp) :: jHiRHSLoc(8,13) = 0          ! upper bound of the block owned by this process
integer(kind=iwp) :: nRHSLocSz(8,13) = 0          ! length of that stripe, 0 if this process owns none
logical(kind=iwp) :: RHSLocOnDisk(8,13) = .false. ! true once the block has been flushed, so it may be read back

real(kind=wp), allocatable, target :: RHSLoc(:)   ! the buffer, the in-core blocks concatenated

integer(kind=iwp), parameter :: IGRP_A = 1, IGRP_B = 2, IGRP_E = 3, IGRP_H = 4, IGRP_D1 = 5, &
                                IGRP_G = 6, IGRP_F = 7, IGRP_C = 8, IGRP_D2 = 9, NRHSGRP = 9
integer(kind=iwp), parameter :: IRHSGRP(2,NRHSGRP) = reshape([ 1, 0, & ! A
                                                               2, 3, & ! B+ B-
                                                               6, 7, & ! E+ E-
                                                              12,13, & ! H+ H-
                                                               5, 0, & ! D1, the upper half of D
                                                              10,11, & ! G+ G-
                                                               8, 9, & ! F+ F-
                                                               4, 0, & ! C
                                                               5, 0],& ! D2, the lower half of D
                                                               [2,NRHSGRP])

! a-block width of ADDRHSF/G_STRIPED_D, 32 measured best
integer(kind=iwp), parameter :: NAMXCAP = 32

! weights of the symmetrized combinations, as in the replicated MKRHS routines
real(kind=wp), parameter :: SQ2 = sqrt(Two), SQ3 = sqrt(Three), SQH = sqrt(Half), SQ32 = sqrt(OneHalf)

public :: RHSLOC_SIZES, RHSLOC_ALLOCATE, RHSLOC_LOAD, RHSLOC_FINALIZE, RHSLOC_FREE
public :: ADDRHSA_STRIPED, ADDRHSB_STRIPED, ADDRHSC_STRIPED, ADDRHSD1_STRIPED, ADDRHSD2_STRIPED, &
          ADDRHSE_STRIPED, ADDRHSF_STRIPED, ADDRHSG_STRIPED, ADDRHSH_STRIPED
public :: IGRP_A, IGRP_B, IGRP_C, IGRP_D1, IGRP_D2, IGRP_E, IGRP_F, IGRP_G, IGRP_H

contains

!-----------------------------------------------------------------------
!
! The first part manages the local RHS buffer
!
!-----------------------------------------------------------------------

subroutine RHSLOC_SIZES(NALL,NGRP,NSYM1)

  use caspt2_module, only: NASUP, NISUP, NSYM

  integer(kind=iwp), intent(out) :: NALL, NGRP, NSYM1

  integer(kind=iwp) :: ICASE, IGRP, IHI, ILO, IPAIR, ISYM, JHI, JLO, NAS, NBLK, NIS

  ! Sizes of the local RHS:
  !   NALL for all blocks at once
  !   NGRP for the largest case group
  !   NSYM1 for the largest single symmetry of one, i.e. the buffer sizes of tiers 0, 1 and 2
  ! Also sets the column ranges and block sizes that the kernels index with.

  NALL = 0
  do ICASE=1,13
    do ISYM=1,NSYM
      iOffRHSLoc(ISYM,ICASE) = 0
      jLoRHSLoc(ISYM,ICASE) = 1
      jHiRHSLoc(ISYM,ICASE) = 0
      nRHSLocSz(ISYM,ICASE) = 0
      RHSLocOnDisk(ISYM,ICASE) = .false.
      NAS = NASUP(ISYM,ICASE)
      NIS = NISUP(ISYM,ICASE)
      if (NAS*NIS == 0) cycle
      ! the same column distribution the RHS global arrays get, from GA_CREATE_STRIPED
      call RHS_DISTRIBUTION(NAS,NIS,ILO,IHI,JLO,JHI)
      if (JHI < JLO) cycle ! this process owns no column of this block, leave the 1:0 set above
      jLoRHSLoc(ISYM,ICASE) = JLO
      jHiRHSLoc(ISYM,ICASE) = JHI
      nRHSLocSz(ISYM,ICASE) = NAS*(JHI-JLO+1)
      NALL = NALL+nRHSLocSz(ISYM,ICASE)
    end do
  end do

  NGRP = 0
  NSYM1 = 0
  do IGRP=1,NRHSGRP
    NBLK = 0
    do IPAIR=1,2
      if (IRHSGRP(IPAIR,IGRP) == 0) cycle
      NBLK = NBLK+sum(nRHSLocSz(1:NSYM,IRHSGRP(IPAIR,IGRP)))
    end do
    NGRP = max(NGRP,NBLK)
    do ISYM=1,NSYM
      NBLK = 0
      do IPAIR=1,2
        if (IRHSGRP(IPAIR,IGRP) == 0) cycle
        NBLK = NBLK+nRHSLocSz(ISYM,IRHSGRP(IPAIR,IGRP))
      end do
      NSYM1 = max(NSYM1,NBLK)
    end do
  end do

end subroutine RHSLOC_SIZES


!-----------------------------------------------------------------------

subroutine RHSLOC_ALLOCATE(ITIER,NBUF)

  use caspt2_module, only: NSYM
  use stdalloc, only: mma_allocate

  integer(kind=iwp), intent(in) :: ITIER, NBUF

  integer(kind=iwp) :: ICASE, IOFF, ISYM

  ! Allocate the buffer depending on the tier (ITIER).
  ! At tier 0 the whole layout is fixed here and RHSLOC_LOAD does nothing
  ! The lower tiers lay out one group at a time

  iRHSLocTier = ITIER
  iRHSLocGrp = 0
  iRHSLocSym = 0
  ISYMTGT = 0

  call mma_allocate(RHSLoc,max(NBUF,1),Label='RHSLoc')
  RHSLoc(:) = Zero

  if (iRHSLocTier == 0) then
    IOFF = 1
    do ICASE=1,13
      do ISYM=1,NSYM
        if (nRHSLocSz(ISYM,ICASE) == 0) cycle
        iOffRHSLoc(ISYM,ICASE) = IOFF
        IOFF = IOFF+nRHSLocSz(ISYM,ICASE)
      end do
    end do
  end if

end subroutine RHSLOC_ALLOCATE

!-----------------------------------------------------------------------

subroutine RHSLOC_LOAD(IVEC,IGRP,ISYMT)

  use caspt2_module, only: NSYM

  integer(kind=iwp), intent(in) :: IVEC, IGRP, ISYMT

  integer(kind=iwp) :: ICASE, IOFF, IPAIR, ISYM, NBLK

  ! Make the blocks of case group IGRP in memory, writing back the previous ones

  ISYMTGT = 0
  if (iRHSLocTier == 2) ISYMTGT = ISYMT

  if (iRHSLocTier == 0) return
  if ((IGRP == iRHSLocGrp) .and. (ISYMT == iRHSLocSym)) return

  call RHSLOC_SAVE(IVEC)

  iOffRHSLoc(:,:) = 0
  IOFF = 1
  do IPAIR=1,2
    ICASE = IRHSGRP(IPAIR,IGRP)
    if (ICASE == 0) cycle
    do ISYM=1,NSYM
      if ((ISYMTGT /= 0) .and. (ISYM /= ISYMTGT)) cycle
      NBLK = nRHSLocSz(ISYM,ICASE)
      if (NBLK == 0) cycle
      iOffRHSLoc(ISYM,ICASE) = IOFF
      if (RHSLocOnDisk(ISYM,ICASE)) then
        call RHSLOC_IO(IVEC,ISYM,ICASE,2,IOFF,NBLK)
      else
        RHSLoc(IOFF:IOFF+NBLK-1) = Zero
      end if
      IOFF = IOFF+NBLK
    end do
  end do
  iRHSLocGrp = IGRP
  iRHSLocSym = ISYMT

end subroutine RHSLOC_LOAD

!-----------------------------------------------------------------------

subroutine RHSLOC_SAVE(IVEC)

  use caspt2_module, only: NSYM

  integer(kind=iwp), intent(in) :: IVEC

  integer(kind=iwp) :: ICASE, ISYM, NBLK

  ! Write the in-core blocks back to their slots on LURHS
  ! At tier 0 this is reached once, at the end of RHSALL2_STRIPED

  do ICASE=1,13
    do ISYM=1,NSYM
      if (iOffRHSLoc(ISYM,ICASE) == 0) cycle
      NBLK = nRHSLocSz(ISYM,ICASE)
      if (NBLK == 0) cycle
      call RHSLOC_IO(IVEC,ISYM,ICASE,1,iOffRHSLoc(ISYM,ICASE),NBLK)
      RHSLocOnDisk(ISYM,ICASE) = .true.
    end do
  end do

end subroutine RHSLOC_SAVE

!-----------------------------------------------------------------------

subroutine RHSLOC_IO(IVEC,ISYM,ICASE,IOPT,IOFF,NBLK)

  use caspt2_global, only: LURHS
  use caspt2_module, only: IOFFRHS

  integer(kind=iwp), intent(in) :: IVEC, ISYM, ICASE, IOPT, IOFF, NBLK

  integer(kind=iwp) :: IDISK

  ! One block of RHSLoc to (IOPT = 1) or from (IOPT = 2) its slot on LURHS.

  IDISK = IOFFRHS(ISYM,ICASE)
  call DDAFILE(LURHS(IVEC),IOPT,RHSLoc(IOFF),NBLK,IDISK)

end subroutine RHSLOC_IO

!-----------------------------------------------------------------------

subroutine RHSLOC_FINALIZE(IVEC)

  use caspt2_module, only: NSYM

  integer(kind=iwp), intent(in) :: IVEC

  integer(kind=iwp) :: IGRP, ISYMT, NSYMT

  ! Write out what is still in memory
  ! All blocks at tier 0, the last group below that
  ! iRHSLocGrp = 0 means no group was ever loaded, i.e. no Cholesky batches

  if ((iRHSLocTier /= 0) .and. (iRHSLocGrp == 0)) then
    NSYMT = 1
    if (iRHSLocTier == 2) NSYMT = NSYM
    do IGRP=1,NRHSGRP
      do ISYMT=1,NSYMT
        call RHSLOC_LOAD(IVEC,IGRP,ISYMT)
      end do
    end do
  end if
  call RHSLOC_SAVE(IVEC)

end subroutine RHSLOC_FINALIZE

!-----------------------------------------------------------------------

subroutine RHSLOC_FREE()

  use stdalloc, only: mma_deallocate

  call mma_deallocate(RHSLoc)

end subroutine RHSLOC_FREE

!-----------------------------------------------------------------------

subroutine RHSLOC_BOUNDS(ISYM,ICASE,ACTIVE,MOFF,JLO,JHI)

  integer(kind=iwp), intent(in) :: ISYM, ICASE
  logical(kind=iwp), intent(in) :: ACTIVE
  integer(kind=iwp), intent(out) :: MOFF, JLO, JHI

  ! Where the local stripe of one block sits in RHSLoc, for the kernels to index with
  ! An absent stripe is reported as JLO:JHI = 1:0 with MOFF = 1
  ! RHSLoc(MOFF) can then be passed as a zero-column array

  MOFF = 1
  JLO = 1
  JHI = 0
  if (.not. ACTIVE) return
  if (jHiRHSLoc(ISYM,ICASE) < jLoRHSLoc(ISYM,ICASE)) return

  JLO = jLoRHSLoc(ISYM,ICASE)
  JHI = jHiRHSLoc(ISYM,ICASE)
  MOFF = iOffRHSLoc(ISYM,ICASE)
  if (MOFF == 0) then
    ! the caller asked for a block that RHSLOC_LOAD did not bring in
    write(u6,'(1X,A,2I4)') 'RHSLOC_BOUNDS: block not in memory, ISYM, ICASE =',ISYM,ICASE
    call AbEnd()
  end if

end subroutine RHSLOC_BOUNDS

!-----------------------------------------------------------------------

function SKIP_SYM(ISYM)

  logical(kind=iwp) :: SKIP_SYM
  integer(kind=iwp), intent(in) :: ISYM

  ! True if the block of symmetry ISYM is not in memory and the kernel should return.
  ! Only tier 2 sets ISYMTGT; below that no symmetry is ever skipped.

  SKIP_SYM = (ISYMTGT /= 0) .and. (ISYM /= ISYMTGT)

end function SKIP_SYM

!-----------------------------------------------------------------------
!
! The rest is for Cases A ~ H
!
! Naming used throughout the kernels below.
!
! Passed in:
!   W...            the local stripe of an RHS block
!   JLO/JHI, MOFF   the column range of that stripe and where it starts in RHSLoc, from RHSLOC_BOUNDS
!                   MOFF counts elements, not columns; an absent stripe is JLO:JHI = 1:0
!   SCR, NSCR       the scratch the integral block of the current batch is built in, and its size
!
! Indices and sizes:
!   N<x>MX          the number of values of index x one batch may hold at most
!   N<x>SZ          the number of values x holds in the current batch
!   NACMX           the largest NAMX*NCMX the scratch allows, the kernel splits it (F, G)
!   I<x>STA/I<x>END first and last index of x in the current batch, respectively
!   I<x>LO/I<x>HI   first and last index of x this process needs at all
!   IQLO/IQHI       first and last index of q, the first member of an unordered superindex pair
!   IQ, IQSTA/IQEND that index and its batch, in the off-diagonal kernels of F and G
!   ...P/...M       the plus (t>=u, i>=j, a>=b) and minus (t>u, ...) combination of a case, wherever both have to be kept apart
!   ...A            an array holding one such quantity per a of an a-block
!
! Positions within a block:
!   JBASP/JBASM     index of the column of the pair (x,1)
!   IRBASP/IRBASM   index of the first row a run of DAXPYs writes
!   IROFF           row offset of such a run, i.e. its i-th row is IROFF+i
!   LDY             the leading dimension of the integral block in SCR
!   IY1/IY2         the two integral blocks in SCR when a kernel needs both
!   I<x>L           position of x within the current batch, i.e. x = I<x>STA+I<x>L-1
!
! A case with plus/minus combinations is split into two kernels,
! selected by whether the bra and ket symmetry labels agree:
!   ..._D  the diagonal block, one integral block gives both terms
!   ..._O  the off-diagonal one, one term per iteration and a uniform weight
! A superindex pair is written (q,r) with q >= r; q is the first argument of KAGEB/KAGTB and KIGEJ/KIGTJ.
! The bra index is q when its symmetry is the higher one.
!
!-----------------------------------------------------------------------

subroutine ADDRHSA_STRIPED(JSYM,ISYJ,ISYX,NT,NJ,NV,NX, &
                           SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)

  use SUPERINDEX, only: KTUV
  use caspt2_module, only: NAES, NINDEP, NISH, NTUV, NTUVES
  use Symmetry_Info, only: Mul
  use SC_NEVPT2, only: Do_SC

  integer(kind=iwp), intent(in) :: JSYM, ISYJ, ISYX, NT, NJ, NV, NX, NSCR, NCHO
  real(kind=wp), intent(out), target :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NT*NJ,NCHO), Cho_Ket(NV*NX,NCHO)

  integer(kind=iwp) :: IJ, IJEND, IJSTA, IROFF, ISYM, ISYT, ISYV, JHI, JLO, LDY, MOFF, NAS, NIS, NJMX, NJSZ

  real(kind=wp), pointer, contiguous :: WBLK(:,:)   ! the local stripe of the block, as (row,column)
  real(kind=wp), pointer, contiguous :: SCR3(:,:,:) ! the DGEMM output, viewed as (t,j,vx)
  real(kind=wp), pointer, contiguous :: WCOL(:,:)   ! one WBLK column's (t,vx) sub-block

  ! Case A: W(tvx,j) = (tj,vx)
  !   SCR((t,j),(v,x)) = sum_P Cho_Bra(t,j)^P Cho_Ket(v,x)^P
  ! j in chunks of NJMX columns, Cho_Bra is (t,j,P)

  ISYT = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYX)
  ISYM = ISYJ
  if (SKIP_SYM(ISYM)) return
  if ((NINDEP(ISYM,1) == 0) .and. (.not. Do_SC)) return
  NAS = NTUV(ISYM)
  NIS = NISH(ISYM)
  if (NAS*NIS == 0) return
  NJMX = NSCR/(NT*NV*NX)
  if (NJMX < 1) then
    write(u6,*) 'Not enough memory in ADDRHSA_STRIPED, I give up'
    call Abend()
  end if
  IROFF = KTUV(1+NAES(ISYT),1+NAES(ISYV),1+NAES(ISYX))-NTUVES(ISYM)-1

  call RHSLOC_BOUNDS(ISYM,1,NAS*NIS>0,MOFF,JLO,JHI)
  if (JHI >= JLO) then
    WBLK(1:NAS,JLO:JHI) => RHSLoc(MOFF:MOFF+nRHSLocSz(ISYM,1)-1)
    NJMX = min(NJMX,NJ)
    do IJSTA=JLO,JHI,NJMX
      IJEND = min(IJSTA+NJMX-1,JHI)
      NJSZ = IJEND-IJSTA+1
      LDY = NT*NJSZ
      call DGEMM_('N','T',LDY,NV*NX,NCHO,One,Cho_Bra(1+NT*(IJSTA-1),1),NT*NJ,Cho_Ket,NV*NX,Zero,SCR,LDY)
      SCR3(1:NT,1:NJSZ,1:NV*NX) => SCR(1:LDY*NV*NX)
      do IJ=IJSTA,IJEND
        WCOL(1:NT,1:NV*NX) => WBLK(IROFF+1:IROFF+NT*NV*NX,IJ)
        WCOL(1:NT,1:NV*NX) = WCOL(1:NT,1:NV*NX)+SCR3(1:NT,IJ-IJSTA+1,1:NV*NX)
      end do
    end do
    nullify(WBLK,SCR3,WCOL)
  end if

end subroutine ADDRHSA_STRIPED

!-----------------------------------------------------------------------

subroutine ADDRHSB_STRIPED(JSYM,ISYJ,ISYL,NT,NJ,NV,NL, &
                           SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)

  use caspt2_module, only: NIGEJ, NIGTJ, NINDEP, NTGEU, NTGTU
  use Symmetry_Info, only: Mul
  use SC_NEVPT2, only: Do_SC

  integer(kind=iwp), intent(in) :: JSYM, ISYJ, ISYL, NT, NJ, NV, NL, NSCR, NCHO
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NT*NJ*NCHO), Cho_Ket(NV*NL*NCHO)

  integer(kind=iwp) :: ISYM, ISYT, ISYV, JHIM, JHIP, JLOM, JLOP, MOFFM, MOFFP, NASM, NASP, NISM, NISP, NLMX

  ! Case B: both combinations of (tj,vl)
  !   rows    the active pair (t,v)
  !   columns the inactive pair (j,l)

  ISYT = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYL)
  if (ISYT < ISYV) return
  ISYM = Mul(ISYJ,ISYL)
  if (SKIP_SYM(ISYM)) return

  NASP = NTGEU(ISYM)
  NISP = NIGEJ(ISYM)
  NASM = NTGTU(ISYM)
  NISM = NIGTJ(ISYM)
  ! a combination with no independent parameters is not built, as in ADDRHSB
  if ((NINDEP(ISYM,2) == 0) .and. (.not. Do_SC)) NISP = 0
  if ((NINDEP(ISYM,3) == 0) .and. (.not. Do_SC)) NISM = 0
  if (NASP*NISP+NASM*NISM == 0) return

  NLMX = NSCR/(NT*NV)
  if (NLMX < 1) then
    write(u6,*) 'Not enough memory in ADDRHSB_STRIPED, I give up'
    call Abend()
  end if

  call RHSLOC_BOUNDS(ISYM,2,NASP*NISP>0,MOFFP,JLOP,JHIP)
  call RHSLOC_BOUNDS(ISYM,3,NASM*NISM>0,MOFFM,JLOM,JHIM)

  ! the bra and the ket are the same block of Cholesky vectors when their symmetry labels agree
  if (ISYJ == ISYL) then
    call ADDRHSB_STRIPED_D(RHSLoc(MOFFP),RHSLoc(MOFFM),NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYT,ISYM,NT,NJ, &
                           SCR,NSCR,NLMX,Cho_Bra,NCHO)
  else
    call ADDRHSB_STRIPED_O(RHSLoc(MOFFP),RHSLoc(MOFFM),NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYL,ISYT,ISYV,ISYM,NT,NJ,NV,NL, &
                           SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)
  end if

end subroutine ADDRHSB_STRIPED

!-----------------------------------------------------------------------

subroutine ADDRHSB_STRIPED_D(WBP,WBM,NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYT,ISYM,NT,NJ, &
                             SCR,NSCR,NLMX,Cho_Bra,NCHO)

  use SUPERINDEX, only: KIGEJ, KIGTJ, KTGEU, KTGTU
  use caspt2_module, only: NAES, NIES, NIGEJES, NIGTJES, NTGEUES, NTGTUES
  use Constants, only: Quart

  integer(kind=iwp), intent(in) :: NASP, NASM, JLOP, JHIP, JLOM, JHIM, ISYJ, ISYT, ISYM, NT, NJ, NSCR, NLMX, NCHO
  real(kind=wp), intent(inout) :: WBP(NASP,JLOP:JHIP), WBM(NASM,JLOM:JHIM)
  real(kind=wp), intent(out), target :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NT*NJ,NCHO)

  integer(kind=iwp) :: ICOL, IJ, IJABS, IL, ILEND, ILHI, ILHIM, ILHIP, ILL, ILLO, ILLOM, ILLOP, ILSTA, IRBASM, IRBASP, IV, &
                       IVABS, JBASM, JBASP, LDY, NLSZ

  real(kind=wp), pointer, contiguous :: Y(:,:,:) ! the integral block of this batch, as (v,l,t)

  ! Case B, diagonal block.
  !   WP(t>=v,j>=l) = ((tj,vl)+(tl,vj))/SQRT((1+Kron(jl))*...), the Half/Quart weights of
  !                   the replicated code
  !   WM(t>v,j>l)   = ((tj,vl)-(tl,vj))/2
  ! same shape as case H, one integral block per (j, l-range) gives both terms
  ! KTGEU/KTGTU is (t,v) with t fastest, rows t = v..NT are contiguous

  do IJ=1,NJ
    IJABS = IJ+NIES(ISYJ)
    ! the columns of l = 1..j (WBP) or 1..j-1 (WBM) are contiguous from JBASP/JBASM
    JBASP = KIGEJ(IJABS,1+NIES(ISYJ))-NIGEJES(ISYM)
    ILLOP = max(1,JLOP-JBASP+1)
    ILHIP = min(IJ,JHIP-JBASP+1)
    JBASM = 0
    ILLOM = 1
    ILHIM = 0
    if ((IJ >= 2) .and. (NASM > 0)) then
      JBASM = KIGTJ(IJABS,1+NIES(ISYJ))-NIGTJES(ISYM)
      ILLOM = max(1,JLOM-JBASM+1)
      ILHIM = min(IJ-1,JHIM-JBASM+1)
    end if

    ! ILLO:ILHI: union of the l ranges the two combinations need
    if ((ILHIP < ILLOP) .and. (ILHIM < ILLOM)) cycle
    if (ILHIP < ILLOP) then
      ILLO = ILLOM
      ILHI = ILHIM
    else if (ILHIM < ILLOM) then
      ILLO = ILLOP
      ILHI = ILHIP
    else
      ILLO = min(ILLOP,ILLOM)
      ILHI = max(ILHIP,ILHIM)
    end if

    do ILSTA=ILLO,ILHI,NLMX
      ILEND = min(ILSTA+NLMX-1,ILHI)
      NLSZ = ILEND-ILSTA+1
      LDY = NT*NLSZ
      ! SCR((v,l),t) = sum_P Cho_Bra(v,l)^P Cho_Bra(t,j)^P = (tj,vl)
      call DGEMM_('N','T',LDY,NT,NCHO,One,Cho_Bra(1+NT*(ILSTA-1),1),NT*NJ,Cho_Bra(1+NT*(IJ-1),1),NT*NJ,Zero,SCR,LDY)
      Y(1:NT,1:NLSZ,1:NT) => SCR(1:LDY*NT)

      do IL=ILSTA,ILEND
        ILL = IL-ILSTA+1

        ! term1 = Y(v,l,t) = (tj,vl), term2 = Y(t,l,v) = (tl,vj); at t == v the two are the same element
        ! plus combination
        if ((IL >= ILLOP) .and. (IL <= ILHIP)) then
          ICOL = JBASP+IL-1
          do IV=1,NT
            IVABS = IV+NAES(ISYT)
            ! rows KTGEU(t,v) for t = v..NT are contiguous
            IRBASP = KTGEU(IVABS,IVABS)-NTGEUES(ISYM)
            if (IL == IJ) then
              WBP(IRBASP,ICOL) = WBP(IRBASP,ICOL)+SQ2*Quart*Y(IV,ILL,IV)
              if (IV < NT) WBP(IRBASP+1:IRBASP+NT-IV,ICOL) = WBP(IRBASP+1:IRBASP+NT-IV,ICOL) &
                                                             +SQ2*Half*Y(IV,ILL,IV+1:NT)
            else
              WBP(IRBASP,ICOL) = WBP(IRBASP,ICOL)+Half*Y(IV,ILL,IV)
              if (IV < NT) WBP(IRBASP+1:IRBASP+NT-IV,ICOL) = WBP(IRBASP+1:IRBASP+NT-IV,ICOL) &
                                                             +Half*(Y(IV,ILL,IV+1:NT)+Y(IV+1:NT,ILL,IV))
            end if
          end do
        end if

        ! minus combination
        if ((IL >= ILLOM) .and. (IL <= ILHIM)) then
          ICOL = JBASM+IL-1
          do IV=1,NT-1
            IVABS = IV+NAES(ISYT)
            ! rows KTGTU(t,v) for t = v+1..NT are contiguous
            IRBASM = KTGTU(IVABS+1,IVABS)-NTGTUES(ISYM)
            WBM(IRBASM:IRBASM+NT-IV-1,ICOL) = WBM(IRBASM:IRBASM+NT-IV-1,ICOL) &
                                              +Half*(Y(IV,ILL,IV+1:NT)-Y(IV+1:NT,ILL,IV))
          end do
        end if

      end do
    end do
  end do
  nullify(Y)

end subroutine ADDRHSB_STRIPED_D

!-----------------------------------------------------------------------

subroutine ADDRHSB_STRIPED_O(WBP,WBM,NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYL,ISYT,ISYV,ISYM,NT,NJ,NV,NL, &
                             SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)

  use SUPERINDEX, only: KIGEJ, KIGTJ, KTGEU, KTGTU
  use caspt2_module, only: NAES, NIES, NIGEJES, NIGTJES, NTGEUES, NTGTUES

  integer(kind=iwp), intent(in) :: NASP, NASM, JLOP, JHIP, JLOM, JHIM, ISYJ, ISYL, ISYT, ISYV, ISYM, NT, NJ, NV, NL, NSCR, NCHO
  real(kind=wp), intent(inout) :: WBP(NASP,JLOP:JHIP), WBM(NASM,JLOM:JHIM)
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NT*NJ,NCHO), Cho_Ket(NV*NL,NCHO)

  integer(kind=iwp) :: ICOLM, ICOLP, IJ, IJABS, IL, ILABS, IRBASM, IRBASP, IV, IVABS
  real(kind=wp) :: SGN

  ! Case B, off-diagonal block (ISYM /= 1, hence ISYT > ISYV and ISYJ /= ISYL)
  ! one term, full rectangles, uniform Half weight
  ! KTGEU(t,v) is contiguous in t for a fixed v, each (j,l) is one column

  ! j is q when ISYJ > ISYL; that affects only the columns, no range of q is selected as in F and G
  SGN = One
  if (ISYJ <= ISYL) SGN = -One

  do IJ=1,NJ
    IJABS = IJ+NIES(ISYJ)
    do IL=1,NL
      ILABS = IL+NIES(ISYL)
      if (ISYJ > ISYL) then
        ICOLP = KIGEJ(IJABS,ILABS)-NIGEJES(ISYM)
        ICOLM = KIGTJ(IJABS,ILABS)-NIGTJES(ISYM)
      else
        ICOLP = KIGEJ(ILABS,IJABS)-NIGEJES(ISYM)
        ICOLM = KIGTJ(ILABS,IJABS)-NIGTJES(ISYM)
      end if
      if (((ICOLP < JLOP) .or. (ICOLP > JHIP)) .and. ((NASM == 0) .or. (ICOLM < JLOM) .or. (ICOLM > JHIM))) cycle
      ! SCR(v,t) = sum_P Cho_Ket(v,l)^P Cho_Bra(t,j)^P = (tj,vl)
      call DGEMM_('N','T',NV,NT,NCHO,One,Cho_Ket(1+NV*(IL-1),1),NV*NL,Cho_Bra(1+NT*(IJ-1),1),NT*NJ,Zero,SCR,NV)
      do IV=1,NV
        IVABS = IV+NAES(ISYV)
        ! plus combination
        if ((ICOLP >= JLOP) .and. (ICOLP <= JHIP)) then
          IRBASP = KTGEU(1+NAES(ISYT),IVABS)-NTGEUES(ISYM)
          call DAXPY_(NT,Half,SCR(IV),NV,WBP(IRBASP,ICOLP),1)
        end if
        ! minus combination
        if ((NASM > 0) .and. (ICOLM >= JLOM) .and. (ICOLM <= JHIM)) then
          IRBASM = KTGTU(1+NAES(ISYT),IVABS)-NTGTUES(ISYM)
          call DAXPY_(NT,SGN*Half,SCR(IV),NV,WBM(IRBASM,ICOLM),1)
        end if
      end do
    end do
  end do

end subroutine ADDRHSB_STRIPED_O

!-----------------------------------------------------------------------

subroutine ADDRHSC_STRIPED(JSYM,ISYU,ISYX,NA,NU,NV,NX, &
                           SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)

  use SUPERINDEX, only: KTUV
  use stdalloc, only: mma_allocate, mma_deallocate
  use caspt2_module, only: NAES, NINDEP, NSSH, NTUV, NTUVES
  use Symmetry_Info, only: Mul
  use SC_NEVPT2, only: Do_SC

  integer(kind=iwp), intent(in) :: JSYM, ISYU, ISYX, NA, NU, NV, NX, NSCR, NCHO
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NV*NX,NCHO)

  integer(kind=iwp) :: IA, IP, IR, IROFF, ISYM, ISYV, IU, JHI, JLO, MOFF, NAS, NIS

  real(kind=wp), allocatable :: CHOBA(:)
  real(kind=wp), pointer, contiguous :: WBLK(:,:) ! the local stripe of the block, as (row,column)

  ! Case C: W(uvx,a) = (au,vx)
  ! same shape as case A, but the bra is (a,u,P), packed and contracted per a

  ISYV = Mul(JSYM,ISYX)
  ISYM = Mul(JSYM,ISYU)
  if (SKIP_SYM(ISYM)) return
  if ((NINDEP(ISYM,4) == 0) .and. (.not. Do_SC)) return
  NAS = NTUV(ISYM)
  NIS = NSSH(ISYM)
  if (NAS*NIS == 0) return
  if (NSCR < NU*NV*NX) then
    write(u6,*) 'Not enough memory in ADDRHSC_STRIPED, I give up'
    call Abend()
  end if
  IROFF = KTUV(1+NAES(ISYU),1+NAES(ISYV),1+NAES(ISYX))-NTUVES(ISYM)-1

  call RHSLOC_BOUNDS(ISYM,4,NAS*NIS>0,MOFF,JLO,JHI)
  if (JHI >= JLO) then
    WBLK(1:NAS,JLO:JHI) => RHSLoc(MOFF:MOFF+nRHSLocSz(ISYM,4)-1)
    call mma_allocate(CHOBA,NU*NCHO,Label='CHOBA')
    do IA=JLO,JHI
      do IP=1,NCHO
        do IU=1,NU
          CHOBA(IU+NU*(IP-1)) = Cho_Bra(IA,IU,IP) ! the (u,P) slice of this a
        end do
      end do
      ! SCR(u,(v,x)) = sum_P CHOBA(u)^P Cho_Ket(v,x)^P = (au,vx)
      call DGEMM_('N','T',NU,NV*NX,NCHO,One,CHOBA,NU,Cho_Ket,NV*NX,Zero,SCR,NU)
      do IR=1,NV*NX
        call DAXPY_(NU,One,SCR(1+NU*(IR-1)),1,WBLK(IROFF+NU*(IR-1)+1,IA),1)
      end do
    end do
    call mma_deallocate(CHOBA)
    nullify(WBLK)
  end if

end subroutine ADDRHSC_STRIPED

!-----------------------------------------------------------------------

subroutine ADDRHSD1_STRIPED(JSYM,ISYJ,ISYX,NA,NJ,NV,NX, &
                            SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)

  use SUPERINDEX, only: KTU
  use stdalloc, only: mma_allocate, mma_deallocate
  use caspt2_module, only: NAES, NINDEP, NISH, NISUP, NSSH, NSYM, NTU, NTUES
  use Symmetry_Info, only: Mul
  use SC_NEVPT2, only: Do_SC

  integer(kind=iwp), intent(in) :: JSYM, ISYJ, ISYX, NA, NJ, NV, NX, NSCR, NCHO
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA*NJ,NCHO), Cho_Ket(NV*NX,NCHO)

  integer(kind=iwp) :: IA, ICOL, IJ, IJHI, IJLO, IOFF, IOFFD, IP, IROFF, ISA, ISI, ISYA, ISYM, ISYV, JHI, JLO, MOFF, &
                       NAS, NIS, NJSZ

  real(kind=wp), allocatable :: CHOBA(:)
  real(kind=wp), pointer, contiguous :: WBLK(:,:) ! the local stripe of the block, as (row,column)

  ! Case D1: W(vx, IOFFD+j+NJ*(a-1)) = (aj,vx)
  ! no pair index, no symmetrization
  ! KTU(v,x) is contiguous with v fastest, one column = one run of NV*NX rows
  ! the bra is (a,j,P), packed per a

  ISYA = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYX)
  ISYM = JSYM
  if (SKIP_SYM(ISYM)) return
  if ((NINDEP(ISYM,5) == 0) .and. (.not. Do_SC)) return
  NAS = 2*NTU(ISYM)
  NIS = NISUP(ISYM,5)
  if (NAS*NIS == 0) return
  if (NSCR < NV*NX*NJ) then
    write(u6,*) 'Not enough memory in ADDRHSD1_STRIPED, I give up'
    call Abend()
  end if

  IOFFD = 0
  do ISA=1,NSYM
    if (ISA == ISYA) exit
    ISI = Mul(ISA,ISYM)
    IOFFD = IOFFD+NSSH(ISA)*NISH(ISI)
  end do
  IROFF = KTU(1+NAES(ISYV),1+NAES(ISYX))-NTUES(ISYM)-1

  call RHSLOC_BOUNDS(ISYM,5,NAS*NIS>0,MOFF,JLO,JHI)
  if (JHI >= JLO) then
    WBLK(1:NAS,JLO:JHI) => RHSLoc(MOFF:MOFF+nRHSLocSz(ISYM,5)-1)
    call mma_allocate(CHOBA,NJ*NCHO,Label='CHOBA')
    do IA=1,NA
      ! IJLO:IJHI: the j range of this a whose columns fall inside the local stripe JLO:JHI
      IJLO = max(1,JLO-IOFFD-NJ*(IA-1))
      IJHI = min(NJ,JHI-IOFFD-NJ*(IA-1))
      if (IJLO > IJHI) cycle
      NJSZ = IJHI-IJLO+1
      do IP=1,NCHO
        do IJ=IJLO,IJHI
          CHOBA(IJ-IJLO+1+NJSZ*(IP-1)) = Cho_Bra(IA+NA*(IJ-1),IP) ! the (j,P) slice of this a, local j only
        end do
      end do
      ! SCR((v,x),j)
      call DGEMM_('N','T',NV*NX,NJSZ,NCHO,One,Cho_Ket,NV*NX,CHOBA,NJSZ,Zero,SCR,NV*NX)
      do IJ=IJLO,IJHI
        ICOL = IOFFD+NJ*(IA-1)+IJ
        IOFF = NV*NX*(IJ-IJLO)
        WBLK(IROFF+1:IROFF+NV*NX,ICOL) = WBLK(IROFF+1:IROFF+NV*NX,ICOL)+SCR(IOFF+1:IOFF+NV*NX)
      end do
    end do
    call mma_deallocate(CHOBA)
    nullify(WBLK)
  end if

end subroutine ADDRHSD1_STRIPED

!-----------------------------------------------------------------------

subroutine ADDRHSD2_STRIPED(JSYM,ISYU,ISYL,NA,NU,NV,NL, &
                            SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)

  use SUPERINDEX, only: KTU
  use stdalloc, only: mma_allocate, mma_deallocate
  use caspt2_module, only: NAES, NINDEP, NISH, NISUP, NSSH, NSYM, NTU, NTUES
  use Symmetry_Info, only: Mul
  use SC_NEVPT2, only: Do_SC

  integer(kind=iwp), intent(in) :: JSYM, ISYU, ISYL, NA, NU, NV, NL, NSCR, NCHO
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NV*NL,NCHO)

  integer(kind=iwp) :: IA, ICOL, IL, ILHI, ILLO, IOFF, IOFFD, IP, IROFF, IROFFU, ISA, ISI, ISYA, ISYM, ISYV, IU, JHI, JLO, &
                       MOFF, NAS, NAS1, NIS, NLSZ

  real(kind=wp), allocatable :: CHOBA(:)
  real(kind=wp), pointer, contiguous :: WBLK(:,:) ! the local stripe of the block, as (row,column)

  ! Case D2: W(NAS1+vu, IOFFD+l+NL*(a-1)) = (au,vl)
  ! the fast column index l comes from the ket, one contraction per a over the owned l range
  ! SCR((v,l),u), one column = NU runs of NV rows, the u-th starting at IROFF+NV*(u-1)

  ISYA = Mul(JSYM,ISYU)
  ISYV = Mul(JSYM,ISYL)
  ISYM = Mul(ISYU,ISYV)
  if (SKIP_SYM(ISYM)) return
  if ((NINDEP(ISYM,5) == 0) .and. (.not. Do_SC)) return
  NAS1 = NTU(ISYM)
  NAS = 2*NAS1
  NIS = NISUP(ISYM,5)
  if (NAS*NIS == 0) return
  if (NSCR < NV*NL*NU) then
    write(u6,*) 'Not enough memory in ADDRHSD2_STRIPED, I give up'
    call Abend()
  end if

  IOFFD = 0
  do ISA=1,NSYM
    if (ISA == ISYA) exit
    ISI = Mul(ISA,ISYM)
    IOFFD = IOFFD+NSSH(ISA)*NISH(ISI)
  end do
  IROFF = NAS1+KTU(1+NAES(ISYV),1+NAES(ISYU))-NTUES(ISYM)-1

  call RHSLOC_BOUNDS(ISYM,5,NAS*NIS>0,MOFF,JLO,JHI)
  if (JHI >= JLO) then
    WBLK(1:NAS,JLO:JHI) => RHSLoc(MOFF:MOFF+nRHSLocSz(ISYM,5)-1)
    call mma_allocate(CHOBA,NU*NCHO,Label='CHOBA')
    do IA=1,NA
      ! ILLO:ILHI: the l range of this a whose columns fall inside the local stripe JLO:JHI
      ILLO = max(1,JLO-IOFFD-NL*(IA-1))
      ILHI = min(NL,JHI-IOFFD-NL*(IA-1))
      if (ILLO > ILHI) cycle
      NLSZ = ILHI-ILLO+1
      do IP=1,NCHO
        do IU=1,NU
          CHOBA(IU+NU*(IP-1)) = Cho_Bra(IA,IU,IP) ! the (u,P) slice of this a
        end do
      end do
      ! SCR((v,l),u)
      call DGEMM_('N','T',NV*NLSZ,NU,NCHO,One,Cho_Ket(1+NV*(ILLO-1),1),NV*NL,CHOBA,NU,Zero,SCR,NV*NLSZ)
      do IL=ILLO,ILHI
        ICOL = IOFFD+NL*(IA-1)+IL
        do IU=1,NU
          IROFFU = IROFF+NV*(IU-1)
          IOFF = NV*(IL-ILLO)+NV*NLSZ*(IU-1)
          WBLK(IROFFU+1:IROFFU+NV,ICOL) = WBLK(IROFFU+1:IROFFU+NV,ICOL)+SCR(IOFF+1:IOFF+NV)
        end do
      end do
    end do
    call mma_deallocate(CHOBA)
    nullify(WBLK)
  end if

end subroutine ADDRHSD2_STRIPED

!-----------------------------------------------------------------------

subroutine ADDRHSE_STRIPED(JSYM,ISYJ,ISYL,NA,NJ,NV,NL, &
                           SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)

  use caspt2_module, only: NISUP
  use general_data, only: NASH
  use Symmetry_Info, only: Mul

  integer(kind=iwp), intent(in) :: JSYM, ISYJ, ISYL, NA, NJ, NV, NL, NSCR, NCHO
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA*NJ,NCHO), Cho_Ket(NV*NL,NCHO)

  integer(kind=iwp) :: ISYA, ISYJL, ISYM, JHIM, JHIP, JLOM, JLOP, MOFFM, MOFFP, NAMX, NAS, NISM, NISP

  ! Case E: both combinations of (aj,vl)
  !   rows    the active index alone
  !   columns a secondary index within an inactive pair

  ISYA = Mul(JSYM,ISYJ)
  ISYM = Mul(JSYM,ISYL)
  if (SKIP_SYM(ISYM)) return
  ISYJL = Mul(ISYJ,ISYL)
  NAS = NASH(ISYM)
  NISP = NISUP(ISYM,6)
  NISM = NISUP(ISYM,7)
  if (NAS*NISP == 0) return

  NAMX = NSCR/max(2*NV,NV*NL)
  if (NAMX < 1) then
    write(u6,*) 'Not enough memory in ADDRHSE_STRIPED, I give up'
    call Abend()
  end if
  NAMX = min(NAMX,NA)

  call RHSLOC_BOUNDS(ISYM,6,NAS*NISP>0,MOFFP,JLOP,JHIP)
  call RHSLOC_BOUNDS(ISYM,7,NAS*NISM>0,MOFFM,JLOM,JHIM)

  ! the bra and the ket are the same block of Cholesky vectors when their symmetry labels agree
  if (ISYJ == ISYL) then
    call ADDRHSE_STRIPED_D(RHSLoc(MOFFP),RHSLoc(MOFFM),NAS,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYA,ISYM,ISYJL,NA,NJ,NV, &
                           SCR,NSCR,NAMX,Cho_Bra,Cho_Ket,NCHO)
  else
    call ADDRHSE_STRIPED_O(RHSLoc(MOFFP),RHSLoc(MOFFM),NAS,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYL,ISYA,ISYM,ISYJL,NA,NJ,NV,NL, &
                           SCR,NSCR,NAMX,Cho_Bra,Cho_Ket,NCHO)
  end if

end subroutine ADDRHSE_STRIPED

!-----------------------------------------------------------------------

subroutine ADDRHSE_STRIPED_D(WEP,WEM,NAS,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYA,ISYM,ISYJL,NA,NJ,NV, &
                             SCR,NSCR,NAMX,Cho_Bra,Cho_Ket,NCHO)

  use SUPERINDEX, only: KIGEJ, KIGTJ
  use caspt2_module, only: NIES, NIGEJ, NIGEJES, NIGTJ, NIGTJES, NSSH, NSYM
  use Symmetry_Info, only: Mul

  integer(kind=iwp), intent(in) :: NAS, JLOP, JHIP, JLOM, JHIM, ISYJ, ISYA, ISYM, ISYJL, NA, NJ, NV, NSCR, NAMX, NCHO
  real(kind=wp), intent(inout) :: WEP(NAS,JLOP:JHIP), WEM(NAS,JLOM:JHIM)
  real(kind=wp), intent(out), target :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA*NJ,NCHO), Cho_Ket(NV*NJ,NCHO)

  integer(kind=iwp) :: IA, IAEND, IAHI, IAL, IALO, IASTA, ICOL, IJ, IJABS, IL, ILHI, ILHIM, ILHIP, ILLO, ILLOM, ILLOP, IOFFM, &
                       IOFFP, ISA, ISIJ, IY1, IY2, JBASM, JBASP, JCOLHI, JCOLLO, NASZ

  real(kind=wp), pointer, contiguous :: Y1(:,:), Y2(:,:) ! the two integral blocks of this batch, as (v,a)

  ! Case E, diagonal block.
  !   WP(v,a,j>=l) = ((aj,vl)+(al,vj))/SQRT(2+2*Kron(jl))
  !   WM(v,a,j>l)  = ((aj,vl)-(al,vj))*SQRT(OneHalf)
  ! a column is one array assignment
  ! columns run a fastest within an inactive pair, the pairs of one j are consecutive
  ! the stripe is a range of (l,a) per j, both operands are contiguous as they stand

  IOFFP = 0
  IOFFM = 0
  do ISA=1,NSYM
    if (ISA == ISYA) exit
    ISIJ = Mul(ISA,ISYM)
    IOFFP = IOFFP+NSSH(ISA)*NIGEJ(ISIJ)
    IOFFM = IOFFM+NSSH(ISA)*NIGTJ(ISIJ)
  end do

  do IJ=1,NJ
    IJABS = IJ+NIES(ISYJ)
    ! the columns of this j are (a, pair JBASP+l-1), l = 1..j, a fastest
    ! JCOLLO:JCOLHI: that flat column range clipped to the local stripe
    ! ILLOP:ILHIP (ILLOM:ILHIM): the same range converted to l bounds
    JBASP = KIGEJ(IJABS,1+NIES(ISYJ))-NIGEJES(ISYJL)
    JCOLLO = max(IOFFP+NA*(JBASP-1)+1,JLOP)
    JCOLHI = min(IOFFP+NA*(JBASP+IJ-1),JHIP)
    if (JCOLLO <= JCOLHI) then
      ILLOP = (JCOLLO-IOFFP-NA*(JBASP-1)-1)/NA+1
      ILHIP = (JCOLHI-IOFFP-NA*(JBASP-1)-1)/NA+1
    else
      ILLOP = 1
      ILHIP = 0
    end if
    JBASM = 0
    ILLOM = 1
    ILHIM = 0
    if ((IJ >= 2) .and. (JHIM >= JLOM)) then
      JBASM = KIGTJ(IJABS,1+NIES(ISYJ))-NIGTJES(ISYJL)
      JCOLLO = max(IOFFM+NA*(JBASM-1)+1,JLOM)
      JCOLHI = min(IOFFM+NA*(JBASM+IJ-2),JHIM)
      if (JCOLLO <= JCOLHI) then
        ILLOM = (JCOLLO-IOFFM-NA*(JBASM-1)-1)/NA+1
        ILHIM = (JCOLHI-IOFFM-NA*(JBASM-1)-1)/NA+1
      end if
    end if

    ! ILLO:ILHI: union of the l ranges the two combinations need
    if ((ILHIP < ILLOP) .and. (ILHIM < ILLOM)) cycle
    if (ILHIP < ILLOP) then
      ILLO = ILLOM
      ILHI = ILHIM
    else if (ILHIM < ILLOM) then
      ILLO = ILLOP
      ILHI = ILHIP
    else
      ILLO = min(ILLOP,ILLOM)
      ILHI = max(ILHIP,ILHIM)
    end if

    do IL=ILLO,ILHI
      ! IALO:IAHI: the a range needed for this (j,l), as the union of the two stripes
      IALO = NA+1
      IAHI = 0
      if ((IL >= ILLOP) .and. (IL <= ILHIP)) then
        IALO = max(1,JLOP-IOFFP-NA*(JBASP+IL-2))
        IAHI = min(NA,JHIP-IOFFP-NA*(JBASP+IL-2))
      end if
      if ((IL >= ILLOM) .and. (IL <= ILHIM)) then
        IALO = min(IALO,max(1,JLOM-IOFFM-NA*(JBASM+IL-2)))
        IAHI = max(IAHI,min(NA,JHIM-IOFFM-NA*(JBASM+IL-2)))
      end if
      if (IAHI < IALO) cycle

      do IASTA=IALO,IAHI,NAMX
        IAEND = min(IASTA+NAMX-1,IAHI)
        NASZ = IAEND-IASTA+1
        IY1 = 1
        IY2 = 1+NV*NASZ
        ! Y1(v,a) = (aj,vl), Y2(v,a) = (al,vj)
        call DGEMM_('N','T',NV,NASZ,NCHO,One,Cho_Ket(1+NV*(IL-1),1),NV*NJ,Cho_Bra(IASTA+NA*(IJ-1),1),NA*NJ,Zero,SCR(IY1),NV)
        if (IL /= IJ) call DGEMM_('N','T',NV,NASZ,NCHO,One,Cho_Ket(1+NV*(IJ-1),1),NV*NJ,Cho_Bra(IASTA+NA*(IL-1),1),NA*NJ, &
          Zero,SCR(IY2),NV)
        Y1(1:NV,1:NASZ) => SCR(IY1:IY1+NV*NASZ-1)
        Y2(1:NV,1:NASZ) => SCR(IY2:IY2+NV*NASZ-1) ! it is used only when IL /= IJ, but a compiler complains...
        do IA=IASTA,IAEND
          IAL = IA-IASTA+1
          ! plus combination
          if ((IL >= ILLOP) .and. (IL <= ILHIP)) then
            ICOL = IOFFP+NA*(JBASP+IL-2)+IA
            if ((ICOL >= JLOP) .and. (ICOL <= JHIP)) then
              if (IL == IJ) then
                WEP(1:NV,ICOL) = WEP(1:NV,ICOL)+Y1(:,IAL)
              else
                WEP(1:NV,ICOL) = WEP(1:NV,ICOL)+SQH*(Y1(:,IAL)+Y2(:,IAL))
              end if
            end if
          end if
          ! minus combination
          if ((IL /= IJ) .and. (IL >= ILLOM) .and. (IL <= ILHIM)) then
            ICOL = IOFFM+NA*(JBASM+IL-2)+IA
            if ((ICOL >= JLOM) .and. (ICOL <= JHIM)) WEM(1:NV,ICOL) = WEM(1:NV,ICOL)+SQ32*(Y1(:,IAL)-Y2(:,IAL))
          end if
        end do
      end do
    end do
  end do
  nullify(Y1,Y2)

end subroutine ADDRHSE_STRIPED_D

!-----------------------------------------------------------------------

subroutine ADDRHSE_STRIPED_O(WEP,WEM,NAS,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYL,ISYA,ISYM,ISYJL,NA,NJ,NV,NL, &
                             SCR,NSCR,NAMX,Cho_Bra,Cho_Ket,NCHO)

  use SUPERINDEX, only: KIGEJ, KIGTJ
  use caspt2_module, only: NIES, NIGEJ, NIGEJES, NIGTJ, NIGTJES, NSSH, NSYM
  use Symmetry_Info, only: Mul

  integer(kind=iwp), intent(in) :: NAS, JLOP, JHIP, JLOM, JHIM, ISYJ, ISYL, ISYA, ISYM, ISYJL, NA, NJ, NV, NL, NSCR, NAMX, NCHO
  real(kind=wp), intent(inout) :: WEP(NAS,JLOP:JHIP), WEM(NAS,JLOM:JHIM)
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA*NJ,NCHO), Cho_Ket(NV*NL,NCHO)

  integer(kind=iwp) :: IA, IAEND, IAHI, IALO, IASTA, IJ, IJABS, IL, ILABS, IOFFM, IOFFP, ISA, ISIJ, JBASM, JBASP, JGEL, JGTL, &
                       LDY, NASZ
  real(kind=wp) :: SGN

  ! Case E, off-diagonal block (ISYJ /= ISYL): one term, uniform weight.
  ! all l of a given j in one contraction, the GEMM is NV*NL wide instead of NV
  ! a is batched on the n side

  LDY = NV*NL

  ! j is q when ISYJ > ISYL; that affects only the columns, no range of q is selected as in F and G
  SGN = One
  if (ISYJ <= ISYL) SGN = -One

  IOFFP = 0
  IOFFM = 0
  do ISA=1,NSYM
    if (ISA == ISYA) exit
    ISIJ = Mul(ISA,ISYM)
    IOFFP = IOFFP+NSSH(ISA)*NIGEJ(ISIJ)
    IOFFM = IOFFM+NSSH(ISA)*NIGTJ(ISIJ)
  end do

  do IJ=1,NJ
    IJABS = IJ+NIES(ISYJ)
    ! IALO:IAHI: union of the a ranges this j needs, over all its l
    IALO = NA+1
    IAHI = 0
    do IL=1,NL
      ILABS = IL+NIES(ISYL)
      if (ISYJ > ISYL) then
        JGEL = KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
        JGTL = KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
      else
        JGEL = KIGEJ(ILABS,IJABS)-NIGEJES(ISYJL)
        JGTL = KIGTJ(ILABS,IJABS)-NIGTJES(ISYJL)
      end if
      JBASP = IOFFP+NA*(JGEL-1)
      JBASM = IOFFM+NA*(JGTL-1)
      if (max(1,JLOP-JBASP) <= min(NA,JHIP-JBASP)) then
        IALO = min(IALO,max(1,JLOP-JBASP))
        IAHI = max(IAHI,min(NA,JHIP-JBASP))
      end if
      if (max(1,JLOM-JBASM) <= min(NA,JHIM-JBASM)) then
        IALO = min(IALO,max(1,JLOM-JBASM))
        IAHI = max(IAHI,min(NA,JHIM-JBASM))
      end if
    end do
    if (IAHI < IALO) cycle

    do IASTA=IALO,IAHI,NAMX
      IAEND = min(IASTA+NAMX-1,IAHI)
      NASZ = IAEND-IASTA+1
      ! SCR((v,l),a) = sum_P Cho_Ket(v,l)^P Cho_Bra(a,j)^P = (aj,vl)
      call DGEMM_('N','T',LDY,NASZ,NCHO,One,Cho_Ket,LDY,Cho_Bra(IASTA+NA*(IJ-1),1),NA*NJ,Zero,SCR,LDY)

      do IL=1,NL
        ILABS = IL+NIES(ISYL)
        if (ISYJ > ISYL) then
          JGEL = KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
          JGTL = KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
        else
          JGEL = KIGEJ(ILABS,IJABS)-NIGEJES(ISYJL)
          JGTL = KIGTJ(ILABS,IJABS)-NIGTJES(ISYJL)
        end if
        JBASP = IOFFP+NA*(JGEL-1)
        JBASM = IOFFM+NA*(JGTL-1)
        do IA=IASTA,IAEND
          ! plus combination
          if ((JBASP+IA >= JLOP) .and. (JBASP+IA <= JHIP)) &
            call DAXPY_(NV,SQH,SCR(1+NV*(IL-1)+LDY*(IA-IASTA)),1,WEP(1,JBASP+IA),1)
          ! minus combination
          if ((JBASM+IA >= JLOM) .and. (JBASM+IA <= JHIM)) &
            call DAXPY_(NV,SGN*SQ32,SCR(1+NV*(IL-1)+LDY*(IA-IASTA)),1,WEM(1,JBASM+IA),1)
        end do
      end do
    end do
  end do

end subroutine ADDRHSE_STRIPED_O

!-----------------------------------------------------------------------

subroutine ADDRHSF_STRIPED(JSYM,ISYU,ISYX,NA,NU,NC,NX, &
                           SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)

  use stdalloc, only: mma_allocate, mma_deallocate
  use caspt2_module, only: NAGEB, NAGTB, NINDEP, NTGEU, NTGTU
  use Symmetry_Info, only: Mul
  use SC_NEVPT2, only: Do_SC

  integer(kind=iwp), intent(in) :: JSYM, ISYU, ISYX, NA, NU, NC, NX, NSCR, NCHO
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NC,NX,NCHO)

  integer(kind=iwp) :: IA, IP, ISYA, ISYC, ISYM, IT, JHIM, JHIP, JLOM, JLOP, MOFFM, MOFFP, NACMX, NASM, NASP, NISM, NISP

  real(kind=wp), allocatable :: CHOBT(:)

  ! Case F: both combinations of (au,cx)
  !   rows    the active pair (u,x)
  !   columns the secondary pair (a,c)

  if (ISYU < ISYX) return
  ISYA = Mul(JSYM,ISYU)
  ISYC = Mul(JSYM,ISYX)
  ISYM = Mul(ISYU,ISYX)
  if (SKIP_SYM(ISYM)) return

  NASP = NTGEU(ISYM)
  NISP = NAGEB(ISYM)
  NASM = NTGTU(ISYM)
  NISM = NAGTB(ISYM)
  ! a combination with no independent parameters is not built, as in ADDRHSF
  if ((NINDEP(ISYM,8) == 0) .and. (.not. Do_SC)) NISP = 0
  if ((NINDEP(ISYM,9) == 0) .and. (.not. Do_SC)) NISM = 0
  if (NASP*NISP+NASM*NISM == 0) return

  ! the largest (a-block width) x (c columns) whose two integral blocks fit
  ! the scratch; the diagonal kernel splits it between the two dimensions
  NACMX = NSCR/(NU*NX+NU*NX)
  if ((NACMX < 1) .and. (ISYU == ISYX)) then
    write(u6,*) 'Not enough memory in ADDRHSF_STRIPED, I give up'
    call Abend()
  end if

  ! CHOBT(NU,NA,NCHO): transpose of Cho_Bra, a on the n side of the contraction
  ! without it the GEMM is only NU wide
  ! it is the ket operand as well, at ISYM == 1 the caller passes the same block as both
  call mma_allocate(CHOBT,NU*NA*NCHO,Label='CHOBT')
  do IP=1,NCHO
    do IA=1,NA
      do IT=1,NU
        CHOBT(IT+NU*(IA-1)+NU*NA*(IP-1)) = Cho_Bra(IA,IT,IP)
      end do
    end do
  end do

  call RHSLOC_BOUNDS(ISYM,8,NASP*NISP>0,MOFFP,JLOP,JHIP)
  call RHSLOC_BOUNDS(ISYM,9,NASM*NISM>0,MOFFM,JLOM,JHIM)

  ! the bra and the ket are the same block of Cholesky vectors when their symmetry labels agree
  if (ISYU == ISYX) then
    call ADDRHSF_STRIPED_D(RHSLoc(MOFFP),RHSLoc(MOFFM),NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYU,ISYX,ISYA,ISYM,NA,NU,NX, &
                           SCR,NSCR,NACMX,CHOBT,CHOBT,NCHO)
  else
    call ADDRHSF_STRIPED_O(RHSLoc(MOFFP),RHSLoc(MOFFM),NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYU,ISYX,ISYA,ISYC,ISYM,NA,NU,NC,NX, &
                           SCR,NSCR,CHOBT,Cho_Ket,NCHO)
  end if

  call mma_deallocate(CHOBT)

end subroutine ADDRHSF_STRIPED

!-----------------------------------------------------------------------

subroutine ADDRHSF_STRIPED_D(WFP,WFM,NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYU,ISYX,ISYA,ISYM,NA,NU,NX, &
                             SCR,NSCR,NACMX,CHOBT,CHOKT,NCHO)

  use SUPERINDEX, only: KAGEB, KAGTB, KTGEU, KTGTU
  use caspt2_module, only: NAES, NAGEBES, NAGTBES, NSES, NTGEUES, NTGTUES
  use Constants, only: Quart

  integer(kind=iwp), intent(in) :: NASP, NASM, JLOP, JHIP, JLOM, JHIM, ISYU, ISYX, ISYA, ISYM, NA, NU, NX, NSCR, NACMX, NCHO
  real(kind=wp), intent(inout) :: WFP(NASP,JLOP:JHIP), WFM(NASM,JLOM:JHIM)
  real(kind=wp), intent(out), target :: SCR(NSCR)
  real(kind=wp), intent(in) :: CHOBT(NU*NA,NCHO), CHOKT(NX*NA,NCHO)

  integer(kind=iwp) :: IA, IAABS, IAEND, IAL, IASTA, IC, ICBHI, ICBLO, ICEND, ICHI, ICHIM, ICHIMA(NAMXCAP), ICHIP, &
                       ICHIPA(NAMXCAP), ICL, ICLO, ICLOM, ICLOMA(NAMXCAP), ICLOP, ICLOPA(NAMXCAP), ICOLM, ICOLP, ICSTA, IRBASM, &
                       IRBASP, IX, IXABS, IY1, IY2, JBASM, JBASMA(NAMXCAP), JBASP, JBASPA(NAMXCAP), NAMX, NASZ, NCMX, NCSZ

  real(kind=wp), pointer, contiguous :: Y1(:,:,:,:) ! the (au,cx) block of this batch, as (u,a,x,c)
  real(kind=wp), pointer, contiguous :: Y2(:,:,:,:) ! the (cu,ax) block of this batch, as (u,c,x,a)

  ! Case F, diagonal block.
  !   WP(u>=x,a>=c) = ((au,cx)+(cu,ax))*(1-Kron(ux)/2)/2 * SQRT(1+Kron(ac))
  !   WM(u>x,a>c)   = ((cu,ax)-(au,cx))/2
  ! columns are secondary pairs, those of one a contiguous, the stripe is a range of c per a
  ! KTGEU/KTGTU(u,x) is contiguous in u for a fixed x
  ! Cho_Bra is (a,u,P), transposed to make "fixed a" a BLAS operand, the ket is the same block

  ! ISYU is kept for symmetry with the other _D kernels, the block is already
  ! selected by the caller
  unused_var(ISYU)

  ! a is blocked as in case G, one GEMM per block over the bounding c-rectangle
  ! the (x,c) panel is then read once per block, not once per a. See ADDRHSG_STRIPED_D
  NAMX = min(NA,NAMXCAP,NACMX)
  NCMX = max(NACMX/NAMX,1)

  do IASTA=1,NA,NAMX
    IAEND = min(IASTA+NAMX-1,NA)
    NASZ = IAEND-IASTA+1

    ! local c ranges of every a in the block (columns of a are the pairs
    ! (a,c), c = 1..a, consecutive) and their union
    ICBLO = NA+1
    ICBHI = 0
    do IA=IASTA,IAEND
      IAABS = IA+NSES(ISYA)
      ! the columns of c = 1..a (WFP) or 1..a-1 (WFM) are contiguous from JBASP/JBASM
      JBASP = KAGEB(IAABS,1+NSES(ISYA))-NAGEBES(ISYM)
      ICLOP = max(1,JLOP-JBASP+1)
      ICHIP = min(IA,JHIP-JBASP+1)
      if (ICLOP <= ICHIP) then
        ICBLO = min(ICBLO,ICLOP)
        ICBHI = max(ICBHI,ICHIP)
      end if
      JBASM = 0
      ICLOM = 1
      ICHIM = 0
      if ((IA >= 2) .and. (NASM > 0)) then
        JBASM = KAGTB(IAABS,1+NSES(ISYA))-NAGTBES(ISYM)
        ICLOM = max(1,JLOM-JBASM+1)
        ICHIM = min(IA-1,JHIM-JBASM+1)
        if (ICLOM <= ICHIM) then
          ICBLO = min(ICBLO,ICLOM)
          ICBHI = max(ICBHI,ICHIM)
        end if
      end if
      IAL = IA-IASTA+1
      JBASPA(IAL) = JBASP
      JBASMA(IAL) = JBASM
      ICLOPA(IAL) = ICLOP
      ICHIPA(IAL) = ICHIP
      ICLOMA(IAL) = ICLOM
      ICHIMA(IAL) = ICHIM
    end do
    if (ICBHI < ICBLO) cycle

    do ICSTA=ICBLO,ICBHI,NCMX
      ICEND = min(ICSTA+NCMX-1,ICBHI)
      NCSZ = ICEND-ICSTA+1
      ! SCR(1:) = Y1((u,a),(x,c)) = (au,cx); SCR(IY2:) = Y2((u,c),(x,a)) =
      ! (cu,ax)
      IY1 = 1
      IY2 = 1+NU*NASZ*NX*NCSZ
      call DGEMM_('N','T',NU*NASZ,NX*NCSZ,NCHO,One,CHOBT(1+NU*(IASTA-1),1),NU*NA,CHOKT(1+NX*(ICSTA-1),1),NX*NA, &
                  Zero,SCR(IY1),NU*NASZ)
      call DGEMM_('N','T',NU*NCSZ,NX*NASZ,NCHO,One,CHOBT(1+NU*(ICSTA-1),1),NU*NA,CHOKT(1+NX*(IASTA-1),1),NX*NA, &
                  Zero,SCR(IY2),NU*NCSZ)

      Y1(1:NU,1:NASZ,1:NX,1:NCSZ) => SCR(IY1:IY1+NU*NASZ*NX*NCSZ-1)
      Y2(1:NU,1:NCSZ,1:NX,1:NASZ) => SCR(IY2:IY2+NU*NCSZ*NX*NASZ-1)

      do IA=IASTA,IAEND
        IAL = IA-IASTA+1
        ICLOP = max(ICLOPA(IAL),ICSTA)
        ICHIP = min(ICHIPA(IAL),ICEND)
        ICLOM = max(ICLOMA(IAL),ICSTA)
        ICHIM = min(ICHIMA(IAL),ICEND)
        if ((ICHIP < ICLOP) .and. (ICHIM < ICLOM)) cycle
        ! ICLO:ICHI: union of the c ranges the two combinations need
        if (ICHIP < ICLOP) then
          ICLO = ICLOM
          ICHI = ICHIM
        else if (ICHIM < ICLOM) then
          ICLO = ICLOP
          ICHI = ICHIP
        else
          ICLO = min(ICLOP,ICLOM)
          ICHI = max(ICHIP,ICHIM)
        end if
        JBASP = JBASPA(IAL)
        JBASM = JBASMA(IAL)

        do IC=ICLO,ICHI
          ICL = IC-ICSTA+1
          ICOLP = JBASP+IC-1
          ICOLM = JBASM+IC-1
          do IX=1,NX
            IXABS = IX+NAES(ISYX)
            ! rows KTGEU(u,x) for u = x..NU are contiguous (first index
            ! fastest)
            IRBASP = KTGEU(IXABS,IXABS)-NTGEUES(ISYM)

            ! plus combination
            if ((IC >= ICLOP) .and. (IC <= ICHIP)) then
              if (IC == IA) then
                ! only one contribution, with the extra SQRT(2)
                WFP(IRBASP,ICOLP) = WFP(IRBASP,ICOLP)+SQ2*Quart*Y1(IX,IAL,IX,ICL)
                if (IX < NU) WFP(IRBASP+1:IRBASP+NU-IX,ICOLP) = WFP(IRBASP+1:IRBASP+NU-IX,ICOLP) &
                                                                +SQ2*Half*Y1(IX+1:NU,IAL,IX,ICL)
              else
                WFP(IRBASP,ICOLP) = WFP(IRBASP,ICOLP)+Quart*(Y1(IX,IAL,IX,ICL)+Y2(IX,ICL,IX,IAL))
                if (IX < NU) WFP(IRBASP+1:IRBASP+NU-IX,ICOLP) = WFP(IRBASP+1:IRBASP+NU-IX,ICOLP) &
                                                                +Half*(Y1(IX+1:NU,IAL,IX,ICL)+Y2(IX+1:NU,ICL,IX,IAL))
              end if
            end if

            ! minus combination
            if ((IC >= ICLOM) .and. (IC <= ICHIM) .and. (IX < NU)) then
              ! rows KTGTU(u,x) for u = x+1..NU are contiguous
              IRBASM = KTGTU(IXABS+1,IXABS)-NTGTUES(ISYM)
              WFM(IRBASM:IRBASM+NU-IX-1,ICOLM) = WFM(IRBASM:IRBASM+NU-IX-1,ICOLM) &
                                                 +Half*(Y2(IX+1:NU,ICL,IX,IAL)-Y1(IX+1:NU,IAL,IX,ICL))
            end if
          end do
        end do
      end do
    end do
  end do
  nullify(Y1,Y2)

end subroutine ADDRHSF_STRIPED_D

!-----------------------------------------------------------------------

subroutine ADDRHSF_STRIPED_O(WFP,WFM,NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYU,ISYX,ISYA,ISYC,ISYM,NA,NU,NC,NX, &
                             SCR,NSCR,CHOBT,Cho_Ket,NCHO)

  use SUPERINDEX, only: KAGEB, KAGTB, KTGEU, KTGTU
  use caspt2_module, only: NAES, NAGEBES, NAGTBES, NSES, NTGEUES, NTGTUES
  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), intent(in) :: NASP, NASM, JLOP, JHIP, JLOM, JHIM, ISYU, ISYX, ISYA, ISYC, ISYM, NA, NU, NC, NX, NSCR, NCHO
  real(kind=wp), intent(inout) :: WFP(NASP,JLOP:JHIP), WFM(NASM,JLOM:JHIM)
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: CHOBT(NU*NA,NCHO), Cho_Ket(NC*NX,NCHO)

  integer(kind=iwp) :: IA, IAABS, IAEND, IASTA, IC, ICABS, ICOLM, ICOLP, IP, IQ, IQHI, IQLO, IRBASM, IRBASP, IX, &
                       IXABS, JBASM, JBASP, LDY, NAMX, NASZ, NQ, NQSZ, NR
  real(kind=wp) :: SGN

  real(kind=wp), allocatable :: CHOKC(:,:)

  ! Case F, off-diagonal block: full rectangles, uniform Half weight
  ! a is batched into the n dimension via the transposed bra
  ! KTGEU(u,x) is contiguous in u for a fixed x
  ! q varies slowest in the columns, the stripe is a range of q against every r:
  !   NQ, NR       how many values q and r take
  !   IQLO:IQHI    the part of q the local stripe covers, NQSZ long
  ! a as q: a row range of the transposed bra
  ! c as q: the ket rows are not contiguous in (c,x), gathered into CHOKC

  SGN = -One
  if (ISYA <= ISYC) SGN = One

  if (ISYA > ISYC) then
    NQ = NA
    NR = NC
  else
    NQ = NC
    NR = NA
  end if

  IQLO = NQ+1
  IQHI = 0
  do IQ=1,NQ
    if (ISYA > ISYC) then
      JBASP = KAGEB(IQ+NSES(ISYA),1+NSES(ISYC))-NAGEBES(ISYM)
      JBASM = KAGTB(IQ+NSES(ISYA),1+NSES(ISYC))-NAGTBES(ISYM)
    else
      JBASP = KAGEB(IQ+NSES(ISYC),1+NSES(ISYA))-NAGEBES(ISYM)
      JBASM = KAGTB(IQ+NSES(ISYC),1+NSES(ISYA))-NAGTBES(ISYM)
    end if
    if ((JBASP+NR-1 >= JLOP) .and. (JBASP <= JHIP)) then
      IQLO = min(IQLO,IQ)
      IQHI = max(IQHI,IQ)
    end if
    if (NASM > 0) then
      if ((JBASM+NR-1 >= JLOM) .and. (JBASM <= JHIM)) then
        IQLO = min(IQLO,IQ)
        IQHI = max(IQHI,IQ)
      end if
    end if
  end do
  if (IQHI < IQLO) return

  if (ISYA > ISYC) then
    ! a is first, so its range is a row range of CHOBT: clamp the a loop to it
    LDY = NC*NX
    if (NSCR < LDY*NU) then
      write(u6,*) 'Not enough memory in ADDRHSF_STRIPED_O, I give up'
      call Abend()
    end if
    NAMX = NSCR/(LDY*NU)
    do IASTA=IQLO,IQHI,NAMX
      IAEND = min(IASTA+NAMX-1,IQHI)
      NASZ = IAEND-IASTA+1

      ! SCR((c,x),(u,a)) = sum_P Cho_Ket(c,x)^P CHOBT(u,a)^P = (au,cx)
      call DGEMM_('N','T',LDY,NASZ*NU,NCHO,One,Cho_Ket,LDY,CHOBT(1+NU*(IASTA-1),1),NU*NA,Zero,SCR,LDY)

      do IA=IASTA,IAEND
        IAABS = IA+NSES(ISYA)
        do IC=1,NC
          ICABS = IC+NSES(ISYC)
          ICOLP = KAGEB(IAABS,ICABS)-NAGEBES(ISYM)
          ICOLM = 0
          if (NASM > 0) ICOLM = KAGTB(IAABS,ICABS)-NAGTBES(ISYM)
          do IX=1,NX
            IXABS = IX+NAES(ISYX)
            ! plus combination
            if ((ICOLP >= JLOP) .and. (ICOLP <= JHIP)) then
              IRBASP = KTGEU(1+NAES(ISYU),IXABS)-NTGEUES(ISYM)
              call DAXPY_(NU,Half,SCR(IC+NC*(IX-1)+LDY*NU*(IA-IASTA)),LDY,WFP(IRBASP,ICOLP),1)
            end if
            ! minus combination
            if ((NASM > 0) .and. (ICOLM >= JLOM) .and. (ICOLM <= JHIM)) then
              IRBASM = KTGTU(1+NAES(ISYU),IXABS)-NTGTUES(ISYM)
              call DAXPY_(NU,SGN*Half,SCR(IC+NC*(IX-1)+LDY*NU*(IA-IASTA)),LDY,WFM(IRBASM,ICOLM),1)
            end if
          end do
        end do
      end do
    end do
  else
    ! c is first: gather the ket rows of its range into CHOKC((c,x),P), c
    ! fastest
    NQSZ = IQHI-IQLO+1
    call mma_allocate(CHOKC,NQSZ*NX,NCHO,Label='CHOKC')
    do IP=1,NCHO
      do IX=1,NX
        CHOKC(1+NQSZ*(IX-1):NQSZ*IX,IP) = Cho_Ket(IQLO+NC*(IX-1):IQHI+NC*(IX-1),IP)
      end do
    end do

    LDY = NQSZ*NX
    if (NSCR < LDY*NU) then
      write(u6,*) 'Not enough memory in ADDRHSF_STRIPED_O, I give up'
      call Abend()
    end if
    NAMX = NSCR/(LDY*NU)
    do IASTA=1,NA,NAMX
      IAEND = min(IASTA+NAMX-1,NA)
      NASZ = IAEND-IASTA+1

      ! SCR((c,x),(u,a)) = sum_P CHOKC(c,x)^P CHOBT(u,a)^P = (au,cx)
      call DGEMM_('N','T',LDY,NASZ*NU,NCHO,One,CHOKC,LDY,CHOBT(1+NU*(IASTA-1),1),NU*NA,Zero,SCR,LDY)

      do IA=IASTA,IAEND
        IAABS = IA+NSES(ISYA)
        do IC=1,NQSZ
          ICABS = IQLO+IC-1+NSES(ISYC)
          ICOLP = KAGEB(ICABS,IAABS)-NAGEBES(ISYM)
          ICOLM = 0
          if (NASM > 0) ICOLM = KAGTB(ICABS,IAABS)-NAGTBES(ISYM)
          do IX=1,NX
            IXABS = IX+NAES(ISYX)
            ! plus combination
            if ((ICOLP >= JLOP) .and. (ICOLP <= JHIP)) then
              IRBASP = KTGEU(1+NAES(ISYU),IXABS)-NTGEUES(ISYM)
              call DAXPY_(NU,Half,SCR(IC+NQSZ*(IX-1)+LDY*NU*(IA-IASTA)),LDY,WFP(IRBASP,ICOLP),1)
            end if
            ! minus combination
            if ((NASM > 0) .and. (ICOLM >= JLOM) .and. (ICOLM <= JHIM)) then
              IRBASM = KTGTU(1+NAES(ISYU),IXABS)-NTGTUES(ISYM)
              call DAXPY_(NU,SGN*Half,SCR(IC+NQSZ*(IX-1)+LDY*NU*(IA-IASTA)),LDY,WFM(IRBASM,ICOLM),1)
            end if
          end do
        end do
      end do
    end do

    call mma_deallocate(CHOKC)
  end if

end subroutine ADDRHSF_STRIPED_O

!-----------------------------------------------------------------------

subroutine ADDRHSG_STRIPED(JSYM,ISYU,ISYL,NA,NU,NC,NL, &
                           SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)

  use stdalloc, only: mma_allocate, mma_deallocate
  use caspt2_module, only: NISUP
  use general_data, only: NASH
  use Symmetry_Info, only: Mul

  integer(kind=iwp), intent(in) :: JSYM, ISYU, ISYL, NA, NU, NC, NL, NSCR, NCHO
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NC*NL,NCHO)

  integer(kind=iwp) :: IA, IP, ISYA, ISYAC, ISYC, ISYM, JHIM, JHIP, JLOM, JLOP, MOFFM, MOFFP, NACMX, NAS, NISM, NISP

  real(kind=wp), allocatable :: CHOBT(:)

  ! Case G: both combinations of (au,cl)
  !   rows    the active index alone
  !   columns an inactive index within a secondary pair

  ISYA = Mul(JSYM,ISYU)
  ISYC = Mul(JSYM,ISYL)
  ISYM = ISYU
  if (SKIP_SYM(ISYM)) return
  ISYAC = Mul(ISYA,ISYC)
  NAS = NASH(ISYM)
  NISP = NISUP(ISYM,10)
  NISM = NISUP(ISYM,11)
  if (NAS*(NISP+NISM) == 0) return

  ! the largest (a-block width) x (c columns) whose two integral blocks fit
  ! the scratch; the diagonal kernel splits it between the two dimensions
  NACMX = NSCR/(2*NU*NL)
  if ((NACMX < 1) .and. (ISYU == ISYL)) then
    write(u6,*) 'Not enough memory in ADDRHSG_STRIPED, I give up'
    call Abend()
  end if

  ! CHOBT puts a on the n side of the contraction; without it the GEMM is only NU wide.
  call mma_allocate(CHOBT,NU*NA*NCHO,Label='CHOBT')
  do IP=1,NCHO
    do IA=1,NA
      CHOBT(1+NU*(IA-1)+NU*NA*(IP-1):NU+NU*(IA-1)+NU*NA*(IP-1)) = Cho_Bra(IA,1:NU,IP)
    end do
  end do

  call RHSLOC_BOUNDS(ISYM,10,NAS*NISP>0,MOFFP,JLOP,JHIP)
  call RHSLOC_BOUNDS(ISYM,11,NAS*NISM>0,MOFFM,JLOM,JHIM)

  ! the bra and the ket are the same block of Cholesky vectors when their symmetry labels agree
  if (ISYU == ISYL) then
    call ADDRHSG_STRIPED_D(RHSLoc(MOFFP),RHSLoc(MOFFM),NAS,JLOP,JHIP,JLOM,JHIM,ISYL,ISYA,ISYM,ISYAC,NA,NU,NL, &
                           SCR,NSCR,NACMX,CHOBT,Cho_Ket,NCHO)
  else
    call ADDRHSG_STRIPED_O(RHSLoc(MOFFP),RHSLoc(MOFFM),NAS,JLOP,JHIP,JLOM,JHIM,ISYL,ISYA,ISYC,ISYM,ISYAC,NA,NU,NC,NL, &
                           SCR,NSCR,CHOBT,Cho_Ket,NCHO)
  end if

  call mma_deallocate(CHOBT)

end subroutine ADDRHSG_STRIPED

!-----------------------------------------------------------------------

subroutine ADDRHSG_STRIPED_D(WGP,WGM,NAS,JLOP,JHIP,JLOM,JHIM,ISYL,ISYA,ISYM,ISYAC,NA,NU,NL, &
                             SCR,NSCR,NACMX,CHOBT,CHOKT,NCHO)

  use SUPERINDEX, only: KAGEB, KAGTB
  use caspt2_module, only: NAGEB, NAGEBES, NAGTB, NAGTBES, NISH, NSES, NSYM
  use Symmetry_Info, only: Mul

  integer(kind=iwp), intent(in) :: NAS, JLOP, JHIP, JLOM, JHIM, ISYL, ISYA, ISYM, ISYAC, NA, NU, NL, NSCR, NACMX, NCHO
  real(kind=wp), intent(inout) :: WGP(NAS,JLOP:JHIP), WGM(NAS,JLOM:JHIM)
  real(kind=wp), intent(out), target :: SCR(NSCR)
  real(kind=wp), intent(in) :: CHOBT(NU*NA,NCHO), CHOKT(NL*NA,NCHO)

  integer(kind=iwp) :: IA, IAABS, IAEND, IAL, IASTA, IC, ICBHI, ICBLO, ICEND, ICHI, ICHIM, ICHIMA(NAMXCAP), ICHIP, &
                       ICHIPA(NAMXCAP), ICL, ICLO, ICLOM, ICLOMA(NAMXCAP), ICLOP, ICLOPA(NAMXCAP), ICOL, ICSTA, IL, IOFFM, IOFFP, &
                       ISAB, ISI, IY1, IY2, JBASM, JBASMA(NAMXCAP), JBASP, JBASPA(NAMXCAP), JCOLHI, JCOLLO, NAMX, NASZ, NCMX, NCSZ

  real(kind=wp), pointer, contiguous :: Y1(:,:,:,:) ! the (au,cl) block of this batch, as (u,a,l,c)
  real(kind=wp), pointer, contiguous :: Y2(:,:,:,:) ! the (cu,al) block of this batch, as (u,c,l,a)

  ! Case G, diagonal block.
  !   WP(u,l,a>=c) = ((au,cl)+(cu,al))/SQRT(2+2*Kron(ac))
  !   WM(u,l,a>c)  = ((au,cl)-(cu,al))*SQRT(OneHalf)
  ! columns run l fastest within a secondary pair, the pairs of one a are consecutive
  ! one a is then a contiguous column range to intersect with the stripe
  ! Cho_Bra is (a,u,P), transposed here
  ! Cho_Ket arrives as (l,c,P) from TRANSPOSE_KET, which makes a range of c contiguous
  ! that is done there and not here, the same ket block is revisited for every bra symmetry

  ! offsets of this inactive symmetry within the two column index sets
  IOFFP = 0
  IOFFM = 0
  do ISI=1,NSYM
    if (ISI == ISYL) exit
    ISAB = Mul(ISI,ISYM)
    IOFFP = IOFFP+NISH(ISI)*NAGEB(ISAB)
    IOFFM = IOFFM+NISH(ISI)*NAGTB(ISAB)
  end do

  ! a is blocked, the (l,c) panel is read once per block instead of once per a
  ! the block spans [ICBLO,ICBHI] though row a needs only a >= c
  ! the extra c > a columns computed this way are about NAMX/(2a) of the block
  ! NACMX = NSCR/(2*NU*NL) is the largest NAMX*NCMX that fits
  NAMX = min(NA,NAMXCAP,NACMX)
  NCMX = max(NACMX/NAMX,1)

  do IASTA=1,NA,NAMX
    IAEND = min(IASTA+NAMX-1,NA)
    NASZ = IAEND-IASTA+1

    ! columns of a in WGP: (p0+c-1, l), c = 1..IA, l fastest; the stripe
    ! [JLOP,JHIP] selects the local c range of every a, ICBLO:ICBHI is the
    ! union over the block
    ICBLO = NA+1
    ICBHI = 0
    do IA=IASTA,IAEND
      IAABS = IA+NSES(ISYA)
      ! the columns of this a are (pair JBASP+c-1, l), c = 1..a, l fastest
      ! JCOLLO:JCOLHI: that flat column range clipped to the local stripe
      ! ICLOP:ICHIP (ICLOM:ICHIM): the same range converted to c bounds
      JBASP = KAGEB(IAABS,1+NSES(ISYA))-NAGEBES(ISYAC)
      JCOLLO = max(IOFFP+(JBASP-1)*NL+1,JLOP)
      JCOLHI = min(IOFFP+(JBASP+IA-1)*NL,JHIP)
      if (JCOLLO <= JCOLHI) then
        ICLOP = (JCOLLO-IOFFP-(JBASP-1)*NL-1)/NL+1
        ICHIP = (JCOLHI-IOFFP-(JBASP-1)*NL-1)/NL+1
        ICBLO = min(ICBLO,ICLOP)
        ICBHI = max(ICBHI,ICHIP)
      else
        ICLOP = 1
        ICHIP = 0
      end if
      if (IA >= 2) then
        JBASM = KAGTB(IAABS,1+NSES(ISYA))-NAGTBES(ISYAC)
        JCOLLO = max(IOFFM+(JBASM-1)*NL+1,JLOM)
        JCOLHI = min(IOFFM+(JBASM+IA-2)*NL,JHIM)
        if (JCOLLO <= JCOLHI) then
          ICLOM = (JCOLLO-IOFFM-(JBASM-1)*NL-1)/NL+1
          ICHIM = (JCOLHI-IOFFM-(JBASM-1)*NL-1)/NL+1
          ICBLO = min(ICBLO,ICLOM)
          ICBHI = max(ICBHI,ICHIM)
        else
          ICLOM = 1
          ICHIM = 0
        end if
      else
        JBASM = 0
        ICLOM = 1
        ICHIM = 0
      end if
      IAL = IA-IASTA+1
      JBASPA(IAL) = JBASP
      JBASMA(IAL) = JBASM
      ICLOPA(IAL) = ICLOP
      ICHIPA(IAL) = ICHIP
      ICLOMA(IAL) = ICLOM
      ICHIMA(IAL) = ICHIM
    end do
    if (ICBHI < ICBLO) cycle

    do ICSTA=ICBLO,ICBHI,NCMX
      ICEND = min(ICSTA+NCMX-1,ICBHI)
      NCSZ = ICEND-ICSTA+1
      ! SCR(1:) = Y1((u,a),(l,c)) = (au,cl); SCR(IY2:) = Y2((u,c),(l,a)) =
      ! (cu,al)
      IY1 = 1
      IY2 = 1+NU*NASZ*NL*NCSZ
      call DGEMM_('N','T',NU*NASZ,NL*NCSZ,NCHO,One,CHOBT(1+NU*(IASTA-1),1),NU*NA,CHOKT(1+NL*(ICSTA-1),1),NL*NA, &
                  Zero,SCR(IY1),NU*NASZ)
      call DGEMM_('N','T',NU*NCSZ,NL*NASZ,NCHO,One,CHOBT(1+NU*(ICSTA-1),1),NU*NA,CHOKT(1+NL*(IASTA-1),1),NL*NA, &
                  Zero,SCR(IY2),NU*NCSZ)

      Y1(1:NU,1:NASZ,1:NL,1:NCSZ) => SCR(IY1:IY1+NU*NASZ*NL*NCSZ-1)
      Y2(1:NU,1:NCSZ,1:NL,1:NASZ) => SCR(IY2:IY2+NU*NCSZ*NL*NASZ-1)

      do IA=IASTA,IAEND
        IAL = IA-IASTA+1
        ICLOP = max(ICLOPA(IAL),ICSTA)
        ICHIP = min(ICHIPA(IAL),ICEND)
        ICLOM = max(ICLOMA(IAL),ICSTA)
        ICHIM = min(ICHIMA(IAL),ICEND)
        if ((ICHIP < ICLOP) .and. (ICHIM < ICLOM)) cycle
        ! ICLO:ICHI: union of the c ranges the two combinations need
        if (ICHIP < ICLOP) then
          ICLO = ICLOM
          ICHI = ICHIM
        else if (ICHIM < ICLOM) then
          ICLO = ICLOP
          ICHI = ICHIP
        else
          ICLO = min(ICLOP,ICLOM)
          ICHI = max(ICHIP,ICHIM)
        end if
        JBASP = JBASPA(IAL)
        JBASM = JBASMA(IAL)

        do IC=ICLO,ICHI
          ICL = IC-ICSTA+1
          do IL=1,NL
            ! plus combination
            ICOL = IOFFP+(JBASP+IC-2)*NL+IL
            if ((IC >= ICLOP) .and. (IC <= ICHIP) .and. (ICOL >= JLOP) .and. (ICOL <= JHIP)) then
              if (IC == IA) then
                WGP(1:NU,ICOL) = WGP(1:NU,ICOL)+Y1(:,IAL,IL,ICL)
              else
                WGP(1:NU,ICOL) = WGP(1:NU,ICOL)+SQH*(Y1(:,IAL,IL,ICL)+Y2(:,ICL,IL,IAL))
              end if
            end if
            ! minus combination, a > c only
            if (IC < IA) then
              ICOL = IOFFM+(JBASM+IC-2)*NL+IL
              if ((IC >= ICLOM) .and. (IC <= ICHIM) .and. (ICOL >= JLOM) .and. (ICOL <= JHIM)) &
                WGM(1:NU,ICOL) = WGM(1:NU,ICOL)+SQ32*(Y1(:,IAL,IL,ICL)-Y2(:,ICL,IL,IAL))
            end if
          end do
        end do
      end do

    end do
  end do
  nullify(Y1,Y2)

end subroutine ADDRHSG_STRIPED_D

!-----------------------------------------------------------------------

subroutine ADDRHSG_STRIPED_O(WGP,WGM,NAS,JLOP,JHIP,JLOM,JHIM,ISYL,ISYA,ISYC,ISYM,ISYAC,NA,NU,NC,NL, &
                             SCR,NSCR,CHOBT,CHOKT,NCHO)

  use SUPERINDEX, only: KAGEB, KAGTB
  use caspt2_module, only: NAGEB, NAGEBES, NAGTB, NAGTBES, NISH, NSES, NSYM
  use Symmetry_Info, only: Mul

  integer(kind=iwp), intent(in) :: NAS, JLOP, JHIP, JLOM, JHIM, ISYL, ISYA, ISYC, ISYM, ISYAC, NA, NU, NC, NL, NSCR, NCHO
  real(kind=wp), intent(inout) :: WGP(NAS,JLOP:JHIP), WGM(NAS,JLOM:JHIM)
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: CHOBT(NU*NA,NCHO), CHOKT(NL*NC,NCHO)

  integer(kind=iwp) :: IA, IAABS, IAGEC, IAGTC, IAHI, IALO, IC, ICABS, ICEND, ICOLM, ICOLP, ICSTA, IL, IOFFM, &
                       IOFFP, IQ, IQEND, IQHI, IQLO, IQSTA, ISAB, ISI, JBASM, JBASP, LDY, NASZ, NCSZ, NQ, NQMX, NR
  real(kind=wp) :: SGN

  ! Case G, off-diagonal block: uniform weight, same flow as case H once both blocks are transposed
  ! q varies slowest in the columns, the stripe is a range of it
  ! the roles of the two operands swap with q
  !   NQ, NR       how many values q and r take
  !   IQLO:IQHI    the part of q the local stripe covers

  SGN = One
  if (ISYA <= ISYC) SGN = -One

  IOFFP = 0
  IOFFM = 0
  do ISI=1,NSYM
    if (ISI == ISYL) exit
    ISAB = Mul(ISI,ISYM)
    IOFFP = IOFFP+NISH(ISI)*NAGEB(ISAB)
    IOFFM = IOFFM+NISH(ISI)*NAGTB(ISAB)
  end do

  if (ISYA > ISYC) then
    NQ = NA
    NR = NC
  else
    NQ = NC
    NR = NA
  end if
  IQLO = NQ+1
  IQHI = 0
  do IQ=1,NQ
    if (ISYA > ISYC) then
      JBASP = KAGEB(IQ+NSES(ISYA),1+NSES(ISYC))-NAGEBES(ISYAC)
      JBASM = KAGTB(IQ+NSES(ISYA),1+NSES(ISYC))-NAGTBES(ISYAC)
    else
      JBASP = KAGEB(IQ+NSES(ISYC),1+NSES(ISYA))-NAGEBES(ISYAC)
      JBASM = KAGTB(IQ+NSES(ISYC),1+NSES(ISYA))-NAGTBES(ISYAC)
    end if
    if ((IOFFP+NL*(JBASP+NR-1) >= JLOP) .and. (IOFFP+NL*(JBASP-1)+1 <= JHIP)) then
      IQLO = min(IQLO,IQ)
      IQHI = max(IQHI,IQ)
    end if
    if ((IOFFM+NL*(JBASM+NR-1) >= JLOM) .and. (IOFFM+NL*(JBASM-1)+1 <= JHIM)) then
      IQLO = min(IQLO,IQ)
      IQHI = max(IQHI,IQ)
    end if
  end do
  if (IQHI < IQLO) return

  if (NSCR < NL*NR*NU) then
    write(u6,*) 'Not enough memory in ADDRHSG_STRIPED_O, I give up'
    call Abend()
  end if
  NQMX = NSCR/(NL*NR*NU)

  do IQSTA=IQLO,IQHI,NQMX
    IQEND = min(IQSTA+NQMX-1,IQHI)
    ! IALO:IAHI and ICSTA:ICEND: the a and c ranges this block covers. Only the
    ! q member is blocked, the other one is taken whole
    if (ISYA > ISYC) then
      IALO = IQSTA
      IAHI = IQEND
      ICSTA = 1
      ICEND = NC
      NASZ = IAHI-IALO+1
      NCSZ = NC
      LDY = NL*NCSZ
      call DGEMM_('N','T',LDY,NASZ*NU,NCHO,One,CHOKT,NL*NC,CHOBT(1+NU*(IQSTA-1),1),NU*NA,Zero,SCR,LDY)
    else
      IALO = 1
      IAHI = NA
      ICSTA = IQSTA
      ICEND = IQEND
      NASZ = IAHI-IALO+1
      NCSZ = ICEND-ICSTA+1
      LDY = NL*NCSZ
      call DGEMM_('N','T',LDY,NASZ*NU,NCHO,One,CHOKT(1+NL*(ICSTA-1),1),NL*NC,CHOBT,NU*NA,Zero,SCR,LDY)
    end if

    do IA=IALO,IAHI
      IAABS = IA+NSES(ISYA)
      do IC=ICSTA,ICEND
        ICABS = IC+NSES(ISYC)
        if (ISYA > ISYC) then
          IAGEC = KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
          IAGTC = KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
        else
          IAGEC = KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
          IAGTC = KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
        end if
        do IL=1,NL
          ! plus combination
          ICOLP = IOFFP+NL*(IAGEC-1)+IL
          if ((ICOLP >= JLOP) .and. (ICOLP <= JHIP)) &
            call DAXPY_(NU,SQH,SCR(IL+NL*(IC-ICSTA)+LDY*NU*(IA-IALO)),LDY,WGP(1,ICOLP),1)
          ! minus combination
          ICOLM = IOFFM+NL*(IAGTC-1)+IL
          if ((ICOLM >= JLOM) .and. (ICOLM <= JHIM)) &
            call DAXPY_(NU,SGN*SQ32,SCR(IL+NL*(IC-ICSTA)+LDY*NU*(IA-IALO)),LDY,WGM(1,ICOLM),1)
        end do
      end do
    end do
  end do

end subroutine ADDRHSG_STRIPED_O

!-----------------------------------------------------------------------

subroutine ADDRHSH_STRIPED(JSYM,ISYJ,ISYL,NA,NJ,NC,NL, &
                           SCR,NSCR,Cho_Bra,Cho_Ket,NCHO)

  use caspt2_module, only: NAGEB, NAGTB, NIGEJ, NIGTJ
  use Symmetry_Info, only: Mul

  integer(kind=iwp), intent(in) :: JSYM, ISYJ, ISYL, NA, NJ, NC, NL, NSCR, NCHO
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA*NJ*NCHO), Cho_Ket(NC*NL*NCHO)

  integer(kind=iwp) :: ISYA, ISYC, ISYM, JHIM, JHIP, JLOM, JLOP, MOFFM, MOFFP, NASM, NASP, NISM, NISP, NLMX

  ! Case H: both combinations of (aj,cl)
  !   rows    the secondary pair (a,c)
  !   columns the inactive pair (j,l)

  if (ISYJ < ISYL) return
  ISYA = Mul(JSYM,ISYJ)
  ISYC = Mul(JSYM,ISYL)
  ISYM = Mul(ISYA,ISYC)
  if (SKIP_SYM(ISYM)) return
  NASP = NAGEB(ISYM)
  NISP = NIGEJ(ISYM)
  NASM = NAGTB(ISYM)
  NISM = NIGTJ(ISYM)
  if (NASP*NISP == 0) return

  ! largest number of ket columns (l) whose integral block fits the scratch
  NLMX = NSCR/(NA*NC)
  if (NLMX < 1) then
    write(u6,*) 'Not enough memory in ADDRHSH_STRIPED, I give up'
    call Abend()
  end if

  call RHSLOC_BOUNDS(ISYM,12,NASP*NISP>0,MOFFP,JLOP,JHIP)
  call RHSLOC_BOUNDS(ISYM,13,NASM*NISM>0,MOFFM,JLOM,JHIM)

  ! the bra and the ket are the same block of Cholesky vectors when their symmetry labels agree
  if (ISYJ == ISYL) then
    call ADDRHSH_STRIPED_D(RHSLoc(MOFFP),RHSLoc(MOFFM),NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYA,ISYM,NA,NJ, &
                           SCR,NSCR,NLMX,Cho_Bra,NCHO)
  else
    call ADDRHSH_STRIPED_O(RHSLoc(MOFFP),RHSLoc(MOFFM),NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYL,ISYA,ISYC,ISYM,NA,NJ,NC,NL, &
                           SCR,NSCR,NLMX,Cho_Bra,Cho_Ket,NCHO)
  end if

end subroutine ADDRHSH_STRIPED

!-----------------------------------------------------------------------

subroutine ADDRHSH_STRIPED_D(WHP,WHM,NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYA,ISYM,NA,NJ, &
                             SCR,NSCR,NLMX,Cho_Bra,NCHO)

  use SUPERINDEX, only: KAGEB, KAGTB, KIGEJ, KIGTJ
  use caspt2_module, only: NAGEBES, NAGTBES, NIES, NIGEJES, NIGTJES, NSES

  integer(kind=iwp), intent(in) :: NASP, NASM, JLOP, JHIP, JLOM, JHIM, ISYJ, ISYA, ISYM, NA, NJ, NSCR, NLMX, NCHO
  real(kind=wp), intent(inout) :: WHP(NASP,JLOP:JHIP), WHM(NASM,JLOM:JHIM)
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA*NJ,NCHO)

  integer(kind=iwp) :: IA, IAABS, ICOL, IJ, IJABS, IL, ILEND, ILHI, ILHIM, ILHIP, ILLO, ILLOM, ILLOP, ILSTA, IOFF, IRBASM, &
                       IRBASP, JBASM, JBASP, LDY, NLSZ
  real(kind=wp) :: SCL

  ! Case H, diagonal block
  !   WP(a>=c,j>=l) = ((aj,cl)+(al,cj))/SQRT((1+Kron(jl))*(1+Kron(ac)))
  !   WM(a>c,j>l)   = ((aj,cl)-(al,cj))*SQRT(Three)
  ! both combinations come from one integral block, over the union of the columns
  ! their two (differently distributed) global arrays own locally
  ! only the l <= j half is formed, the replicated algorithm computes the other half and drops it

  do IJ=1,NJ
    IJABS = IJ+NIES(ISYJ)
    ! the columns of l = 1..j (WHP) or 1..j-1 (WHM) are contiguous from JBASP/JBASM
    JBASP = KIGEJ(IJABS,1+NIES(ISYJ))-NIGEJES(ISYM)
    ILLOP = max(1,JLOP-JBASP+1)
    ILHIP = min(IJ,JHIP-JBASP+1)
    JBASM = 0
    ILLOM = 1
    ILHIM = 0
    if ((IJ >= 2) .and. (NASM > 0)) then
      JBASM = KIGTJ(IJABS,1+NIES(ISYJ))-NIGTJES(ISYM)
      ILLOM = max(1,JLOM-JBASM+1)
      ILHIM = min(IJ-1,JHIM-JBASM+1)
    end if

    if ((ILHIP < ILLOP) .and. (ILHIM < ILLOM)) cycle

    ! ILLO:ILHI: union of the l ranges the two combinations need
    if (ILHIP < ILLOP) then ! only minus
      ILLO = ILLOM
      ILHI = ILHIM
    else if (ILHIM < ILLOM) then ! only plus
      ILLO = ILLOP
      ILHI = ILHIP
    else ! both
      ILLO = min(ILLOP,ILLOM)
      ILHI = max(ILHIP,ILHIM)
    end if

    do ILSTA=ILLO,ILHI,NLMX
      ILEND = min(ILSTA+NLMX-1,ILHI)
      NLSZ = ILEND-ILSTA+1
      LDY = NA*NLSZ
      ! SCR(cl,a) = sum_P Cho_Bra(c,l)^P Cho_Bra(a,j)^P = (aj,cl)
      call DGEMM_('N','T',LDY,NA,NCHO,One,Cho_Bra(1+NA*(ILSTA-1),1),NA*NJ,Cho_Bra(1+NA*(IJ-1),1),NA*NJ,Zero,SCR,LDY)

      do IL=ILSTA,ILEND
        IOFF = NA*(IL-ILSTA)

        ! plus combination
        if ((IL >= ILLOP) .and. (IL <= ILHIP)) then
          ICOL = JBASP+IL-1
          SCL = One
          if (IL == IJ) SCL = SQH
          do IA=1,NA
            IAABS = IA+NSES(ISYA)
            ! the rows KAGEB(a,c), c = 1..a, are contiguous
            IRBASP = KAGEB(IAABS,1+NSES(ISYA))-NAGEBES(ISYM)
            if (IA >= 2) then
              call DAXPY_(IA-1,SCL,SCR(IOFF+1+LDY*(IA-1)),1,WHP(IRBASP,ICOL),1)
              call DAXPY_(IA-1,SCL,SCR(IOFF+IA),LDY,WHP(IRBASP,ICOL),1)
            end if
            WHP(IRBASP+IA-1,ICOL) = WHP(IRBASP+IA-1,ICOL)+SQ2*SCL*SCR(IOFF+IA+LDY*(IA-1))
          end do
        end if

        ! minus combination
        if ((IL >= ILLOM) .and. (IL <= ILHIM)) then
          ICOL = JBASM+IL-1
          do IA=2,NA
            IAABS = IA+NSES(ISYA)
            ! the rows KAGTB(a,c), c = 1..a-1, are contiguous
            IRBASM = KAGTB(IAABS,1+NSES(ISYA))-NAGTBES(ISYM)
            call DAXPY_(IA-1, SQ3,SCR(IOFF+1+LDY*(IA-1)),1,WHM(IRBASM,ICOL),1)
            call DAXPY_(IA-1,-SQ3,SCR(IOFF+IA),LDY,WHM(IRBASM,ICOL),1)
          end do
        end if

      end do
    end do
  end do

end subroutine ADDRHSH_STRIPED_D

!-----------------------------------------------------------------------

subroutine ADDRHSH_STRIPED_O(WHP,WHM,NASP,NASM,JLOP,JHIP,JLOM,JHIM,ISYJ,ISYL,ISYA,ISYC,ISYM,NA,NJ,NC,NL, &
                             SCR,NSCR,NLMX,Cho_Bra,Cho_Ket,NCHO)

  use SUPERINDEX, only: KAGEB, KAGTB, KIGEJ, KIGTJ
  use caspt2_module, only: NAGEBES, NAGTBES, NIES, NIGEJES, NIGTJES, NSES

  integer(kind=iwp), intent(in) :: NASP, NASM, JLOP, JHIP, JLOM, JHIM, ISYJ, ISYL, ISYA, ISYC, ISYM, NA, NJ, NC, NL, &
                                   NSCR, NLMX, NCHO
  real(kind=wp), intent(inout) :: WHP(NASP,JLOP:JHIP), WHM(NASM,JLOM:JHIM)
  real(kind=wp), intent(out) :: SCR(NSCR)
  real(kind=wp), intent(in) :: Cho_Bra(NA*NJ,NCHO), Cho_Ket(NC*NL,NCHO)

  integer(kind=iwp) :: IA, IAABS, IC, ICABS, ICOL, IJ, IJABS, IL, ILEND, ILHI, ILHIM, ILHIP, ILLO, ILLOM, ILLOP, ILSTA, IOFF, &
                       IRBASM, IRBASP, JBASM, JBASP, LDY, NLSZ
  real(kind=wp) :: SGN

  ! Case H, off-diagonal block
  ! the partner term comes from the JSYM = JSYM*ISYM iteration with the secondary symmetries exchanged
  !   rows  KAGEB(a,c) or KAGEB(c,a), contiguous in the ket/bra secondary
  !   cols  KIGEJ(j,l), contiguous in l for a fixed j
  ! ISYA > ISYC picks between the two row forms
  ! unlike F and G the rows are not striped, only the run direction changes

  SGN = One
  if (ISYA <= ISYC) SGN = -One

  do IJ=1,NJ
    IJABS = IJ+NIES(ISYJ)
    ! the columns of l = 1..NL are contiguous from JBASP/JBASM
    JBASP = KIGEJ(IJABS,1+NIES(ISYL))-NIGEJES(ISYM)
    ILLOP = max(1,JLOP-JBASP+1)
    ILHIP = min(NL,JHIP-JBASP+1)
    JBASM = 0
    ILLOM = 1
    ILHIM = 0
    if (NASM > 0) then
      JBASM = KIGTJ(IJABS,1+NIES(ISYL))-NIGTJES(ISYM)
      ILLOM = max(1,JLOM-JBASM+1)
      ILHIM = min(NL,JHIM-JBASM+1)
    end if

    if ((ILHIP < ILLOP) .and. (ILHIM < ILLOM)) cycle

    ! ILLO:ILHI: union of the l ranges the two combinations need
    if (ILHIP < ILLOP) then
      ILLO = ILLOM
      ILHI = ILHIM
    else if (ILHIM < ILLOM) then
      ILLO = ILLOP
      ILHI = ILHIP
    else
      ILLO = min(ILLOP,ILLOM)
      ILHI = max(ILHIP,ILHIM)
    end if

    do ILSTA=ILLO,ILHI,NLMX
      ILEND = min(ILSTA+NLMX-1,ILHI)
      NLSZ = ILEND-ILSTA+1
      LDY = NC*NLSZ
      ! SCR(cl,a) = sum_P Cho_Ket(c,l)^P Cho_Bra(a,j)^P = (aj,cl)
      call DGEMM_('N','T',LDY,NA,NCHO,One,Cho_Ket(1+NC*(ILSTA-1),1),NC*NL,Cho_Bra(1+NA*(IJ-1),1),NA*NJ,Zero,SCR,LDY)

      do IL=ILSTA,ILEND
        IOFF = NC*(IL-ILSTA)
        if (ISYA > ISYC) then
          ! rows KAGEB(a,c): contiguous in c for a fixed a
          do IA=1,NA
            IAABS = IA+NSES(ISYA)
            ! plus combination
            if ((IL >= ILLOP) .and. (IL <= ILHIP)) then
              ICOL = JBASP+IL-1
              IRBASP = KAGEB(IAABS,1+NSES(ISYC))-NAGEBES(ISYM)
              call DAXPY_(NC,One,SCR(IOFF+1+LDY*(IA-1)),1,WHP(IRBASP,ICOL),1)
            end if
            ! minus combination
            if ((NASM > 0) .and. (IL >= ILLOM) .and. (IL <= ILHIM)) then
              ICOL = JBASM+IL-1
              IRBASM = KAGTB(IAABS,1+NSES(ISYC))-NAGTBES(ISYM)
              call DAXPY_(NC,SGN*SQ3,SCR(IOFF+1+LDY*(IA-1)),1,WHM(IRBASM,ICOL),1)
            end if
          end do
        else
          ! rows KAGEB(c,a): contiguous in a for a fixed c
          do IC=1,NC
            ICABS = IC+NSES(ISYC)
            ! plus combination
            if ((IL >= ILLOP) .and. (IL <= ILHIP)) then
              ICOL = JBASP+IL-1
              IRBASP = KAGEB(ICABS,1+NSES(ISYA))-NAGEBES(ISYM)
              call DAXPY_(NA,One,SCR(IOFF+IC),LDY,WHP(IRBASP,ICOL),1)
            end if
            ! minus combination
            if ((NASM > 0) .and. (IL >= ILLOM) .and. (IL <= ILHIM)) then
              ICOL = JBASM+IL-1
              IRBASM = KAGTB(ICABS,1+NSES(ISYA))-NAGTBES(ISYM)
              call DAXPY_(NA,SGN*SQ3,SCR(IOFF+IC),LDY,WHM(IRBASM,ICOL),1)
            end if
          end do
        end if
      end do
    end do
  end do

end subroutine ADDRHSH_STRIPED_O

end module ADDRHS_STRIPED

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(ADDRHS_STRIPED)

#endif
