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

subroutine RHSALL2_STRIPED(IVEC)

! Striped construction of the CASPT2 right-hand side, selected with PRHS = 4 or STRIPED (iParRHS = 4).
! This algorithm can be considered an improved version of the direct RHS (iParRHS = 3)
!
! With iParRHS = 2, each process builds whole RHS blocks from its own Cholesky vectors and the blocks are reduced with GADGOP.
! Here the Cholesky vectors are gathered instead, and each process evaluates only the RHS columns of its own stripe,
! so the RHS is never communicated.
! This is basically similar to what the direct RHS does.
! The kernels that do the accumulation are in ADDRHS_STRIPED.

use ADDRHS_STRIPED, only: IGRP_A, IGRP_B, IGRP_C, IGRP_D1, IGRP_D2, IGRP_E, IGRP_F, IGRP_G, IGRP_H
use ADDRHS_STRIPED, only: RHSLOC_SIZES, RHSLOC_ALLOCATE, RHSLOC_LOAD, RHSLOC_FINALIZE, RHSLOC_FREE
use Symmetry_Info, only: Mul
use CHOVEC_IO, only: NPQ_CHOTYPE, NVLOC_CHOBATCH
use PrintLevel, only: USUAL, VERBOSE
use caspt2_global, only: Buff, FIMO, idxb, iParRHS, iPrGlb, PIQK
use general_data, only: NASH
use caspt2_module, only: NAES, NASHT, NBTCH, NBTCHES, NISH, NSSH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVEC

integer(kind=iwp) :: IB, IB1, IB2, IBEND, IBGRP, IBSTA, iOffi, iOffK, iOffp, iOffQ, ISYI, ISYK, ISYMT, ISYP, ISYQ, ITIER, JSYM, &
                     LBRASM, LKETSM, MINREQ, MXAVAIL, MXBGRP, MXPIQK, nAA, NBGRP, nBra, NBRABUF, NBRASM, NI, NK, nKet, NKETBUF, &
                     NKETSM, &
                     NP, NPI, NQ, NQK, NRHSGRPLOC, NRHSLOC, NRHSSYMLOC, nSh(8,3), NSYMT, NTUVX, NUMERR = 0, NV, NVLOC
! buffer sizes per JSYM, set by STRIPED_SIZES
integer(kind=iwp) :: MAXPIQK_J(8), MINPIQK_J(8), MXBATCH_J(8), MXVLOC_J(8), NAAPI_J(8), NPIBRA_J(8), NPIKET_J(8), NPIKTR_J(8), &
                     NPIPCK_J(8), NPITRA_J(8), NVECTOT_J(8)

integer(kind=iwp), allocatable :: BGRP(:,:)
real(kind=wp), allocatable :: BRA(:), CHOAA(:), KET(:), TUVX(:)

! the buffers are not used here, fixed to 1
integer(kind=iwp), parameter :: NADDBUF = 1
integer(kind=iwp), parameter :: Inactive = 1, Active = 2, Virtual = 3

  ! for test: MOLCAS_RHSTIER, see where ITIER is selected
! integer(kind=iwp) :: ISTAT, ITIERF
! character(len=8) :: TIERSTR

!                                                                      *
!***********************************************************************
!                                                                      *
nSh(1:NSYM,Inactive) = NISH(1:NSYM)
nSh(1:NSYM,Active)   = NASH(1:NSYM)
nSh(1:NSYM,Virtual)  = NSSH(1:NSYM)
!                                                                      *
!***********************************************************************
!                                                                      *

if (IPRGLB >= VERBOSE) write(u6,'(1X,A)') ' Using RHSALL2+ADDRHS (STRIPED) algorithm'

! stripe sizes for the three tiers
call RHSLOC_SIZES(NRHSLOC,NRHSGRPLOC,NRHSSYMLOC)

! ITIER (iRHSLocTier in ADDRHS_STRIPED) decides the memory usage:
!   tier 0  all blocks in memory, one write, no reads
!   tier 1  one case group, RHSLOC_LOAD reads it and writes it back
!   tier 2  one symmetry of one case group

NTUVX = NASHT**4
call STRIPED_SIZES()
MINREQ = NTUVX
do JSYM=1,NSYM
  if (NBTCH(JSYM) <= 0) cycle
  ! The same as MINSLOW in MEMORY_ESTIMATE_STRIPED:
  !   BRA, KET, CHOAA  each for the largest orbital-type pair it holds
  !   + the largest of the kernels' transposes and the gather buffer, over one batch
  !   + the block the ket transpose works on
  ! MINREQ and MINSLOW have to agree, or a tier chosen here may not fit there
  MINREQ = max(MINREQ,NTUVX+MINPIQK_J(JSYM)+2*NADDBUF+NPIKTR_J(JSYM)+ &
                      (NPIBRA_J(JSYM)+NPIKET_J(JSYM)+NAAPI_J(JSYM)+ &
                       max(2*NPITRA_J(JSYM),NPIBRA_J(JSYM),NPIKET_J(JSYM),NAAPI_J(JSYM), &
                           NPIPCK_J(JSYM)))*MXBATCH_J(JSYM))
end do
call mma_MaxDBLE(MXAVAIL)
if (NRHSLOC <= MXAVAIL-MINREQ) then
  ITIER = 0
else if (NRHSGRPLOC <= MXAVAIL-MINREQ) then
  ITIER = 1
else if (NRHSSYMLOC <= MXAVAIL-MINREQ) then
  ITIER = 2
else
  ITIER = 3
end if
call GAIGOP_SCAL(ITIER,'max')

  ! for test: force a higher tier than the one selected
! call get_environment_variable('MOLCAS_RHSTIER',TIERSTR,STATUS=ISTAT)
! if (ISTAT == 0) then
!   read(TIERSTR,*,IOSTAT=ISTAT) ITIERF
!   if ((ISTAT == 0) .and. (ITIERF > ITIER) .and. (ITIERF <= 3)) ITIER = ITIERF
! end if

if (IPRGLB > VERBOSE) then
  write(u6,*)
  write(u6,'(A)') '  Memory for the striped RHS'
  write(u6,'(A,2X,I16)') '   allocatable:    ',MXAVAIL
  write(u6,'(A,2X,I16)') '   rest, at least: ',MINREQ
  write(u6,'(A,2X,I16)') '   whole local RHS:',NRHSLOC
  write(u6,'(A,2X,I16)') '   largest group:  ',NRHSGRPLOC
  write(u6,'(A,2X,I16)') '   one symmetry:   ',NRHSSYMLOC
  write(u6,'(A,2X,I16)') '   selected tier:  ',ITIER
  write(u6,*)
end if

NSYMT = 1
if (ITIER == 1) then
  NRHSLOC = NRHSGRPLOC
  if (IPRGLB >= VERBOSE) write(u6,'(1X,A)') ' The striped RHS does not fit. One case group is held at a time'
else if (ITIER == 2) then
  NRHSLOC = NRHSSYMLOC
  NSYMT = NSYM
  if (IPRGLB >= VERBOSE) write(u6,'(1X,A)') ' Not even a whole case group fits. One symmetry of one is held at a time'
end if
if (ITIER == 3) then
  if (IPRGLB >= USUAL) then
    write(u6,*)
    write(u6,'(1X,A)') ' Not enough memory for the striped RHS. The replicated algorithm is used instead'
  end if
  if (IPRGLB >= VERBOSE) then
    write(u6,'(2X,A8,2X,I14)') 'MXAVAIL ',MXAVAIL
    write(u6,'(2X,A8,2X,I14)') 'NRHSSYM ',NRHSSYMLOC
    write(u6,'(2X,A8,2X,I14)') 'MINREQ  ',MINREQ
  end if
  if (IPRGLB >= USUAL) write(u6,*)
  iParRHS = 2
  call RHS_ZERO(IVEC)
  call RHSALL2(IVEC)
  return
end if

call RHSLOC_ALLOCATE(ITIER,NRHSLOC)

! TUVX RHSX        Na^4
! TJVX RHSA        Na^3 Ni
! TJVL RHSB        Na^2 Ni^2
! AJVX RHSD1       Na^2 Ni   Ns
! AJCL RHSH             Ni^2 Ns^2     N^4
! AUVX RHSC        Na^3      Ns
! AUCX RHSF        Na^2      Ns^2
! AUVL RHSD2       Na^2 Ni   Ns
! AUCL RHSG        Na   Ni   Ns^2     N^3
! AJVL RHSE        Na   Ni^2 Ns       N^3

!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate and clear TUVX, two-electron integrals for active
! orbital indices only: Simple storage, same as for GAMMA2.
! TUVX is kept allocated until end of subroutine. NTUVX is set above, it is part of the memory estimate.
call mma_allocate(TUVX,NTUVX,Label='TUVX')
TUVX(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
do JSYM=1,NSYM

  IB1 = NBTCHES(JSYM)+1
  IB2 = NBTCHES(JSYM)+NBTCH(JSYM)

  MXBGRP = IB2-IB1+1
  if (MXBGRP <= 0) cycle
  call mma_allocate(BGRP,2,MXBGRP,Label='BGRP')
  BGRP(1,1:IB2-IB1+1) = [(IB,IB=IB1,IB2)]
  BGRP(2,1:IB2-IB1+1) = [(IB,IB=IB1,IB2)]
  NBGRP = MXBGRP

  call MEMORY_ESTIMATE_STRIPED()
  if (IPRGLB > VERBOSE) then
    write(u6,*)
    write(u6,'(A,I12)') '  Number of Cholesky batches: ',IB2-IB1+1
    write(u6,'(A,I12)') '  Number of batch groups:     ',NBGRP
    write(u6,*)
  end if

  ! buffers are kept allocated until the end of JSYM loop.
  call mma_allocate(PIQK,MXPIQK,Label='PIQK')
  call mma_allocate(BUFF,NADDBUF,Label='BUFF')
  call mma_allocate(IDXB,NADDBUF,Label='IDXB')

  call mma_allocate(BRA,NBRABUF,Label='BRA')
  call mma_allocate(KET,NKETBUF,Label='KET')

  ! Loop over groups of batches of Cholesky vectors

  do IBGRP=1,NBGRP
    IBSTA = BGRP(1,IBGRP)
    IBEND = BGRP(2,IBGRP)

    ! the buffers hold the vectors of every process, NV is the global count
    NVLOC = sum(NVLOC_CHOBATCH(IBSTA:IBEND))
    NV = NVLOC
    call GAIGOP_SCAL(NV,'+')

    if (IPRGLB > VERBOSE) then
      write(u6,'(A,I12)') '  Cholesky vectors in this group = ',NV
      write(u6,*)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read the L(VX) vectors (all symmetries) into their own buffer
    ! cases A, D1 and C use them as kets

    call mma_allocate(CHOAA,max(NAAPI_J(JSYM)*NV,1),Label='CHOAA')
    call Get_Cholesky_Vectors(Active,Active,JSYM,CHOAA,size(CHOAA),nAA,IBSTA,IBEND)
    call Gather_Cholesky_Vectors(Active,Active,JSYM,CHOAA,size(CHOAA),nAA,NVLOC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Assemble contributions to TUVX integrals
    ! Reuse the L(VX) vectors as L(TU) bra vectors

    LBRASM = 1
    do ISYI=1,NSYM
      NI = NASH(ISYI)
      iOffi = NAES(iSYI)
      if (NI == 0) cycle
      ISYP = Mul(ISYI,JSYM)
      NP = NASH(ISYP)
      iOffp = NAES(iSYP)
      if (NP == 0) cycle
      NPI = NP*NI
      NBRASM = NPI*NV
      LKETSM = 1

      do ISYK=1,NSYM
        NK = NASH(ISYK)
        iOffK = NAES(iSYK)
        if (NK == 0) cycle
        ISYQ = Mul(ISYK,JSYM)
        NQ = NASH(ISYQ)
        iOffQ = NAES(iSYQ)
        if (NQ == 0) cycle
        NQK = NQ*NK
        NKETSM = NQK*NV

        if (NPI*NQK > mxPIQK) then
          write(u6,*) 'NPIQK larger than mxPIQK in TUVX, bug?'
          call AbEnd()
        end if
        call DGEMM_('N','T',NPI,NQK,NV,One,CHOAA(LBRASM),NPI,CHOAA(LKETSM),NQK,Zero,PIQK,NPI)

        call ADDTUVX(NP,NI,NQ,NK,NASHT,iOffP,iOffI,iOffQ,iOffK,TUVX,nTUVX,PIQK,NPI*NQK,NUMERR)

        LKETSM = LKETSM+NKETSM
      end do
      LBRASM = LBRASM+NBRASM
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read bra (Cholesky vectors) in the form L(TJ): All symmetries

    call Get_Cholesky_Vectors(Inactive,Active,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    call Gather_Cholesky_Vectors(Inactive,Active,JSYM,BRA,size(BRA),nBra,NVLOC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Assemble contributions to TJVX
    ! Loop over the bras and kets, form <A|0>

    do ISYMT=1,NSYMT
      call RHSLOC_LOAD(IVEC,IGRP_A,ISYMT)
      call Process_RHS_Block(Inactive,Active,Active,Active,'A ',BRA,nBra,CHOAA,nAA,nSh,JSYM,IVEC,NV)
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! TJVL RHSB
    ! TJVL: Use TJ buffer as if it was VL, form <B|0>

    do ISYMT=1,NSYMT
      call RHSLOC_LOAD(IVEC,IGRP_B,ISYMT)
      call Process_RHS_Block(Inactive,Active,Inactive,Active,'B ',BRA,nBra,BRA,nBra,nSh,JSYM,IVEC,NV)
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read the L(AJ) vectors into KET, the L(TJ) vectors stay in BRA

    call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,KET,size(KET),nKet,IBSTA,IBEND)
    call Gather_Cholesky_Vectors(Inactive,Virtual,JSYM,KET,size(KET),nKet,NVLOC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AJVL RHSE
    ! AJVL: L(AJ) are the bras here, L(VL) the kets. Form <E|0>

    do ISYMT=1,NSYMT
      call RHSLOC_LOAD(IVEC,IGRP_E,ISYMT)
      call Process_RHS_Block(Inactive,Virtual,Inactive,Active,'E ',KET,nKet,BRA,nBra,nSh,JSYM,IVEC,NV)
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AJCL RHSH
    ! AJCL: Use the L(AJ) vectors as if they were CL, form <H|0>

    do ISYMT=1,NSYMT
      call RHSLOC_LOAD(IVEC,IGRP_H,ISYMT)
      call Process_RHS_Block(Inactive,Virtual,Inactive,Virtual,'H ',KET,nKet,KET,nKet,nSh,JSYM,IVEC,NV)
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AJVX RHSD1
    ! Loop over the bra and ket vectors, form <D1|0>

    do ISYMT=1,NSYMT
      call RHSLOC_LOAD(IVEC,IGRP_D1,ISYMT)
      call Process_RHS_Block(Inactive,Virtual,Active,Active,'D1',KET,nKet,CHOAA,nAA,nSh,JSYM,IVEC,NV)
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read the L(AU) vectors into BRA, the L(AJ) vectors stay in KET

    call Get_Cholesky_Vectors(Active,Virtual,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    call Gather_Cholesky_Vectors(Active,Virtual,JSYM,BRA,size(BRA),nBra,NVLOC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUCL RHSG
    ! Use the L(AJ) vectors as if they were CL, form <G|0>
    ! ADDRHSG_STRIPED wants them as (l,c,P), so transpos (c,l,P) to (l,c,P) in place
    ! Here is the last place L(AJ) is used

    call TRANSPOSE_KET()

    do ISYMT=1,NSYMT
      call RHSLOC_LOAD(IVEC,IGRP_G,ISYMT)
      call Process_RHS_Block(Active,Virtual,Inactive,Virtual,'G ',BRA,nBra,KET,nKet,nSh,JSYM,IVEC,NV)
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUCX RHSF
    ! AUCX: Use AU buffer still in core as if it was CX, form <F|0>

    do ISYMT=1,NSYMT
      call RHSLOC_LOAD(IVEC,IGRP_F,ISYMT)
      call Process_RHS_Block(Active,Virtual,Active,Virtual,'F ',BRA,nBra,BRA,nBra,nSh,JSYM,IVEC,NV)
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUVX RHSC
    ! AUVX: Loop over the bras and kets, form <C|0>

    do ISYMT=1,NSYMT
      call RHSLOC_LOAD(IVEC,IGRP_C,ISYMT)
      call Process_RHS_Block(Active,Virtual,Active,Active,'C ',BRA,nBra,CHOAA,nAA,nSh,JSYM,IVEC,NV)
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read the L(VL) vectors again, replacing L(AJ).

    call Get_Cholesky_Vectors(Inactive,Active,JSYM,KET,size(KET),nKet,IBSTA,IBEND)
    call Gather_Cholesky_Vectors(Inactive,Active,JSYM,KET,size(KET),nKet,NVLOC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUVL RHSD2
    ! Loop over bras and kets, form <D2|0>.

    do ISYMT=1,NSYMT
      call RHSLOC_LOAD(IVEC,IGRP_D2,ISYMT)
      call Process_RHS_Block(Active,Virtual,Inactive,Active,'D2',BRA,nBra,KET,nKet,nSh,JSYM,IVEC,NV)
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    call mma_deallocate(CHOAA)
    ! End of loop over batch groups, IBGRP
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call mma_deallocate(BRA)
  call mma_deallocate(KET)
  call mma_deallocate(PIQK)
  call mma_deallocate(BUFF)
  call mma_deallocate(IDXB)
  call mma_deallocate(BGRP)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! End of loop over JSYM
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call RHSLOC_FINALIZE(IVEC)
call RHSLOC_FREE()

! The RHS elements of Cases A, C, D1  need a correction:
call MODRHS(IVEC,FIMO,size(FIMO))

! Put TUVX on disk for possible later use:
call PT2_PUT(NTUVX,'TUVX',TUVX)
call mma_deallocate(TUVX)
!                                                                      *
!***********************************************************************
!                                                                      *

contains

subroutine TRANSPOSE_KET()

  use caspt2_module, only: NISH, NSSH, NSYM
  use Symmetry_Info, only: Mul

  integer(kind=iwp) :: IC, IL, IOFF, IP, ISYC, ISYL, LKET, NC, NL, NTMP

  real(kind=wp), allocatable :: TMP(:)

  ! Transpose KET from (c,l,P) to (l,c,P) in place
  ! case G needs no transposed copy of them

  NTMP = 0
  do ISYL=1,NSYM
    ISYC = Mul(ISYL,JSYM)
    NTMP = max(NTMP,NSSH(ISYC)*NISH(ISYL))
  end do
  if (NTMP == 0) return
  call mma_allocate(TMP,NTMP,Label='ChoKetT')

  LKET = 1
  do ISYL=1,NSYM
    NL = NISH(ISYL)
    if (NL == 0) cycle
    ISYC = Mul(ISYL,JSYM)
    NC = NSSH(ISYC)
    if (NC == 0) cycle
    do IP=1,NV
      IOFF = LKET-1+NC*NL*(IP-1)
      TMP(1:NC*NL) = KET(IOFF+1:IOFF+NC*NL)
      do IC=1,NC
        do IL=1,NL
          KET(IOFF+IL+NL*(IC-1)) = TMP(IC+NC*(IL-1))
        end do
      end do
    end do
    LKET = LKET+NC*NL*NV
  end do

  call mma_deallocate(TMP)

end subroutine TRANSPOSE_KET

!-----------------------------------------------------------------------

subroutine STRIPED_SIZES()

  integer(kind=iwp) :: ICASE, ISYI, ISYK, ISYP, ISYQ, JB1, JB2, JSY, MAXPIQK, MINPIQK, MXVLOC, NI, NK, NP, NPI, NPI_AV, NPI_IA, &
                       NPI_IV, NPIKTR, NPIPCK, NPITRA, NPMAX, NPMIN, NQ, NQK

  integer(kind=iwp), allocatable :: NVEFF(:)

  ! orbital types (bra p, bra i, ket q, ket k) of each case,
  ! icase = 1, 2, 3, 4, 5, 6, 7, 8, 9 for A, B, D1, H, C, F, D2, G, E
  integer(kind=iwp), parameter :: ITYPE(4,9) = reshape([Inactive,  Active,  Active,  Active, &
                                                        Inactive,  Active,Inactive,  Active, &
                                                        Inactive, Virtual,  Active,  Active, &
                                                        Inactive, Virtual,Inactive, Virtual, &
                                                          Active, Virtual,  Active,  Active, &
                                                          Active, Virtual,  Active, Virtual, &
                                                          Active, Virtual,Inactive,  Active, &
                                                          Active, Virtual,Inactive, Virtual, &
                                                        Inactive, Virtual,Inactive,  Active],[4,9])

  ! Buffer sizes for all JSYM at once
  ! needed before RHSLoc is allocated, the tier cannot be changed afterwards

  NPIBRA_J(:) = 0
  NPIKET_J(:) = 0
  NPIKTR_J(:) = 0
  NPIPCK_J(:) = 0
  NPITRA_J(:) = 0
  MXVLOC_J(:) = 0
  NAAPI_J(:) = 0
  MINPIQK_J(:) = 0
  MAXPIQK_J(:) = 0
  MXBATCH_J(:) = 0
  NVECTOT_J(:) = 0

  do JSY=1,NSYM

    JB1 = NBTCHES(JSY)+1
    JB2 = NBTCHES(JSY)+NBTCH(JSY)
    if (JB2 < JB1) cycle

    ! Pairs per vector, one count per orbital-type pair:
    !   NPQ_CHOTYPE(1,..) = (inactive,active)  L(TJ) and L(VL)
    !   NPQ_CHOTYPE(2,..) = (active,active)    L(VX)
    !   NPQ_CHOTYPE(3,..) = (active,secondary) L(AU)
    !   NPQ_CHOTYPE(4,..) = (inactive,secondary) L(AJ)
    !   BRA, KET  whole buffers, summed over the symmetry blocks, each for the largest type it holds
    !   CHOBT     the bra transpose of ADDRHSF/G_STRIPED, one symmetry block at a time, so a maximum
    NPI_IA = 0
    NPI_AV = 0
    NPI_IV = 0
    NPIPCK = 0
    NPITRA = 0
    NPIKTR = 0
    do ISYQ=1,NSYM
      NPI_IA = NPI_IA+NPQ_CHOTYPE(1,ISYQ,JSY)
      NAAPI_J(JSY) = NAAPI_J(JSY)+NPQ_CHOTYPE(2,ISYQ,JSY)
      NPI_AV = NPI_AV+NPQ_CHOTYPE(3,ISYQ,JSY)
      NPI_IV = NPI_IV+NPQ_CHOTYPE(4,ISYQ,JSY)
      NPIPCK = max(NPIPCK,NISH(ISYQ),NASH(ISYQ))
      NPITRA = max(NPITRA,NPQ_CHOTYPE(3,ISYQ,JSY))
      NPIKTR = max(NPIKTR,NPQ_CHOTYPE(4,ISYQ,JSY))
    end do
    NPIBRA_J(JSY) = max(NPI_IA,NPI_AV) ! L(TJ) until case C, then L(AU)
    NPIKET_J(JSY) = max(NPI_IV,NPI_IA) ! L(AJ) until case C, then L(VL)
    NPIPCK_J(JSY) = NPIPCK             ! CHOBA, the one-orbital slice C, D1 and D2 pack
    NPITRA_J(JSY) = NPITRA             ! CHOBT, and CHOKC beside it in ADDRHSF_STRIPED_O
    NPIKTR_J(JSY) = NPIKTR             ! the single block TRANSPOSE_KET works on

    ! Size needed to hold the integral matrix, at least the TUVX integrals (NASHT**4)
    ! The whole (pi,qk) block is never formed, the scratch holds one chunk of it
    ! NPMAX: the size at which the outer loop needs no chunking
    ! NPMIN: the smallest chunk that still works
    MAXPIQK = NASHT**4
    MINPIQK = NASHT**4
    do ICASE=1,9
      do ISYI=1,NSYM
        NI = nSh(ISYI,ITYPE(1,ICASE))
        ISYP = Mul(ISYI,JSY)
        NP = nSh(ISYP,ITYPE(2,ICASE))
        NPI = NP*NI
        do ISYK=1,NSYM
          NK = nSh(ISYK,ITYPE(3,ICASE))
          ISYQ = Mul(ISYK,JSY)
          NQ = nSh(ISYQ,ITYPE(4,ICASE))
          NQK = NQ*NK
          select case (ICASE)
            case (1)                                 ! A: Y((t,j),(v,x))
              NPMAX = NPI*NQK
              NPMIN = NP*NQK
            case (2,4)                               ! B, H: Y((.,l),.)
              NPMAX = NP*NQK
              NPMIN = NP*NQ
            case (3,5,7)                             ! D1, C, D2: one a/j
              NPMAX = NI*NQK
              NPMIN = NI*NQK
            case (6,8)                               ! F, G: two blocks
              ! either NP or NQ can be the un-batched one, depending on ISYA vs ISYC
              NPMAX = max(2*NI*NK*NP,NI*NK*NQ)
              NPMIN = max(2*NI*NK,NI*NK*NQ,NI*NK*NP)
            case default                             ! 9 = E: Y(v,a)
              ! the diagonal kernel takes one l at a time, the off-diagonal one the whole range
              NPMAX = NP*max(2*NQ,NQ*NK)
              NPMIN = max(2*NQ,NQ*NK)
          end select
          MAXPIQK = max(MAXPIQK,NPMAX)
          MINPIQK = max(MINPIQK,NPMIN)
        end do
      end do
    end do
    MAXPIQK_J(JSY) = MAXPIQK
    MINPIQK_J(JSY) = MINPIQK

    ! Total number of Cholesky vectors.
    call mma_allocate(NVEFF,[JB1,JB2],Label='NVEFF')
    NVEFF(JB1:JB2) = NVLOC_CHOBATCH(JB1:JB2)
    call GAIGOP(NVEFF,JB2-JB1+1,'+')
    NVECTOT_J(JSY) = sum(NVEFF(JB1:JB2))
    MXBATCH_J(JSY) = maxval(NVEFF(JB1:JB2))
    call mma_deallocate(NVEFF)

    ! the largest over the processes, every process has to reach the same tier
    MXVLOC = sum(NVLOC_CHOBATCH(JB1:JB2))
    call GAIGOP_SCAL(MXVLOC,'max')
    MXVLOC_J(JSY) = MXVLOC
  end do

end subroutine STRIPED_SIZES

!-----------------------------------------------------------------------

subroutine MEMORY_ESTIMATE_STRIPED()

  integer(kind=iwp) :: IB, IBGRP, IBRANCH, MAXPIQK, MINGOOD, MINNICE, MINPIQK, MINSLOW, MXBATCH, MXCHOVEC, MXRHS, MXVLOC, NAABUF, &
                       NAAPI, NCHOVEC, NCHUNK, NPIBRA, NPIGAT, NPIKET, NPIKTR, NPIPCK, NPIPER, NPITRA, NPIXTR, NV, NVECTOT, &
                       NVGRP, NXTRA

  integer(kind=iwp), allocatable :: NVEFF(:)

  ! for test: MOLCAS_RHSBRANCH, see where IBRANCH is set
! integer(kind=iwp) :: ISTAT
! character(len=8) :: BRSTR

  ! Striped counterpart of MEMORY_ESTIMATE

  ! No reservation is needed for the RHS
  ! RHSLoc is allocated before this is reached, mma_MaxDBLE accounts for it
  MXRHS = 0

  NPIBRA = NPIBRA_J(JSYM)
  NPIKET = NPIKET_J(JSYM)
  NPIKTR = NPIKTR_J(JSYM)
  NPIPCK = NPIPCK_J(JSYM)
  NPITRA = NPITRA_J(JSYM)
  MXVLOC = MXVLOC_J(JSYM)
  NAAPI = NAAPI_J(JSYM)
  MAXPIQK = MAXPIQK_J(JSYM)
  MINPIQK = MINPIQK_J(JSYM)
  NVECTOT = NVECTOT_J(JSYM)
  MXBATCH = MXBATCH_J(JSYM)

  ! NVEFF is needed again below to cut the batch groups
  call mma_allocate(NVEFF,[IB1,IB2],Label='NVEFF')
  NVEFF(IB1:IB2) = NVLOC_CHOBATCH(IB1:IB2)
  call GAIGOP(NVEFF,IB2-IB1+1,'+')

  ! Cost per Cholesky vector
  NPIPER = NPIBRA+NPIKET+NAAPI         ! BRA, KET and CHOAA
  NPIGAT = max(NPIBRA,NPIKET,NAAPI)    ! LocBuf in Gather_Cholesky_Vectors, made for every type
  NPIXTR = max(2*NPITRA,NPIGAT,NPIPCK) ! that, CHOBT and the CHOKC beside it, or CHOBA

  ! The branch below and the group cuts have to be the same everywhere, so set the common MXAVAIL
  call mma_MaxDBLE(MXAVAIL)
  call GAIGOP_SCAL(MXAVAIL,'min')

  NXTRA = max(2*NPITRA*NVECTOT,NPIPCK*NVECTOT,NPIGAT*MXVLOC,NPIKTR)
  MINNICE = MXRHS+MAXPIQK+2*NADDBUF+NPIPER*NVECTOT+NXTRA
  MINGOOD = MXRHS+MINPIQK+2*NADDBUF+NPIPER*NVECTOT+NXTRA
  ! MINSLOW can exceed MINGOOD when there is only one batch, clamp it
  MINSLOW = min(MXRHS+MINPIQK+2*NADDBUF+(NPIPER+NPIXTR)*MXBATCH+NPIKTR,MINGOOD)

  IBRANCH = 0
  ! for test: force the convenient (2) or the minimum (3) branch
! call get_environment_variable('MOLCAS_RHSBRANCH',BRSTR,STATUS=ISTAT)
! if (ISTAT == 0) then
!   read(BRSTR,*,IOSTAT=ISTAT) IBRANCH
!   if (ISTAT /= 0) IBRANCH = 0
!   if (IBRANCH == 2) MXAVAIL = min(MXAVAIL,MINGOOD)
!   if (IBRANCH == 3) MXAVAIL = min(MXAVAIL,MINSLOW)
! end if

  if (IPRGLB > VERBOSE) then
    write(u6,*)
    write(u6,'(A,I1)') '  Memory estimates in RHSALL (striped), SYM ',JSYM
    write(u6,'(A,2X,I16)') '   allocatable:    ',MXAVAIL
    write(u6,'(A,2X,I16)') '   recommended:    ',MINNICE
    write(u6,'(A,2X,I16)') '   convenient:     ',MINGOOD
    write(u6,'(A,2X,I16)') '   minimum:        ',MINSLOW
    write(u6,*)
    if ((MXAVAIL >= MINNICE) .and. (IBRANCH == 0)) then
      write(u6,*) ' I can use all cholesky vectors at once'
      write(u6,*) ' as well as the whole integral matrix.'
    else if ((MXAVAIL >= MINGOOD) .and. (IBRANCH /= 3)) then
      write(u6,*) ' I will group batches of cholesky vectors'
      write(u6,*) ' and then maximize use of the integral matrix.'
    else if (MXAVAIL >= MINSLOW) then
      write(u6,*) ' I will at least try to group batches.'
    else
      write(u6,*) ' Do you see my problem?'
    end if
  end if

  NVGRP = 0
  if ((MXAVAIL >= MINNICE) .and. (IBRANCH == 0)) then
    ! group all batches and take the maximum needed for integrals
    NVGRP = NVECTOT
    NBGRP = 1
    BGRP(1,1) = IB1
    BGRP(2,1) = IB2
    MXPIQK = MAXPIQK
  else if ((MXAVAIL >= MINGOOD) .and. (IBRANCH /= 3)) then
    ! group all batches and try to max out the integrals, keeping them larger than the minimum
    NVGRP = NVECTOT
    NBGRP = 1
    BGRP(1,1) = IB1
    BGRP(2,1) = IB2
    NCHUNK = (MXAVAIL-MXRHS-2*NADDBUF-NPIPER*NVECTOT-NXTRA)/MINPIQK
    MXPIQK = MINPIQK*NCHUNK
  else if (MXAVAIL >= MINSLOW) then
    ! As many vectors per group as fit, never fewer than the largest batch
    ! MXAVAIL >= MINSLOW already implies MXCHOVEC >= MXBATCH, the max() is redundant
    MXPIQK = MINPIQK
    MXCHOVEC = max((MXAVAIL-MXRHS-MXPIQK-2*NADDBUF-NPIKTR)/(NPIPER+NPIXTR),MXBATCH)
    NVGRP = MXCHOVEC
    ! create batch groups that have at most MXCHOVEC cholesky vectors
    NCHOVEC = 0
    IBGRP = 1
    BGRP(1,IBGRP) = IB1
    do IB=IB1,IB2
      NV = NVEFF(IB)
      NCHOVEC = NCHOVEC+NV
      if (NCHOVEC > MXCHOVEC) then
        BGRP(2,IBGRP) = IB-1
        IBGRP = IBGRP+1
        BGRP(1,IBGRP) = IB
        NCHOVEC = NV
      end if
    end do
    BGRP(2,IBGRP) = IB2
    NBGRP = IBGRP
  else
    write(u6,*)
    write(u6,*) '  Not enough memory in RHSALL2_STRIPED...'
    write(u6,'(A,I16)') '   allocatable:    ',MXAVAIL
    write(u6,'(A,I16)') '   minimum:        ',MINSLOW
    call AbEnd()
  end if

  ! sized from the group that was actually chosen
  NBRABUF = NPIBRA*NVGRP
  NKETBUF = NPIKET*NVGRP
  NAABUF = NAAPI*NVGRP
  NXTRA = max(2*NPITRA*NVGRP,NPIPCK*NVGRP,NPIGAT*min(MXVLOC,NVGRP),NPIKTR)

  ! sanity check, should not happen.
  if (MXRHS > MXAVAIL-NBRABUF-NKETBUF-NAABUF-NXTRA-MXPIQK-2*NADDBUF) then
    write(u6,*)
    write(u6,*) 'RHSALL2_STRIPED: RHS allocation will starve.'
    write(u6,*) 'Possible bug in memory estimate.'
    write(u6,*) 'This should not happen, please report.'
    write(u6,*)
    write(u6,'(2X,A8,2X,I14)') 'MXAVAIL ',MXAVAIL
    write(u6,'(2X,A8,2X,I14)') 'MXRHS   ',MXRHS
    write(u6,'(2X,A8,2X,I14)') 'NBRABUF ',NBRABUF
    write(u6,'(2X,A8,2X,I14)') 'NKETBUF ',NKETBUF
    write(u6,'(2X,A8,2X,I14)') 'NPIQK   ',MXPIQK
    write(u6,'(2X,A8,2X,I14)') 'NADDBUF  ',NADDBUF
    call AbEnd()
  end if

  if (IPRGLB > VERBOSE) then
    write(u6,*)
    write(u6,'(A16,2A16)') '  Buffer sizes:','           used','          ideal'
    write(u6,'(A16,2I16)') '  ChoVecs:  ',NBRABUF+NKETBUF+NAABUF+NXTRA, &
                                          NPIPER*NVECTOT+max(2*NPITRA*NVECTOT,NPIPCK*NVECTOT,NPIGAT*MXVLOC,NPIKTR)
    write(u6,'(A16,2I16)') '  Integral: ',MXPIQK,MAXPIQK
    write(u6,*)
  end if

  call mma_deallocate(NVEFF)

end subroutine MEMORY_ESTIMATE_STRIPED

end subroutine RHSALL2_STRIPED

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(RHSALL2_STRIPED)

#endif
