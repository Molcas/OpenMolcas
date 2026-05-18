!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine RHSALL2(IVEC)
! ----------------------------------------------------------------
! Code for processing all the cholesky vectors
! in construction of caspt2 right-hand-side array
! Also form the active two-electron integrals 'TUVX'.
! ================================================================

use Symmetry_Info, only: Mul
use CHOVEC_IO, only: NVLOC_CHOBATCH
use PrintLevel, only: VERBOSE
use caspt2_global, only: Buff, FIMO, idxb, iPrGlb, PIQK
use caspt2_module, only: NAES, NASH, NASHT, NBTCH, NBTCHES, NISH, NSSH, NSYM
#ifdef _DEBUGPRINT_
use caspt2_module, only: NASUP, NISUP
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IVEC
integer(kind=iwp) :: IB, IB1, IB2, IBEND, IBGRP, IBSTA, iOffi, iOffK, iOffp, iOffQ, ISYI, ISYK, ISYP, ISYQ, JSYM, LBRASM, LKETSM, &
                     MXBGRP, MXPIQK, NADDBUF, NBGRP, nBra, NBRASM, NCHOBUF, NG1, NG2, NI, NK, nKet, NKETSM, NP, NPI, NQ, NQK, &
                     nSh(8,3), NTUVX, NUMERR = 0, NV
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ICASE, ISYM, lg_W, NAS, NIS
real(kind=wp) :: DNRM2
real(kind=wp), external :: RHS_DDot
#endif
integer(kind=iwp), allocatable :: BGRP(:,:)
real(kind=wp), allocatable :: BRA(:), KET(:), TUVX(:)
integer(kind=iwp), parameter :: Inactive = 1, Active = 2, Virtual = 3

!                                                                      *
!***********************************************************************
!                                                                      *
call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
call ICopy(NSYM,NASH,1,nSh(1,Active),1)
call ICopy(NSYM,NSSH,1,nSh(1,Virtual),1)
!                                                                      *
!***********************************************************************
!                                                                      *

if (IPRGLB >= VERBOSE) write(u6,'(1X,A)') ' Using RHSALL2+ADDRHS algorithm'

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

!     Initialize RHS array as 'vector' nr IVEC=IRHS on LUSOLV:

!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate and clear TUVX, two-electron integrals for active
! orbital indices only: Simple storage, same as for GAMMA2.
! TUVX is kept allocated until end of subroutine.
NG1 = NASHT**2
NG2 = NG1**2
NTUVX = NG2
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
  IBGRP = 1
  do IB=IB1,IB2
    BGRP(1,IBGRP) = IB
    BGRP(2,IBGRP) = IB
    IBGRP = IBGRP+1
  end do
  NBGRP = MXBGRP

  call MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,NCHOBUF,MXPIQK,NADDBUF)
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

  call mma_allocate(BRA,NCHOBUF,Label='BRA')
  call mma_allocate(KET,NCHOBUF,Label='KET')

  ! Loop over groups of batches of Cholesky vectors

  !IBSTEP = 1

  !do IBSTA=IB1,IB2,IBSTEP
  do IBGRP=1,NBGRP
    IBSTA = BGRP(1,IBGRP)
    IBEND = BGRP(2,IBGRP)

    NV = 0
    do IB=IBSTA,IBEND
      NV = NV+NVLOC_CHOBATCH(IB)
    end do

    if (IPRGLB > VERBOSE) then
      write(u6,'(A,I12)') '  Cholesky vectors in this group = ',NV
      write(u6,*)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read kets (Cholesky vectors) in the form L(VX), all symmetries:

    call Get_Cholesky_Vectors(Active,Active,JSYM,KET,size(KET),nKet,IBSTA,IBEND)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Assemble contributions to TUVX integrals
    ! Reuse the ket vectors as L(TU) bra vectors

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
        call DGEMM_('N','T',NPI,NQK,NV,One,KET(LBRASM),NPI,KET(LKETSM),NQK,Zero,PIQK,NPI)

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
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Assemble contributions to TJVX
    ! Loop over the bras and kets, form <A|0>

    call Process_RHS_Block(Inactive,Active,Active,Active,'A ',BRA,nBra,KET,nKet,nSh,JSYM,IVEC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! TJVL RHSB
    ! TJVL: Use TJ buffer as if it was VL, form <B|0>

    call Process_RHS_Block(Inactive,Active,Inactive,Active,'B ',BRA,nBra,BRA,nBra,nSh,JSYM,IVEC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read bra (Cholesky vectors) in the form L(AJ), form <D1|0>
    ! We still have L(VX) vectors in core, at KET.

    call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AJVX RHSD1
    ! Loop over the bra and ket vectors.

    call Process_RHS_Block(Inactive,Virtual,Active,Active,'D1',BRA,nBra,KET,nKet,nSh,JSYM,IVEC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AJCL RHSH
    ! AJCL: Use AJ buffer still in core as if it was CL, form <H|0>

    call Process_RHS_Block(Inactive,Virtual,Inactive,Virtual,'H ',BRA,nBra,BRA,nBra,nSh,JSYM,IVEC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read Bra (Cholesky vectors)= L(AU)

    call Get_Cholesky_Vectors(Active,Virtual,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUVX RHSC
    ! AUVX: Loop over the bras and kets

    call Process_RHS_Block(Active,Virtual,Active,Active,'C ',BRA,nBra,KET,nKet,nSh,JSYM,IVEC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUCX RHSF
    ! AUCX: Use AU buffer still in core as if it was CX, form <F|0>

    call Process_RHS_Block(Active,Virtual,Active,Virtual,'F ',BRA,nBra,BRA,nBra,nSh,JSYM,IVEC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read kets (Cholesky vectors) in the form L(VL), all symmetries:

    call Get_Cholesky_Vectors(Inactive,Active,JSYM,KET,size(KET),nKet,IBSTA,IBEND)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUVL RHSD2
    ! Loop over bras and kets, form <D2|0>.

    call Process_RHS_Block(Active,Virtual,Inactive,Active,'D2',BRA,nBra,KET,nKet,nSh,JSYM,IVEC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read kets (Cholesky vectors) in the form L(CL), all symmetries:

    call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,KET,size(KET),nKet,IBSTA,IBEND)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUCL RHSG
    ! Loop over bras and kets, form  <G|0>

    call Process_RHS_Block(Active,Virtual,Inactive,Virtual,'G ',BRA,nBra,KET,nKet,nSh,JSYM,IVEC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read bra vectors AJ

    call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read kets in the form L(VL)

    call Get_Cholesky_Vectors(Inactive,Active,JSYM,KET,size(KET),nKet,IBSTA,IBEND)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AJVL RHSE
    ! AJVL: Loop over bras and kets. Form <E|0>

    call Process_RHS_Block(Inactive,Virtual,Inactive,Active,'E ',BRA,nBra,KET,nKet,nSh,JSYM,IVEC,NV)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! End of loop over batches, IB
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
! Synchronized add RHS partial arrays from all nodes into each node.

!-SVC: read the DRA's from disk and copy them all to LUSOLV to continue
!      in serial mode.  FIXME: this call has to be removed when we reach
!      full parallel capabilities
!call SYNRHS(IVEC)
!-SVC: at this point, the RHS elements are on disk, both in LUSOLV and
!      as DRAs with the name RHS_XX_XX_XX with XX a number representing
!      the case, symmetry, and rhs vector respectively.

! The RHS elements of Cases A, C, D1  need a correction:
call MODRHS(IVEC,FIMO,size(FIMO))

#ifdef _DEBUGPRINT_
! compute and print RHS fingerprints
write(u6,'(1X,A4,1X,A3,1X,A18)') 'Case','Sym','Fingerprint'
write(u6,'(1X,A4,1X,A3,1X,A18)') '====','===','==========='
do ICASE=1,13
  do ISYM=1,NSYM
    NAS = NASUP(ISYM,ICASE)
    NIS = NISUP(ISYM,ICASE)
    if (NAS*NIS /= 0) then
      call RHS_ALLO(NAS,NIS,lg_W)
      call RHS_READ(NAS,NIS,lg_W,iCASE,iSYM,iVEC)
      DNRM2 = RHS_DDOT(NAS,NIS,lg_W,lg_W)
      write(u6,'(1X,I4,1X,I3,1X,F18.11)') ICASE,ISYM,DNRM2
    end if
  end do
end do
#endif

! Synchronized add tuvx partial arrays from all nodes into each node.
call CHO_GADGOP(TUVX,NTUVX,'+')
! Put TUVX on disk for possible later use:
call PT2_PUT(NTUVX,'TUVX',TUVX)
call mma_deallocate(TUVX)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine RHSALL2
