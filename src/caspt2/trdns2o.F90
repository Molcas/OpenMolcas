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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine TRDNS2O(IVEC,JVEC,DPT2,MDPT2,NDPT2,SCAL)
! Add to the block-off-diag parts of transition density matrix,
!    DPT2(p,q) = Add <IVEC| E(p,q) |JVEC>.
! i.e. inact/act, inact/virt and act/virt submatrices only,
! where IVEC, JVEC stands for the 1st-order perturbed CASPT2
! wave functions stored as vectors nr IVEC, JVEC on LUSOLV.

#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use EQSOLV, only: IfCoup
use fake_GA, only: GA_Arrays
use caspt2_global, only: LISTS
use caspt2_module, only: FockType, G1SecIn, nActEl, nAsh, nASup, nInDep, nIsh, nISup, nOrb, nSsh, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IVEC, JVEC, MDPT2
real(kind=wp), intent(inout) :: DPT2(MDPT2)
integer(kind=iwp), intent(inout) :: NDPT2
real(kind=wp), intent(in) :: SCAL
integer(kind=iwp) :: iCase1, iCase2, idoff, idpq, idqp, iLoop, iMltOp, ip, iq, iSta, iSym, iSym1, iSym2, lScr, lScr2, lVec1, &
                     lVec2, na, nas1, nas2, ni, nis1, nis2, nLoop, no, nVec1, nVec2, nWec1
real(kind=wp) :: Fact
real(kind=wp), allocatable :: SCR(:), WEC1(:)
#ifdef _MOLCAS_MPP_
real(kind=wp), allocatable :: TMP1(:), TMP2(:)
#endif

! If the G1 correction to the Fock matrix is used, then the
! inactive/virtual coupling elements (which are non-zero for the
! case of average CASSCF) cannot be used in the CASPT2 equations.
if ((FOCKTYPE == 'G1') .and. (.not. G1SECIN)) then
  IFCOUP(12,5) = 0
  IFCOUP(13,5) = 0
end if

NDPT2 = sum(NORB(1:NSYM)**2)

! Loop over ordering: First, <IVEC|...|JVEC>, then reverse.
! For each order, compute the upper-triangular blocks.
! For true density matrices, use symmetry of D-matrix.

! Transform to standard representation, contravariant form.
call PTRTOC(0,IVEC,IVEC)
if (IVEC /= JVEC) call PTRTOC(0,JVEC,JVEC)
NLOOP = 2
if (IVEC == JVEC) NLOOP = 1
do ILOOP=1,NLOOP
  !if (ILOOP == 1) then
  !  IBRA = IVEC
  !  IKET = JVEC
  !else
  !  IBRA = JVEC
  !  IKET = IVEC
  !end if

  ! Loop over types and symmetry block of VEC1 vector:
  do ICASE1=1,13
    do ISYM1=1,NSYM
      if (NINDEP(ISYM1,ICASE1) == 0) cycle
      NIS1 = NISUP(ISYM1,ICASE1)
      NAS1 = NASUP(ISYM1,ICASE1)
      NVEC1 = NIS1*NAS1
      if (NVEC1 == 0) cycle
      ! Form VEC1 from the BRA vector, transformed to covariant form.
      call RHS_ALLO(NAS1,NIS1,LVEC1)
      call RHS_SCAL(NAS1,NIS1,LVEC1,Zero)
      if (ICASE1 <= 11) then
        call RHS_ALLO(NAS1,NIS1,LSCR)
        call RHS_READ(NAS1,NIS1,LSCR,ICASE1,ISYM1,IVEC) !! IBRA)
        if ((IVEC /= JVEC) .and. (ILOOP == 1)) then
          if (SCAL /= One) call RHS_SCAL(NAS1,NIS1,LSCR,SCAL)
          call RHS_ALLO(NAS1,NIS1,LSCR2)
          call RHS_READ(NAS1,NIS1,LSCR2,ICASE1,ISYM1,JVEC)
          call RHS_DAXPY(NAS1,NIS1,One,LSCR2,LSCR)
          call RHS_FREE(LSCR2)
        end if
        call RHS_STRANS(NAS1,NIS1,One,LSCR,LVEC1,ICASE1,ISYM1)
        call RHS_FREE(LSCR)
      else
        call RHS_READ(NAS1,NIS1,LVEC1,ICASE1,ISYM1,IVEC) !! IBRA)
        if ((IVEC /= JVEC) .and. (ILOOP == 1)) then
          if (SCAL /= One) call RHS_SCAL(NAS1,NIS1,LVEC1,SCAL)
          call RHS_ALLO(NAS1,NIS1,LSCR2)
          call RHS_READ(NAS1,NIS1,LSCR2,ICASE1,ISYM1,JVEC)
          call RHS_DAXPY(NAS1,NIS1,One,LSCR2,LVEC1)
          call RHS_FREE(LSCR2)
        end if
      end if
      ! Form WEC1 from VEC1, if needed.
      NWEC1 = 0
      FACT = One/real(max(1,NACTEL),kind=wp)
      if (ICASE1 == 1) NWEC1 = NASH(ISYM1)*NISH(ISYM1)
      if (ICASE1 == 4) NWEC1 = NASH(ISYM1)*NSSH(ISYM1)
      if ((ICASE1 == 5) .and. (ISYM1 == 1)) NWEC1 = NIS1
      if (NWEC1 > 0) then
        call mma_allocate(WEC1,NWEC1,Label='WEC1')
        WEC1(:) = Zero
        IMLTOP = 1
#       ifdef _MOLCAS_MPP_
        if (IS_REAL_PAR()) then
          call mma_allocate(TMP1,NVEC1,Label='TMP1')
          call RHS_GET(NAS1,NIS1,LVEC1,TMP1)
          if (ICASE1 == 1) then
            call SPEC1A(IMLTOP,FACT,ISYM1,TMP1,NVEC1,WEC1,NWEC1)
          else if (ICASE1 == 4) then
            call SPEC1C(IMLTOP,FACT,ISYM1,TMP1,NVEC1,WEC1,NWEC1)
          else if ((ICASE1 == 5) .and. (ISYM1 == 1)) then
            call SPEC1D(IMLTOP,FACT,TMP1,NVEC1,WEC1,NWEC1)
          end if
          call mma_deallocate(TMP1)
        else
#       endif
          if (ICASE1 == 1) then
            call SPEC1A(IMLTOP,FACT,ISYM1,GA_Arrays(LVEC1)%A,size(GA_Arrays(LVEC1)%A),WEC1,NWEC1)
          else if (ICASE1 == 4) then
            call SPEC1C(IMLTOP,FACT,ISYM1,GA_Arrays(LVEC1)%A,size(GA_Arrays(LVEC1)%A),WEC1,NWEC1)
          else if ((ICASE1 == 5) .and. (ISYM1 == 1)) then
            call SPEC1D(IMLTOP,FACT,GA_Arrays(LVEC1)%A,size(GA_Arrays(LVEC1)%A),WEC1,NWEC1)
          end if
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      else
        NWEC1 = 1
        call mma_allocate(WEC1,NWEC1,Label='WEC1')
      end if
      ! Note: WEC1 is identical to <IBRA| E(p,q) |0> for the cases
      ! (p,q)=(t,i), (a,t), and (a,i), resp.
      do ICASE2=ICASE1+1,13
        if (IFCOUP(ICASE2,ICASE1) == 0) cycle
        do ISYM2=1,NSYM
          if (NINDEP(ISYM2,ICASE2) == 0) cycle
          NIS2 = NISUP(ISYM2,ICASE2)
          NAS2 = NASUP(ISYM2,ICASE2)
          NVEC2 = NIS2*NAS2
          if (NVEC2 == 0) cycle
          call RHS_ALLO(NAS2,NIS2,LVEC2)
          call RHS_READ(NAS2,NIS2,LVEC2,ICASE2,ISYM2,IVEC) !! IKET)
          if ((IVEC /= JVEC) .and. (ILOOP == 2)) then
            if (SCAL /= One) call RHS_SCAL(NAS2,NIS2,LVEC2,SCAL)
            call RHS_ALLO(NAS2,NIS2,LSCR2)
            call RHS_READ(NAS2,NIS2,LSCR2,ICASE2,ISYM2,JVEC)
            call RHS_DAXPY(NAS2,NIS2,One,LSCR2,LVEC2)
            call RHS_FREE(LSCR2)
          end if
#         ifdef _MOLCAS_MPP_
          if (IS_REAL_PAR()) then
            call mma_allocate(TMP1,NVEC1,Label='TMP1')
            call mma_allocate(TMP2,NVEC2,Label='TMP2')
            call RHS_GET(NAS1,NIS1,LVEC1,TMP1)
            call RHS_GET(NAS2,NIS2,LVEC2,TMP2)
            call OFFDNS(ISYM1,ICASE1,ISYM2,ICASE2,WEC1,NWEC1,TMP1,nVEC1,DPT2,mDPT2,TMP2,nVEC2,LISTS,size(LISTS))
            call mma_deallocate(TMP1)
            call mma_deallocate(TMP2)
          else
#         endif
            call OFFDNS(ISYM1,ICASE1,ISYM2,ICASE2,WEC1,NWEC1,GA_Arrays(LVEC1)%A,size(GA_Arrays(LVEC1)%A),DPT2,MDPT2, &
                        GA_Arrays(LVEC2)%A,size(GA_Arrays(LVEC2)%A),LISTS,size(LISTS))
#         ifdef _MOLCAS_MPP_
          end if
#         endif
          call RHS_FREE(LVEC2)
        end do
      end do
      call RHS_FREE(LVEC1)
      call mma_deallocate(WEC1)
    end do
  end do

end do

call GADGOP(DPT2,NDPT2,'+')

if (IVEC /= JVEC) then
  ! Transpose the density matrix.
  call mma_allocate(SCR,NDPT2,Label='SCR')
  ISTA = 1
  do ISYM=1,NSYM
    NO = NORB(ISYM)
    call TRNSPS(NO,NO,DPT2(ISTA),SCR)
    DPT2(ISTA:ISTA+NO**2-1) = SCR(1:NO**2)
    ISTA = ISTA+NO**2
  end do
  call mma_deallocate(SCR)
end if

! Transform vectors back to eigenbasis of H0(diag).
call PTRTOSR(1,IVEC,IVEC)
if (IVEC /= JVEC) call PTRTOSR(1,JVEC,JVEC)
if (IVEC == JVEC) then
  ! Fill in lower-triangular block elements by symmetry.
  IDOFF = 0
  do ISYM=1,NSYM
    NI = NISH(ISYM)
    NA = NASH(ISYM)
    NO = NORB(ISYM)
    do IP=1,NI+NA
      do IQ=NI+1,NO
        IDPQ = IDOFF+IP+NO*(IQ-1)
        IDQP = IDOFF+IQ+NO*(IP-1)
        DPT2(IDQP) = DPT2(IDPQ)
      end do
    end do
    IDOFF = IDOFF+NO**2
  end do
end if

end subroutine TRDNS2O
