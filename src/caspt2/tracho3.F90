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
! Copyright (C) Per Ake Malmqvist                                      *
!***********************************************************************

subroutine TRACHO3(CMO,NCMO)

use Symmetry_Info, only: Mul
use CHOVEC_IO, only: ChoVec_Coll, ChoVec_load, ChoVec_Save, npq_ChoType, NVLOC_ChoBatch
use Cholesky, only: InfVec
use ChoCASPT2, only: MxNVc, nChSpc, nFtSpc, nHtSpc, NumCho_PT2
use general_data, only: nAsh
use caspt2_module, only: nBas, nBtch, nBtches, nFro, nInaBx, nIsh, nSecBx, nSsh, nSym, RHSDirect
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NCMO
real(kind=wp), intent(in) :: CMO(NCMO)
integer(kind=iwp) :: IAEND, IASTA, IB, IBATCH, IBATCH_TOT, IBEND, IBSTA, IC, ICASE, IIEND, IISTA, ILOC, ip_htspc, ip_HTVec(8), &
                     IP_LHT, IRC, ISTART(8), ISYM, ISYMA, ISYMB, ISYP, ISYQ, JNUM, JRED, JRED1, JRED2, JREDC, JSTART, JSYM, JV1, &
                     JV2, MUSED, N, N1, N2, NA, NASZ, NBATCH, NBUFFY, NCES(8), NHTOFF, NI, NISZ, NPQ, NUMV, NUSE(8), NVECS_RED
real(kind=wp), allocatable :: BUFFY(:), CHSPC(:), FTSPC(:), HTSPC(:)
#include "warnings.h"

!***********************************************************************
! ======================================================================
! This section deals with density matrices and CMO''s
! Offsets into CMO arrays:
IC = 0
do ISYM=1,NSYM
  NCES(ISYM) = IC
  IC = IC+NBAS(ISYM)**2
end do

! ======================================================================
call mma_allocate(CHSPC,NCHSPC,LABEL='CHSPC')
call mma_allocate(HTSPC,NHTSPC,LABEL='HTSPC')
ip_HTSPC = 1
call mma_allocate(FTSPC,NFTSPC,LABEL='FTSPC')
! ======================================================================

!IBATCH_TOT=0
! Loop over JSYM
do JSYM=1,NSYM
  IBATCH_TOT = NBTCHES(JSYM)

  if (NUMCHO_PT2(JSYM) == 0) cycle

  JRED1 = InfVec(1,2,jSym)
  JRED2 = InfVec(NumCho_PT2(jSym),2,jSym)
  !write(u6,*) 'tracho3:  JRED1,JRED2:',JRED1,JRED2

  ! Loop over JRED
  do JRED=JRED1,JRED2

    call Cho_X_nVecRS(JRED,JSYM,JSTART,NVECS_RED)
    if (NVECS_RED == 0) cycle

    ILOC = 3
    call CHO_X_SETRED(IRC,ILOC,JRED)
    ! For a reduced set, the structure is known, including
    ! the mapping between reduced index and basis set pairs.
    ! The reduced set is divided into suitable batches.
    ! First vector is JSTART. Nr of vectors in r.s. is NVECS_RED.

    ! Determine batch length for this reduced set.
    ! Make sure to use the same formula as in the creation of disk
    ! address tables, etc, above:
    NBATCH = 1+(NVECS_RED-1)/MXNVC

    ! Loop over IBATCH
    JV1 = JSTART
    do IBATCH=1,NBATCH
      IBATCH_TOT = IBATCH_TOT+1

      JNUM = NVLOC_CHOBATCH(IBATCH_TOT)
      JV2 = JV1+JNUM-1

      JREDC = JRED
      ! Read a batch of reduced vectors
      call CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,JSYM,NUMV,JREDC,MUSED)
      if (NUMV /= JNUM) then
        write(u6,*) ' Rats! CHO_VECRD was called, assuming it to'
        write(u6,*) ' read JNUM vectors. Instead it returned NUMV'
        write(u6,*) ' vectors: JNUM, NUMV=',JNUM,NUMV
        write(u6,*) ' Back to the drawing board?'
        call QUIT(_RC_INTERNAL_ERROR_)
      end if
      if (JREDC /= JRED) then
        write(u6,*) ' Rats! It was assumed that the Cholesky vectors'
        write(u6,*) ' in HALFTRNSF all belonged to a given reduced'
        write(u6,*) ' set, but they don''t!'
        write(u6,*) ' JRED, JREDC:',JRED,JREDC
        write(u6,*) ' Back to the drawing board?'
        write(u6,*) ' Let the program continue and see what happens.'
      end if

      ! Frozen half-transformation:
      NHTOFF = 0
      do ISYMA=1,NSYM
        ISYMB = Mul(ISYMA,JSYM)
        IP_HTVEC(ISYMA) = IP_HTSPC+NHTOFF
        ISTART(ISYMA) = 1
        NUSE(ISYMA) = NFRO(ISYMA)
        NHTOFF = NHTOFF+NUSE(ISYMA)*NBAS(ISYMB)*JNUM
      end do
      call HALFTRNSF(IRC,CHSPC,NCHSPC,1,JV1,JNUM,JNUM,JSYM,JREDC,CMO,NCMO,ISTART,NUSE,IP_HTVEC,HTSPC,NHTSPC)

      ! Inactive half-transformation:
      ! Vectors of type HALF(K,J,B) = Sum(CHO(AB,J)*CMO(A,K) where
      ! A,B are basis functions of symmetry ISYMA, ISYMB,
      ! K is inactive of symmetry ISYMA, J is vector number in 1..NUMV
      ! numbered within the present batch.
      ! Symmetry block ISYMA,ISYMB is found at HTSPC(IP_HTVEC(ISYMA)
      NHTOFF = 0
      do ISYMA=1,NSYM
        ISYMB = Mul(ISYMA,JSYM)
        IP_HTVEC(ISYMA) = IP_HTSPC+NHTOFF
        ISTART(ISYMA) = NFRO(ISYMA)+1
        NUSE(ISYMA) = NISH(ISYMA)
        NHTOFF = NHTOFF+NUSE(ISYMA)*NBAS(ISYMB)*JNUM
      end do
      call HALFTRNSF(IRC,CHSPC,NCHSPC,1,JV1,JNUM,JNUM,JSYM,JREDC,CMO,NCMO,ISTART,NUSE,IP_HTVEC,HTSPC,NHTSPC)

      ! Loop over ISYQ
      do ISYQ=1,NSYM
        ISYP = Mul(ISYQ,JSYM)

        N = NBAS(ISYP)
        ! ---------------------------------------------------
        N1 = NASH(ISYP)
        N2 = NISH(ISYQ)
        IC = 1+NCES(ISYP)+(NFRO(ISYP)+NISH(ISYP))*N
        IP_LHT = IP_HTVEC(ISYQ)
        ! Compute fully transformed TK
        if (N1*N2 > 0) then
          call FULLTRNSF(N1,N2,N,CMO(IC),JNUM,HTSPC(IP_LHT),FTSPC)
          call CHOVEC_SAVE(FTSPC,NFTSPC,1,ISYQ,JSYM,IBATCH_TOT)
        end if
        ! ---------------------------------------------------
        N1 = NSSH(ISYP)
        N2 = NISH(ISYQ)
        IC = 1+NCES(ISYP)+(NFRO(ISYP)+NISH(ISYP)+NASH(ISYP))*N
        IP_LHT = IP_HTVEC(ISYQ)
        ! Compute fully transformed AK
        if (N1*N2 > 0) then

          !call FULLTRNSF(N1,N2,N,CMO(IC),JNUM,HTSPC(IP_LHT),FTSPC)

          ! =SVC= modified for using boxed ordering of pairs, note that the boxed
          ! routine is less efficient than the original one (loop over J values)
          NA = N1
          NI = N2
          ! Allocate memory for small buffer used in FULLTRNSF_BOXED
          NBUFFY = NA*NI
          call mma_allocate(BUFFY,NBUFFY,LABEL='BUFFY')
          ! Loop over boxes
          do IASTA=1,NA,nSecBX
            IAEND = min(IASTA-1+nSecBX,NA)
            NASZ = 1+IAEND-IASTA
            do IISTA=1,NI,nInaBX
              IIEND = min(IISTA-1+nInaBX,NI)
              NISZ = 1+IIEND-IISTA
              ! =SVC= note that WITHIN this box, the index of the outer box A (P in
              ! the FULLTRNSF subroutine, i.e. the secondary orbital index) is the
              ! fast-running index, as LFT([A],[I],J) = CMO(P,[A])^T * LHT([I],J,P)^T
              ! with P=1,NB.  So if used in e.g. ADDRHS as BRA(c,l,J), making an inner
              ! loop over secondary orbital index c is more efficient.
              call FULLTRNSF_BOXED(IASTA,IISTA,NASZ,NISZ,NA,NI,N,CMO(IC+N*(IASTA-1)),JNUM,HTSPC(IP_LHT),FTSPC,BUFFY)
            end do
          end do
          call mma_deallocate(BUFFY)
          call CHOVEC_SAVE(FTSPC,NFTSPC,4,ISYQ,JSYM,IBATCH_TOT)
        end if
        ! End loop ISYQ
      end do
      ! ---------------------------------------------------

      ! Active scaled natural orbitals half-transformation:
      ! Vectors of type HALF(W,J,B) = Sum(CHO(AB,J)*CMO(A,W) where
      ! A,B are basis functions of symmetry ISYMA, ISYMB,
      ! W is active of symmetry ISYMA, J is vector number in 1..NUMV
      ! numbered within the present batch.
      ! Symmetry block ISYMA,ISYMB is found at HTSPC(IP_HTVEC(ISYMA)
      NHTOFF = 0
      do ISYMA=1,NSYM
        ISYMB = Mul(ISYMA,JSYM)
        IP_HTVEC(ISYMA) = IP_HTSPC+NHTOFF
        ISTART(ISYMA) = NFRO(ISYMA)+NISH(ISYMA)+1
        NUSE(ISYMA) = NASH(ISYMA)
        NHTOFF = NHTOFF+NUSE(ISYMA)*NBAS(ISYMB)*JNUM
      end do
      ! ---------------------------------------------------
      ! Active half-transformation:
      call HALFTRNSF(IRC,CHSPC,NCHSPC,1,JV1,JNUM,JNUM,JSYM,JREDC,CMO,NCMO,ISTART,NUSE,IP_HTVEC,HTSPC,NHTSPC)

      do ISYQ=1,NSYM
        ISYP = Mul(ISYQ,JSYM)

        N = NBAS(ISYP)
        ! ---------------------------------------------------
        ! Loop over ISYQ
        N1 = NASH(ISYP)
        N2 = NASH(ISYQ)
        IC = 1+NCES(ISYP)+(NFRO(ISYP)+NISH(ISYP))*N
        IP_LHT = IP_HTVEC(ISYQ)
        ! Compute fully transformed TV
        if (N1*N2 > 0) then
          call FULLTRNSF(N1,N2,N,CMO(IC),JNUM,HTSPC(IP_LHT),FTSPC)
          call CHOVEC_SAVE(FTSPC,NFTSPC,2,ISYQ,JSYM,IBATCH_TOT)
        end if
        ! ---------------------------------------------------
        N1 = NSSH(ISYP)
        N2 = NASH(ISYQ)
        IC = 1+NCES(ISYP)+(NFRO(ISYP)+NISH(ISYP)+NASH(ISYP))*N
        IP_LHT = IP_HTVEC(ISYQ)
        ! Compute fully transformed AV
        if (N1*N2 > 0) then
          call FULLTRNSF(N1,N2,N,CMO(IC),JNUM,HTSPC(IP_LHT),FTSPC)
          call CHOVEC_SAVE(FTSPC,NFTSPC,3,ISYQ,JSYM,IBATCH_TOT)
        end if
        ! ---------------------------------------------------
        ! End loop ISYQ
      end do

      ! End loop IBATCH
      JV1 = JV1+JNUM
    end do

    ! End loop JRED
  end do

! End loop JSYM
end do

! if using the RHS on-demand, we need all cholesky vectors on each
! process, collect them here
if (RHSDIRECT) then
  do JSYM=1,NSYM
    IBSTA = NBTCHES(JSYM)+1
    IBEND = NBTCHES(JSYM)+NBTCH(JSYM)
    do IB=IBSTA,IBEND
      do ISYQ=1,NSYM
        do ICASE=1,4
          NPQ = NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
          if (NPQ == 0) cycle
          call CHOVEC_LOAD(FTSPC,NFTSPC,ICASE,ISYQ,JSYM,IB)
          call CHOVEC_COLL(FTSPC,NFTSPC,ICASE,ISYQ,JSYM,IB)
        end do
      end do
    end do
  end do
end if

call mma_deallocate(CHSPC)
call mma_deallocate(HTSPC)
call mma_deallocate(FTSPC)

end subroutine TRACHO3
