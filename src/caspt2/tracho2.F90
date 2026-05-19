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

subroutine TRACHO2(CMO,NCMO,DREF,NDREF,FFAO,FIAO,FAAO,IF_TRNSF)

use Symmetry_Info, only: Mul
use CHOVEC_IO, only: chovec_coll, chovec_load, chovec_save, NPQ_CHOTYPE, NVLOC_CHOBATCH
use Cholesky, only: InfVec, nDimRS
use ChoCASPT2, only: MXCHARR, MXNVC, NCHSPC, NFTSPC, NHTSPC, NUMCHO_PT2
use caspt2_module, only: nAsh, nBas, nBasT, nBSqT, nBtch, nBtches, nBTri, nFro, nInaBx, nIsh, nSecBx, nSsh, nSym, RHSDirect
#ifdef _DEBUGPRINT_
use caspt2_module, only: PotNuc
use Definitions, only: u6
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NCMO, NDREF
real(kind=wp), intent(in) :: CMO(NCMO), DREF(NDREF)
real(kind=wp), intent(out) :: FFAO(NBTRI), FIAO(NBTRI), FAAO(NBTRI)
logical(kind=iwp), intent(in) :: IF_TRNSF
integer(kind=iwp) :: I, IA, IAEND, IASTA, IB, IBATCH, IBATCH_TOT, IBEND, IBSTA, IC, ICASE, IDAIJ, IDFIJ, IDIIJ, IIEND, IISTA, &
                     ILOC, ip_htspc, ip_HTVec(8), IP_LHT, IRC, ISFA, ISFF, ISFI, ISTART(8), ISYM, ISYMA, ISYMB, ISYMK, ISYMW, &
                     ISYP, ISYQ, J, JNUM, JRED, JRED1, JRED2, JREDC, JSTART, JSYM, JV1, JV2, LC, LO, LSC, LSO, MUSED, N, N1, N2, &
                     NA, NASZ, NB, NBATCH, NBUFFY, NCES(8), NF, NHTOFF, NI, NISZ, NK, NPQ, NRS, NUMV, NUSE(8), NVECS_RED, NW
real(kind=wp) :: FACTC, FACTXA, FACTXI
real(kind=wp), allocatable :: BUFFY(:), CHSPC(:), CNAT(:), DA(:), DA_RED(:), DF(:), DF_RED(:), DI(:), DI_RED(:), FA_RED(:), &
                              FF_RED(:), FI_RED(:), FTSPC(:), HTSPC(:), OCC(:), VEC(:)
#ifdef _DEBUGPRINT_
real(kind=wp) :: E, ECORE, ECORE1, ECORE2
real(kind=wp), external :: DDOT_
#endif
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

! Compute natural orbitals for the reference wave function:
call mma_allocate(OCC,NBasT,Label='OCC')
call mma_allocate(CNAT,NBSQT,Label='CNAT')
!write(u6,*)' Active/Active density matrix, triangular'
!if (NASHT > 0 ) then
!  call TRIPRT(' ',' ',DREF,NASHT)
!end if
call REF_NATO(DREF,nDREF,CMO,nCMO,OCC,NBasT,CNAT,NBSQT)

! Initialize Fock matrices in AO basis to zero:
FFAO(:) = Zero
FIAO(:) = Zero
FAAO(:) = Zero
! Construct density matrix for frozen orbitals
call mma_allocate(DF,NBTRI,Label='DF')
do ISYM=1,NSYM
  ISTART(ISYM) = 1
  NUSE(ISYM) = NFRO(ISYM)
end do
call GDMAT(NSYM,NBAS,ISTART,NUSE,CNAT,NBSQT,OCC,NBasT,DF,NBTRI)
! Construct density matrix for inactive orbitals
call mma_allocate(DI,NBTRI,Label='DI')
do ISYM=1,NSYM
  ISTART(ISYM) = NFRO(ISYM)+1
  NUSE(ISYM) = NISH(ISYM)
end do
call GDMAT(NSYM,NBAS,ISTART,NUSE,CNAT,NBSQT,OCC,NBasT,DI,NBTRI)
! Same, for active density:
call mma_allocate(DA,NBTRI,Label='DA')
do ISYM=1,NSYM
  ISTART(ISYM) = NFRO(ISYM)+NISH(ISYM)+1
  NUSE(ISYM) = NASH(ISYM)
end do
call GDMAT(NSYM,NBAS,ISTART,NUSE,CNAT,NBSQT,OCC,NBasT,DA,NBTRI)
! The Cholesky routines want density matrices in a particular storage, and
! also the off-diagonal elements should be doubled. Double them:
IDFIJ = 1
IDIIJ = 1
IDAIJ = 1
do ISYM=1,NSYM
  do I=1,NBAS(ISYM)
    do J=1,I-1
      DF(IDFIJ) = Two*DF(IDFIJ)
      DI(IDIIJ) = Two*DI(IDIIJ)
      DA(IDAIJ) = Two*DA(IDAIJ)
      IDFIJ = IDFIJ+1
      IDIIJ = IDIIJ+1
      IDAIJ = IDAIJ+1
    end do
    IDFIJ = IDFIJ+1
    IDIIJ = IDIIJ+1
    IDAIJ = IDAIJ+1
  end do
end do

! Scale natural orbitals by multiplying with square root of half the
! occupation number -- This allows computing the exchange contribution
! to Fock matrices using the same formula as for closed shells.
LSO = 1
LSC = 1
do ISYM=1,NSYM
  NF = NFRO(ISYM)
  NI = NISH(ISYM)
  NA = NASH(ISYM)
  NB = NBAS(ISYM)
  LO = LSO+NF+NI
  LC = LSC+NB*(NF+NI)
  do IA=1,NA
    CNAT(LC:LC+NB-1) = sqrt(Half*OCC(LO))*CNAT(LC:LC+NB-1)
    LO = LO+1
    LC = LC+NB
  end do
  LSO = LSO+NB
  LSC = LSC+NB**2
end do
! ======================================================================
call mma_allocate(CHSPC,NCHSPC,LABEL='CHSPC')
call mma_allocate(HTSPC,NHTSPC,LABEL='HTSPC')
IP_HTSPC = 1
if (IF_TRNSF) call mma_allocate(FTSPC,NFTSPC,LABEL='FTSPC')
! ======================================================================

!IBATCH_TOT = 0
! Loop over JSYM
do JSYM=1,NSYM
  IBATCH_TOT = NBTCHES(JSYM)

  !write(u6,*) ' Tracho2 JSYM=',JSYM
  !write(u6,*) '    NUMCHO_PT2(JSYM)=',NUMCHO_PT2(JSYM)
  if (NUMCHO_PT2(JSYM) == 0) cycle

  JRED1 = InfVec(1,2,jSym)
  JRED2 = InfVec(NumCho_PT2(jSym),2,jSym)
  !write(u6,*) 'tracho2:  JRED1,JRED2:',JRED1,JRED2

  if (JSYM == 1) then
    ! Allocate space for temporary vector 'Vec' used for Coulomb contrib to
    ! Fock matrices:
    call mma_allocate(VEC,MXCHARR,LABEL='VEC')
    ! Local density matrices, which will be needed if JSYM=1. At the same time,
    ! allocate Fock matrices with the same structure and initialize to zero.
    call mma_allocate(DF_RED,MXCHARR,LABEL='DF_RED')
    call mma_allocate(DI_RED,MXCHARR,LABEL='DI_RED')
    call mma_allocate(DA_RED,MXCHARR,LABEL='DA_RED')
    call mma_allocate(FF_RED,MXCHARR,LABEL='FF_RED')
    call mma_allocate(FI_RED,MXCHARR,LABEL='FI_RED')
    call mma_allocate(FA_RED,MXCHARR,LABEL='FA_RED')
  end if

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
    !JEND = JSTART+NVECS_RED-1
    !write(u6,*) '  JRED:  JSTART,JEND:',JRED,JSTART,JEND

    if (JSYM == 1) then
      NRS = NDIMRS(JSYM,JRED)
      call full2red(DF,NBTRI,DF_Red,nRS)
      call full2red(DI,NBTRI,DI_Red,nRS)
      call full2red(DA,NBTRI,DA_Red,nRS)
      FF_RED(1:nRS) = Zero
      FI_RED(1:nRS) = Zero
      FA_RED(1:nRS) = Zero
    end if

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

      if (JSYM == 1) then
        ! Coulomb contribution to Fock arrays.
        ! V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * D(rs)
        ! Starting at CHSPC is now an array of vectors, conceptually
        ! L(rs,J), where temporarily we can regard J as ranging 1..JNUM, and
        ! the layout of pair indices rs is unknown ('reduced storage', a secret
        ! inside cholesky.) Compute array V(J) at temporary space VEC:
        call DGEMV_('T',NRS,JNUM,One,CHSPC,NRS,DF_RED,1,Zero,VEC,1)
        ! F(rs){#J} <- F(rs){#J} + FactC * sum_J L(rs,{#J})*V{#J}
        FactC = One
        call DGEMV_('N',NRS,JNUM,FactC,CHSPC,NRS,VEC,1,One,FF_RED,1)
        ! The same thing, now for the inactive and active density matrices:
        call DGEMV_('T',NRS,JNUM,One,CHSPC,NRS,DI_RED,1,Zero,VEC,1)
        call DGEMV_('N',NRS,JNUM,FactC,CHSPC,NRS,VEC,1,One,FI_RED,1)
        call DGEMV_('T',NRS,JNUM,One,CHSPC,NRS,DA_RED,1,Zero,VEC,1)
        call DGEMV_('N',NRS,JNUM,FactC,CHSPC,NRS,VEC,1,One,FA_RED,1)
        !write(u6,*) ' Finished Coulomb contributions to Fock matrix.'
        !write(u6,*) ' Frozen Fock mat at FF_RED'
        !write(u6,'(1x,8f10.4)') (FF_RED(i),i=1,nRS)
        !write(u6,*) ' Inactive Fock mat at FI_RED'
        !write(u6,'(1x,8f10.4)') (FI_RED(i),i=1,nRS)
        !write(u6,*) ' Active Fock matrix at FA_RED.'
        !write(u6,'(1x,8f10.4)') (FA_RED(i),i=1,nRS)
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
      ! Frozen contributions to exchange:
      FactXI = -One
      ISFF = 1
      do ISYMB=1,NSYM
        iSymk = Mul(jSym,iSymb)
        ! --------------------------------------------------------------
        ! *** Compute the LT part of the FROZEN exchange matrix ********
        !     FF(ab) = FF(ab) + FactXI * sum_Jk  LkJ,a * LkJ,b
        ! --------------------------------------------------------------
        NK = NFRO(ISYMK)
        NB = NBAS(ISYMB)
        if (NB*NK /= 0) &
          call DGEMM_TRI('T','N',NB,NB,NK*JNUM,FactXI,HTSPC(ip_HTVec(iSymk)),NK*JNUM,HTSPC(ip_HTVec(iSymk)),NK*JNUM,One, &
                         FFAO(ISFF),NB)
        ISFF = ISFF+(NB*(NB+1))/2
      end do

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
      ! Inactive contributions to exchange:
      FactXI = -One
      ISFI = 1
      do ISYMB=1,NSYM
        iSymk = Mul(jSym,iSymb)
        ! --------------------------------------------------------------
        ! *** Compute the LT part of the INACTIVE exchange matrix ******
        !     FI(ab) = FI(ab) + FactXI * sum_Jk  LkJ,a * LkJ,b
        ! --------------------------------------------------------------
        NK = NISH(iSymk)
        NB = NBAS(ISYMB)
        if (NB*NK /= 0) &
          call DGEMM_TRI('T','N',NB,NB,NK*JNUM,FactXI,HTSPC(ip_HTVec(iSymk)),NK*JNUM,HTSPC(ip_HTVec(iSymk)),NK*JNUM,One, &
                         FIAO(ISFI),NB)
        ISFI = ISFI+(NB*(NB+1))/2
      end do
      !write(u6,*) ' Inactive Fock mat in FIAO'
      !write(u6,'(1x,8f10.4)') (FIAO(i),i=1,nbtri)

      if (IF_TRNSF) then
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
            call mma_allocate(BUFFY,NBUFFY,Label='BUFFY')
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
          ! ---------------------------------------------------
          ! End loop ISYQ
        end do
      end if
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
      call HALFTRNSF(IRC,CHSPC,NCHSPC,1,JV1,JNUM,JNUM,JSYM,JREDC,CNAT,NBSQT,ISTART,NUSE,IP_HTVEC,HTSPC,NHTSPC)
      ! Active (scaled) contributions to exchange:
      FactXA = -One
      ISFA = 1
      do ISYMB=1,NSYM
        iSymw = Mul(jSym,iSymb)
        ! --------------------------------------------------------------
        ! *** Compute the LT part of the ACTIVE exchange matrix ********
        !     FA(ab) = FA(ab) + FactXA * sum_Jw  LwJ,a * LwJ,b
        ! --------------------------------------------------------------
        NW = NASH(iSymw)
        NB = NBAS(ISYMB)
        if (NB*NW /= 0) &
          call DGEMM_TRI('T','N',NB,NB,NW*JNUM,FactXA,HTSPC(ip_HTVec(iSymw)),NW*JNUM,HTSPC(ip_HTVec(iSymw)),NW*JNUM,One, &
                         FAAO(ISFA),NB)
        ISFA = ISFA+(NB*(NB+1))/2
      end do

      !write(u6,*) ' Active Fock matrix in FAAO.'
      !write(u6,'(1x,8f10.4)') (FAAO(i),i=1,nbtri)

      ! ---------------------------------------------------
      ! Active half-transformation:
      call HALFTRNSF(IRC,CHSPC,NCHSPC,1,JV1,JNUM,JNUM,JSYM,JREDC,CMO,NCMO,ISTART,NUSE,IP_HTVEC,HTSPC,NHTSPC)

      if (IF_TRNSF) then
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
      end if

      ! End loop IBATCH
      JV1 = JV1+JNUM
    end do

    if (jSym == 1) then
      ! Add Coulomb contributions in local Fock matrices (in 'reduced storage')
      ! into the global ones:
      call red2full(FFAO,NBTRI,FF_RED,nRS)
      call red2full(FIAO,NBTRI,FI_RED,nRS)
      call red2full(FAAO,NBTRI,FA_RED,nRS)
    end if
    ! End loop JRED
  end do

  if (jSym == 1) then
    ! Deallocate local density and fock matrices
    call mma_deallocate(VEC)
    call mma_deallocate(DF_RED)
    call mma_deallocate(DI_RED)
    call mma_deallocate(DA_RED)
    call mma_deallocate(FF_RED)
    call mma_deallocate(FI_RED)
    call mma_deallocate(FA_RED)
  end if
  ! End loop JSYM
end do

! if using the RHS on-demand, we need all cholesky vectors on each
! process, collect them here
if (IF_TRNSF .and. RHSDIRECT) then
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

! Synchronize and add the contributions from all nodes into each node:
call GADGOP(FFAO,NBTRI,'+')
call GADGOP(FIAO,NBTRI,'+')
call GADGOP(FAAO,NBTRI,'+')

! Add OneHam to finalize frozen Fock matrix in AO basis.
! (It is in fact an effective one-electron Hamiltonian).
call ADD1HAM(FFAO,NBTRI)
#ifdef _DEBUGPRINT_
! Two-electron contribution to the effective core energy
ECORE2 = Half*DDOT_(NBTRI,DF,1,FFAO,1)
! The contraction of frozen Fock matrix with frozen density:
E = DDOT_(NBTRI,DF,1,FFAO,1)
! Correct for double-counting two-electron part:
E = E-ECORE2
! One-electron part:
ECORE1 = E-ECORE2
! Nuclear repulsion energy:
ECORE = POTNUC+ECORE1+ECORE2

write(u6,'(6X,A,ES20.10)') 'NUCLEAR REPULSION ENERGY:',POTNUC
write(u6,'(6X,A,ES20.10)') 'ONE-ELECTRON CORE ENERGY:',ECORE1
write(u6,'(6X,A,ES20.10)') 'TWO-ELECTRON CORE ENERGY:',ECORE2
write(u6,'(6X,A,ES20.10)') '       TOTAL CORE ENERGY:',ECORE
#endif

call mma_deallocate(OCC)
call mma_deallocate(CNAT)
call mma_deallocate(DF)
call mma_deallocate(DI)
call mma_deallocate(DA)

call mma_deallocate(CHSPC)
call mma_deallocate(HTSPC)
if (IF_TRNSF) call mma_deallocate(FTSPC)

#ifdef _DEBUGPRINT_
write(u6,'(6X,A)') 'TEST PRINT FROM TRACHO2.'
write(u6,'(6X,A)')
write(u6,*) ' NSYM:',NSYM
write(u6,*) ' NBAS:',(NBAS(ISYM),ISYM=1,8)
write(u6,'(6X,A)')
write(u6,'(6X,A)') '***** FROZEN FOCK MATRIX ***** '
ISFF = 1
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  if (NB > 0) then
    write(u6,'(6X,A)')
    write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FFAO(ISFF),NB)
    ISFI = ISFF+(NB*(NB+1))/2
  end if
end do
write(u6,'(6X,A)')
write(u6,'(6X,A)') '***** INACTIVE FOCK MATRIX ***** '
ISFI = 1
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  if (NB > 0) then
    write(u6,'(6X,A)')
    write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FIAO(ISFI),NB)
    ISFI = ISFI+(NB*(NB+1))/2
  end if
end do
write(u6,'(6X,A)')
write(u6,'(6X,A)') '***** ACTIVE FOCK MATRIX ***** '
ISFA = 1
do ISYM=1,NSYM
  NB = NBAS(ISYM)
  if (NB > 0) then
    write(u6,'(6X,A)')
    write(u6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',FAAO(ISFA),NB)
    ISFA = ISFA+(NB*(NB+1))/2
  end if
end do
#endif

end subroutine TRACHO2
