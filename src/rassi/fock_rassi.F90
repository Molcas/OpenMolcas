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

subroutine FOCK_RASSI(DINAO,FOCKAO)
! Purpose:
!  ADD THE TWO-ELECTRON PART OF A FOCK MATRIX USING THE
!  INACTIVE DENSITY MATRIX. THE FOCK MATRIX SHOULD CONTAIN THE
!  ONE-ELECTRON HAMILTONIAN MATRIX BEFORE THE CALL. THE MATRICES
!  ARE STORED IN SYMMETRY-BLOCKED SQUARE FORMAT.

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: MUL, nIrrep
use rassi_data, only: NBASF, NBMX, NBSQ, NBSQPR, NBTRI, NISH
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: DINAO(NBSQ)
real(kind=wp), intent(inout) :: FOCKAO(NBSQ)
integer(kind=iwp) :: IDF, IDPQ, IDRS, IFPQ, IFRS, II, INTBUF, INUSE, IOPT, IPQ, IRC, IRSST, ISP, ISQ, ISR, ISS, ISSMX, ISYM, &
                     KEEP(8), KEEPP, KEEPQ, KEEPR, KEEPS, KEEPT, LPQ, NB, NBP, NBPQ, NBQ, NBR, NBRS, NBS, NBSX(8), NBT, NBTR, &
                     NBTRPR(8), NFTRI, NIP, NIQ, NIR, NIS, NP, NPQ, NPQM, NQ, NQM, NSPQ, NSPQR, NSQBUF, NSYMX
real(kind=wp) :: DF, FPQ
logical(kind=iwp) :: ISQARX
real(kind=wp), allocatable :: DTRI(:), FTRI(:), PQRS(:), SQBUF(:)

! RETRIEVE BASE DATA FROM UNIT LUORD:
IRC = 0
call GETORD(IRC,ISQARX,NSYMX,NBSX,KEEP)
! ALLOCATE WORK SPACE FOR FOLDED DENSITY MATRIX, FOR ONE
! TRIANGULAR SYMMETRY-BLOCK OF THE TWO-ELECTRON CONTRIBUTION,
! FOR AN INTEGRAL BUFFER, AND FOR EXPANDING TRIANGULAR INTEGRAL
! MATRICES INTO SQUARE FORMAT:
NSQBUF = NBMX**2
INTBUF = max(NSQBUF,256*256)
call mma_allocate(PQRS,INTBUF,Label='PQRS')
call mma_allocate(DTRI,NBTRI,Label='DTRI')
call mma_allocate(SQBUF,NSQBUF,Label='SQBUF')
!PAM00 Nr of triangular matrices in previous symmetry blocks
NBTR = 0
do ISYM=1,nIrrep
  NBTRPR(ISYM) = NBTR
  NB = NBASF(ISYM)
  NBTR = NBTR+nTri_Elem(NB)
end do
NFTRI = NBTR
call mma_allocate(FTRI,NFTRI,Label='FTRI')
FTRI(:) = Zero
! FOLD THE D MATRIX:
IDF = 1
do ISYM=1,nIrrep
  NB = NBASF(ISYM)
  NPQ = NBSQPR(ISYM)
  do NP=1,NB
    do NQ=1,NP
      DF = DINAO(NPQ+NB*(NP-1)+NQ)+DINAO(NPQ+NB*(NQ-1)+NP)
      if (NP == NQ) DF = Half*DF
      DTRI(IDF) = DF
      IDF = IDF+1
    end do
  end do
end do
! LOOP OVER ISP:
do ISP=1,nIrrep
  NBP = NBASF(ISP)
  NIP = NISH(ISP)
  KEEPP = KEEP(ISP)
  ! LOOP OVER ISQ:
  do ISQ=1,ISP
    NBQ = NBASF(ISQ)
    NIQ = NISH(ISQ)
    NSPQ = MUL(ISP,ISQ)
    KEEPQ = KEEP(ISQ)
    ! LOOP OVER ISR:
    do ISR=1,ISP
      NBR = NBASF(ISR)
      NIR = NISH(ISR)
      NSPQR = MUL(NSPQ,ISR)
      KEEPR = KEEP(ISR)
      ! LOOP OVER ISS:
      ISSMX = ISR
      if (ISP == ISR) ISSMX = ISQ
      do ISS=1,ISSMX
        if (ISS /= NSPQR) cycle
        NBS = NBASF(ISS)
        NIS = NISH(ISS)
        KEEPS = KEEP(ISS)
        KEEPT = KEEPP+KEEPQ+KEEPR+KEEPS
        NBT = NBP*NBQ*NBR*NBS
        if ((KEEPT /= 0) .and. (NBT /= 0)) call ABEND()
        if (NBT == 0) cycle
        INUSE = 0
        if ((ISP == ISQ) .and. (NBP > 0) .and. (NIR > 0)) INUSE = 1
        if ((ISP == ISQ) .and. (NIP > 0) .and. (NBR > 0)) INUSE = 1
        if ((ISP == ISR) .and. (NBP > 0) .and. (NIQ > 0)) INUSE = 1
        if ((ISP == ISR) .and. (NIP > 0) .and. (NBQ > 0)) INUSE = 1
        if (INUSE == 0) cycle

        ! SIZES OF INTEGRAL MATRICES:
        NBPQ = NBP*NBQ
        if (ISP == ISQ) NBPQ = nTri_Elem(NBP)
        NBRS = NBR*NBS
        if (ISR == ISS) NBRS = nTri_Elem(NBR)
        ! READ THE AO INTEGRALS INTO PQRS WHENEVER NEEDED.
        ! ACCUMULATE CONTRIBUTIONS TO FOCK MATRIX INTO FTRI, OR
        ! DIRECTLY INTO FOCKAO.
        IRC = 0
        IOPT = 1
        IPQ = 0
        LPQ = 0
        NPQ = 0
        IRSST = 1-NBRS
        do NP=1,NBP
          NQM = NBQ
          if (ISP == ISQ) NQM = NP
          do NQ=1,NQM
            IPQ = IPQ+1
            ! READ A NEW BUFFER OF INTEGRALS
            if (LPQ == NPQ) then
              call RDORD(IRC,IOPT,ISP,ISQ,ISR,ISS,PQRS,INTBUF,NPQ)
              IOPT = 2
              LPQ = 0
              IRSST = 1-NBRS
              ! COULOMB CONTRIBUTION:
              if (ISP == ISQ) then
                if (NIR > 0) then
                  IFPQ = NBTRPR(ISP)+IPQ
                  IDRS = NBTRPR(ISR)+1
                  NPQM = min(NPQ,NBPQ-IPQ+1)
                  call DGEMV_('T',NBRS,NPQM,One,PQRS,NBRS,DTRI(IDRS),1,One,FTRI(IFPQ),1)
                end if
                !PAM00 Added code, to obtain RSPQ contributions from the PQRS integrals
                if ((ISP > ISR) .and. (NIP > 0)) then
                  IFRS = NBTRPR(ISR)+1
                  IDPQ = NBTRPR(ISP)+IPQ
                  NPQM = min(NPQ,NBPQ-IPQ+1)
                  call DGEMV_('N',NBRS,NPQM,One,PQRS,NBRS,DTRI(IDPQ),1,One,FTRI(IFRS),1)
                end if
              end if
            end if
            LPQ = LPQ+1
            IRSST = IRSST+NBRS
            ! EXCHANGE CONTRIBUTIONS:
            if (ISP == ISR) then
              if (ISR == ISS) call SQUARE(PQRS(IRSST),SQBUF,1,NBR,NBR)
              if (((ISP /= ISQ) .or. (NP /= NQ)) .and. (NIR > 0)) then
                if (ISR == ISS) then
                  call DGEMV_('N',NBS,NBR,-Half,SQBUF,NBS,DINAO(NBSQPR(ISP)+NBR*(NP-1)+1),1,One,FOCKAO(NBSQPR(ISQ)+NQ),NBQ)
                else
                  call DGEMV_('N',NBS,NBR,-Half,PQRS(IRSST),NBS,DINAO(NBSQPR(ISP)+NBR*(NP-1)+1),1,One,FOCKAO(NBSQPR(ISQ)+NQ),NBQ)
                end if
              end if
              if (ISP == ISS) then
                if (NIR > 0) &
                  call DGEMV_('N',NBS,NBR,-Half,SQBUF,NBS,DINAO(NBSQPR(ISQ)+NBR*(NQ-1)+1),1,One,FOCKAO(NBSQPR(ISP)+NP),NBP)
              else if (NIS > 0) then
                call DGEMV_('T',NBS,NBR,-Half,PQRS(IRSST),NBS,DINAO(NBSQPR(ISQ)+NBS*(NQ-1)+1),1,One,FOCKAO(NBSQPR(ISP)+NP),NBP)
              end if
            end if
          end do  ! nq
        end do  ! np

      end do  ! iss
    end do  ! isr
  end do  ! isq
end do ! isp
! END OF QUADRUPLE SYMMETRY LOOP.
! ADD COULOMB CONTRIBUTIONS FROM TRIANGULAR STORAGE TO FOCKAO:
!PAM00 OK -- Now add Coulomb contributions:
IFPQ = 0
do ISP=1,nIrrep
  NBP = NBASF(ISP)
  do NP=1,NBP
    do NQ=1,NP
      IFPQ = IFPQ+1
      FPQ = FTRI(IFPQ)
      II = NBSQPR(ISP)+NBP*(NP-1)+NQ
      FOCKAO(II) = FOCKAO(II)+FPQ
      if (NQ /= NP) then
        II = NBSQPR(ISP)+NBP*(NQ-1)+NP
        FOCKAO(II) = FOCKAO(II)+FPQ
      end if
    end do
  end do
end do
call mma_deallocate(PQRS)
call mma_deallocate(DTRI)
call mma_deallocate(SQBUF)
call mma_deallocate(FTRI)

call GADGOp(FOCKAO,NBSQ,'+')

end subroutine FOCK_RASSI
