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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************
!  CHO_get_ER
!
!> @brief
!>   Compute the Edmiston--Ruedenberg functional for a given set of occupied MOs
!> @author F. Aquilante
!>
!> @details
!> Computes the Edmiston--Ruedenberg functional
!>
!> \f[ W = \sum_i \mathit{ER}[i] = \sum_i (ii|ii) \f]
!>
!> for a given set of occupied MOs.
!>
!> The functional orbital components \p ER(i) are
!> computed by using the Cholesky representation \f$ L_{ab,J} \f$
!> of the AO two-electron integrals, namely
!>
!> \f[ \mathit{ER}[i] = \sum_J V[i]_J V[i]_J \f]
!>
!> where
!>
!> \f[ V[i]_J = \sum_{ab} D[i]_{ab} L_{ab,J} \f]
!>
!> and
!>
!> \f[ D[i]_{a,b} = C[i]_a C[i]_b \f]
!>
!> @note
!> Requires initialization of the Cholesky information.
!>
!> @param[out] irc     return code
!> @param[in]  CMO     MOs matrix, stored as \p C(a,k)
!> @param[in]  nOcc    number of occupied orbitals in each symmetry
!> @param[out] ER      orbital components of the ER functional}
!> @param[out] W       value of the Edmiston--Ruedenberg functional
!> @param[in]  timings switch on/off timings printout
!***********************************************************************

subroutine CHO_get_ER(irc,CMO,nOcc,ER,W,timings)

use Index_Functions, only: iTri, nTri_Elem
use Cholesky, only: InfVec, nBas, nDimRS, nSym, NumCho
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: irc
real(kind=wp), intent(in) :: CMO(*)
integer(kind=iwp), intent(in) :: nOcc(*)
real(kind=wp), intent(_OUT_) :: ER(*)
real(kind=wp), intent(out) :: W
logical(kind=iwp), intent(in) :: timings
integer(kind=iwp) :: ia, ib, iBatch, ik, iLoc, iOcc(8), ipa, ipab, ipb, isMO(8), IVEC2, iVrs, JNUM, JRED, JRED_, JRED1, JRED2, &
                     JSYM, jv, JVEC, kSym, LREAD, LWORK, MaxB, MaxBB, MUSED, nBatch, nOccT, nRS, NUMV, nVec, nVrs
real(kind=wp) :: TCI1, TCI2, TCR1, TCR2, tintg(2), TOTCPU, TOTCPU1, TOTCPU2, TOTWALL, TOTWALL1, TOTWALL2, tread(2), TWI1, TWI2, &
                 TWR1, TWR2
real(kind=wp), allocatable :: Dab(:), DLT(:), Lab(:), VJ(:)
character(len=*), parameter :: SECNAM = 'CHO_get_ER'

JSYM = 1
if (NumCho(JSYM) < 1) then
  write(u6,*) SECNAM//'No total symmetric vectors present'
  irc = 77
  return
end if

call CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

! 1 --> CPU   2 --> Wall
tread(:) = Zero  !time for reading the vectors
tintg(:) = Zero  !time for computing the functional

! compute some offsets and other quantities
MaxB = nBas(1)
iOcc(1) = 0
isMO(1) = 0
do kSym=2,nSym
  iOcc(kSym) = iOcc(kSym-1)+nOcc(kSym-1)
  isMO(kSym) = isMO(kSym-1)+nBas(kSym-1)**2
  MaxB = max(MaxB,nBas(kSym))
end do

MaxBB = nTri_Elem(MaxB)
nOccT = iOcc(nSym)+nOcc(nSym)

W = Zero ! initialization of the ER-functional value
ER(1:nOccT) = Zero ! and its orbital components

call mma_allocate(DLT,MaxBB,Label='DLT')
DLT(:) = Zero

iLoc = 3 ! use scratch location in reduced index arrays

JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

do JRED=JRED1,JRED2

  ! Memory management section -----------------------------

  call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

  if (nVrs == 0) cycle

  if (nVrs < 0) then
    write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs < 0. STOP!!'
    call abend()
  end if

  call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
  if (irc /= 0) then
    write(u6,*) SECNAM//'cho_X_setred non-zero return code. rc= ',irc
    call abend()
  end if

  nRS = nDimRS(JSYM,JRED)

  call mma_allocate(Dab,nRS,Label='Dab')

  call mma_maxDBLE(LWORK)

  nVec = min(LWORK/(nRS+1),nVrs)

  if (nVec < 1) then
    write(u6,*) SECNAM//': Insufficient memory for batch'
    write(u6,*) 'LWORK= ',LWORK
    write(u6,*) 'min. mem. need= ',nRS+1
    irc = 33
    call Abend()
    nBatch = -9999  ! dummy assignment
  end if

  LREAD = nRS*nVec

  call mma_allocate(Lab,LREAD,Label='Lab')
  call mma_allocate(VJ,nVec,Label='VJ')

  ! BATCH over the vectors in JSYM=1 ----------------------------

  nBatch = (nVrs-1)/nVec+1

  do iBatch=1,nBatch

    if (iBatch == nBatch) then
      JNUM = nVrs-nVec*(nBatch-1)
    else
      JNUM = nVec
    end if

    JVEC = nVec*(iBatch-1)+iVrs
    IVEC2 = JVEC-1+JNUM

    call CWTIME(TCR1,TWR1)

    JRED_ = JRED
    call CHO_VECRD(Lab,LREAD,JVEC,IVEC2,JSYM,NUMV,JRED_,MUSED)

    if ((NUMV <= 0) .or. (NUMV /= JNUM)) then
      irc = 77
      return
    end if

    call CWTIME(TCR2,TWR2)
    tread(1) = tread(1)+(TCR2-TCR1)
    tread(2) = tread(2)+(TWR2-TWR1)

    call CWTIME(TCI1,TWI1)

    do kSym=1,nSym

      do ik=1,nOcc(kSym)

        ! Compute the i-th orbital component D[i](a,b) of the density
        ! D[i](a,b) is stored as upper-triangular
        do ib=1,nBas(kSym)

          ipb = isMO(kSym)+nBas(kSym)*(ik-1)+ib

          do ia=1,ib-1

            ipa = isMO(kSym)+nBas(kSym)*(ik-1)+ia
            ipab = iTri(ib,ia)
            DLT(ipab) = CMO(ipa)*CMO(ipb)

          end do

          ipab = nTri_Elem(ib)
          DLT(ipab) = Half*CMO(ipb)**2  !diagonal scaled

        end do

        ! Transform the density to reduced storage
        call switch_density(iLoc,DLT,Dab,kSym)

        !  ( r >= s )
        !---------------------------------------------------------
        ! V[i]{#J} <- V[i]{#J} + 2 * sum_rs  L(rs,{#J}) * D[i](rs)
        !=========================================================

        call DGEMV_('T',nRS,JNUM,Two,Lab,nRS,Dab,1,Zero,VJ,1)

        !-------------------------------------------------------
        ! ER[i] <- ER[i]  +  sum_J V[i](J)^2
        !=======================================================
        do jv=1,JNUM

          ER(iOcc(kSym)+ik) = ER(iOcc(kSym)+ik)+VJ(jv)**2
        end do

      end do

    end do

    call CWTIME(TCI2,TWI2)
    tintg(1) = tintg(1)+(TCI2-TCI1)
    tintg(2) = tintg(2)+(TWI2-TWI1)

  end do  !end batch loop

  ! free memory
  call mma_deallocate(VJ)
  call mma_deallocate(Lab)
  call mma_deallocate(Dab)

end do   ! loop over red sets

call mma_deallocate(DLT)

! Sync components
call GAdGOp(ER,nOccT,'+')

! Compute the ER-functional from its orbital components
W = sum(ER(1:nOccT))

call CWTIME(TOTCPU2,TOTWALL2)
TOTCPU = TOTCPU2-TOTCPU1
TOTWALL = TOTWALL2-TOTWALL1

! Write out timing information
if (timings) then

  write(u6,*)
  write(u6,*) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*) 'Timing from ',SECNAM,'            CPU      WALL '
  write(u6,*) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,'(2x,A26,2f10.2)') 'READ VECTORS                              ',tread(1),tread(2)
  write(u6,'(2x,A26,2f10.2)') 'COMPUTE (ii|ii)                           ',tintg(1),tintg(2)
  write(u6,'(2x,A26,2f10.2)') 'TOTAL                                     ',TOTCPU,TOTWALL
  write(u6,*) '- - - - - - - - - - - - - - - - - - - - - - - - -'
  write(u6,*)

end if

irc = 0

return

end subroutine CHO_get_ER
