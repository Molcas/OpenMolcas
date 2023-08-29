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
!               Thomas Bondo Pedersen                                  *
!***********************************************************************
!  Cho_X_CalcChoDiag
!
!> @brief
!>   Calculate integral diagonal from Cholesky vectors
!> @author Francesco Aquilante
!> @modified_by Thomas Bondo Pedersen
!>
!> @details
!> This routine calculates the integral diagonal from Cholesky
!> vectors,
!>
!> \f[ (ab|ab) = \sum_J L_{ab,J}^2 \quad (a,b: \text{AO-indices}) \f]
!>
!> The diagonal calculation is parallelized.
!> The diagonal is returned in first reduced set storage and must
!> be allocated before calling this routine.
!> Return code is ``0`` if successful execution.
!>
!> @param[out] rc   Return code
!> @param[out] Diag Array containing diagonal on exit
!***********************************************************************

subroutine Cho_X_CalcChoDiag(rc,Diag)

use Cholesky, only: iiBstR, IndRed, InfVec, IndRed, nDimRS, nnBstRT, nSym, NumCho
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: rc
real(kind=wp), intent(_OUT_) :: Diag(*)
integer(kind=iwp) :: iBatch, iLoc, irc, IREDC, IVEC2, iVrs, JNUM, JRED, JRED1, JRED2, jrs, jSym, JVEC, krs, LWORK, mrs, MUSED, &
                     nBatch, nRS, NUMV, nVec, nVrs
real(kind=wp), allocatable :: Lrs(:,:)
character(len=*), parameter :: SECNAM = 'Cho_X_CalcChoDiag'

Diag(1:nnBstRT(1)) = Zero

IREDC = -1  ! unknown reduced set

iLoc = 3 ! use scratch location in reduced index arrays

! *************** BIG LOOP OVER VECTORS SYMMETRY *******************
do jSym=1,nSym

  if (NumCho(jSym) < 1) cycle

  JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
  JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec
  do JRED=JRED1,JRED2

    call Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

    if (nVrs == 0) cycle  ! no vectors in that (jred,jsym)

    if (nVrs < 0) then
      write(u6,*) SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
      rc = 77
      return
    end if

    call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
    if (irc /= 0) then
      write(u6,*) SECNAM//'cho_X_setred non-zero return code.  rc= ',irc
      rc = irc
      return
    end if

    IREDC = JRED

    nRS = nDimRS(JSYM,JRED)

    call mma_maxDBLE(LWORK)

    nVec = min(LWORK/max(nRS,1),nVrs)

    if (nVec < 1) then
      write(u6,*) SECNAM//': Insufficient memory for batch'
      write(u6,*) ' LWORK= ',LWORK
      write(u6,*) ' jsym= ',jsym
      write(u6,*) ' min. mem. need for reading= ',nRS
      rc = 33
      return
      nBatch = -9999  ! dummy assignment
    end if

    call mma_allocate(Lrs,nRS,nVec,Label='Lrs')

    ! BATCH over the vectors ----------------------------

    nBatch = (nVrs-1)/nVec+1

    do iBatch=1,nBatch

      if (iBatch == nBatch) then
        JNUM = nVrs-nVec*(nBatch-1)
      else
        JNUM = nVec
      end if

      JVEC = nVec*(iBatch-1)+iVrs
      IVEC2 = JVEC-1+JNUM

      call CHO_VECRD(Lrs,size(Lrs),JVEC,IVEC2,JSYM,NUMV,IREDC,MUSED)

      if ((NUMV <= 0) .or. (NUMV /= JNUM)) then
        call mma_deallocate(Lrs)
        rc = 77
        return
      end if

      ! ----------------------------------------------------------------
      ! Compute the diagonals :   D(ab) = D(ab) + sum_J (Lab,J)^2
      !
      ! Stored in the 1st reduced set

      do krs=1,nRS

        mrs = iiBstR(JSYM,iLoc)+krs
        jrs = IndRed(mrs,iLoc) ! address in 1st red set

        Diag(jrs) = Diag(jrs)+sum(Lrs(krs,1:JNUM)**2)

      end do

      ! ----------------------------------------------------------------
      ! ----------------------------------------------------------------

    end do  ! end batch loop

    ! free memory
    call mma_deallocate(Lrs)

  end do   ! loop over red sets

end do   !loop over JSYM

call Cho_GAdGOp(Diag(1),NNBSTRT(1),'+')

rc = 0

end subroutine Cho_X_CalcChoDiag
