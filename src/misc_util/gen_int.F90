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
! Copyright (C) 2004, Francesco Aquilante                              *
!***********************************************************************
!  GEN_INT
!
!> @brief
!>   Generates integrals from Cholesky vectors
!> @author F. Aquilante, May 2004
!> @modified_by F. Aquilante, Sep. 2004
!>
!> @note
!> The transposition ``L(ab,J)`` &rarr; ``L(ba,J)`` of the vectors
!> ``(syma /= symb)`` is necessary because the calling routine
!> requires the integrals in the order \f$ (sr|qp) \f$ which
!> is reversed compared to the order of the symmetries
!> given as input arguments.
!>
!> @param[out] rc
!> @param[in]  iSymp
!> @param[in]  iSymq
!> @param[in]  iSymr
!> @param[in]  iSyms
!> @param[in]  ipq1
!> @param[in]  numpq
!> @param[out] Xint
!***********************************************************************

subroutine GEN_INT(rc,iSymp,iSymq,iSymr,iSyms,ipq1,numpq,Xint)
!***********************************************************************
!
! Modified  September 2004
! Reason:
! the transposition L(ab,J) --> L(ba,J) of the vectors
! (syma /= symb) is necessary because the calling routine
! requires the integrals in the order (sr|qp) which
! is reversed compared to the order of the symmetries
! given as input arguments
!
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use GetInt_mod, only: LuCVec, nBas, NumCho, pq1
use TwoDat, only: rcTwo
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: rc
integer(kind=iwp), intent(in) :: iSymp, iSymq, iSymr, iSyms, ipq1, numpq
real(kind=wp), intent(_OUT_) :: Xint(*)
integer(kind=iwp) :: iBatch, iVec1, J, jp, jpq, jq, jr, js, jSym, jvec, koff1, koff2, LWORK, mBatch, mNeed, Npq, Npqrs, Nrs, NumV, &
                     nVec, pq, pq1_save
real(kind=wp), allocatable :: Vec1(:), Vec2(:), Vec3(:)

jSym = Mul(iSymp,iSymq)

if (NumCho(jSym) < 1) return

! save the value of pq1 because it belongs to a Common block
pq1_save = pq1
pq1 = ipq1

if (iSymp == iSymq) then
  Npq = nTri_Elem(nBas(iSymp))
else
  Npq = nBas(iSymp)*nBas(iSymq)
end if

if (iSymr == iSyms) then
  Nrs = nTri_Elem(nBas(iSymr))
else
  Nrs = nBas(iSymr)*nBas(iSyms)
end if

if (iSymp == iSymr) then
  Npqrs = Npq
else
  Npqrs = max(Npq,Nrs)
end if

! Set up the batch procedure
! --------------------------
call mma_maxDBLE(LWORK)

! Memory management
if (iSymp /= iSymr) then
  mNeed = 2*max(Npq,Nrs)+Nrs
else
  mNeed = 2*Npq
end if

if (mNeed <= 0) then
  ! ***QUIT*** bad initialization
  write(u6,*) 'Gen_Int: bad initialization'
  rc = rcTwo%RD11
  call Abend()
end if
nVec = min(LWORK/mNeed,NumCho(jSym))
if (nVec <= 0) then
  ! ***QUIT*** insufficient memory
  write(u6,*) 'Gen_Int: Insufficient memory for batch'
  write(u6,*) 'LWORK= ',LWORK
  write(u6,*) 'mNeed= ',mNeed
  write(u6,*) 'NumCho= ',NumCho(jsym)
  write(u6,*) 'jsym= ',jsym
  rc = rcTwo%RD05
  call Abend()
end if
mBatch = (NumCho(jSym)-1)/nVec+1

! Start the batch procedure for reading the vectors and computing the integrals

Xint(1:numpq*Nrs) = Zero

! Allocate memory for reading the vectors and do the transposition
call mma_allocate(Vec1,Npqrs*nVec,label='MemC1')
call mma_allocate(Vec2,Npqrs*nVec,label='MemC2')
if (iSymp /= iSymr) call mma_allocate(Vec3,Nrs*nVec,label='MemC3')

do iBatch=1,mBatch

  if (iBatch == mBatch) then
    NumV = NumCho(jSym)-nVec*(mBatch-1)
  else
    NumV = nVec
  end if

  iVec1 = nVec*(iBatch-1)+1

  !--- Copying out (and transpose) the elements of the 1st vector ---!
  !------------- L(pq,J) ---> L(qp,J)  ------------------------------!
  !------------------------------------------------------------------!
  if (iSymp /= iSymq) then  ! transposition needed
    call RdChoVec(Vec1,Npq,NumV,iVec1,LuCVec(1))
    do jvec=1,NumV
      do jq=1,nBas(iSymq)
        do jp=1,nBas(iSymp)

          koff1 = Npq*(jvec-1)+nBas(iSymp)*(jq-1)+jp
          koff2 = Npq*(jvec-1)+nBas(iSymq)*(jp-1)+jq
          Vec2(koff2) = Vec1(koff1)

        end do
      end do
    end do
  else  ! no need to transpose "diagonal" vectors
    call RdChoVec(Vec2,Npq,NumV,iVec1,LuCVec(1))
  end if

  if (numpq == Npq) then

    if (iSymp /= iSymr) then !need to read the 2nd vector also

      if (iSymr /= iSyms) then
        call RdChoVec(Vec1,Nrs,NumV,iVec1,LuCVec(2))
        do jvec=1,NumV
          do js=1,nBas(iSyms)
            do jr=1,nBas(iSymr)

              koff1 = Nrs*(jvec-1)+nBas(iSymr)*(js-1)+jr
              koff2 = Nrs*(jvec-1)+nBas(iSyms)*(jr-1)+js
              Vec3(koff2) = Vec1(koff1)

            end do
          end do
        end do
      else
        call RdChoVec(Vec3,Nrs,NumV,iVec1,LuCVec(2))
      end if

      ! Computing the integrals (II|JJ)
      ! -------------------------------
      ! (sr|{qp}) <- (sr|{qp}) + sum_I L(sr,#I)* L({qp},#I)
      !======================================================
      call DGEMM_('N','T',Nrs,numpq,NumV,One,Vec3,Nrs,Vec2,numpq,One,Xint,Nrs)

    else ! isymp = isymr   (Npq=Nrs)

      ! Computing integrals of the type (II|II) and (IJ|IJ)

      call DGEMM_('N','T',Nrs,numpq,NumV,One,Vec2,Nrs,Vec2,numpq,One,Xint,Nrs)

    end if

  else  ! numpq /= Npq

    !--- Copying out the elements of the 1st vector ---!
    !--------------------------------------------------!
    do J=1,NumV
      do jpq=1,numpq
        pq = ipq1+jpq-1
        ! Address of the matrix element (pq,J) in the full matrix
        kOff1 = Npq*(J-1)+pq
        ! Address of the matrix element (pq,J) in the sub-block matrix
        kOff2 = numpq*(J-1)+jpq
        ! Copy out the elements of the sub-block matrix if not the full matrix
        Vec1(kOff2) = Vec2(kOff1)
      end do
    end do

    if (iSymp /= iSymr) then

      if (iSymr /= iSyms) then  !   L(rs,J) ---> L(sr,J)
        call RdChoVec(Vec2,Nrs,NumV,iVec1,LuCVec(2))
        do jvec=1,NumV
          do js=1,nBas(iSyms)
            do jr=1,nBas(iSymr)

              koff1 = Nrs*(jvec-1)+nBas(iSymr)*(js-1)+jr
              koff2 = Nrs*(jvec-1)+nBas(iSyms)*(jr-1)+js
              Vec3(koff2) = Vec2(koff1)

            end do
          end do
        end do
      else
        call RdChoVec(Vec3,Nrs,NumV,iVec1,LuCVec(2))
      end if

      ! Computing the integrals
      ! -------------------------------
      ! (sr|{qp}) <- (sr|{qp}) + sum_I L(sr,#I)* L({qp},#I)
      !======================================================
      call DGEMM_('N','T',Nrs,numpq,NumV,One,Vec3,Nrs,Vec1,numpq,One,Xint,Nrs)

    else

      call DGEMM_('N','T',Nrs,numpq,NumV,One,Vec2,Nrs,Vec1,numpq,One,Xint,Nrs)

    end if

  end if

end do  ! end of the batch procedure

! Free the memory
call mma_deallocate(Vec1)
call mma_deallocate(Vec2)
if (iSymp /= iSymr) call mma_deallocate(Vec3)

rc = rcTwo%good
pq1 = pq1_save

return

end subroutine GEN_INT
