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

implicit real*8(a-h,o-z)
integer rc
integer iSymp, iSymq, iSymr, iSyms
integer pq, numpq, pq1_save
real*8 Xint(*)
#include "RdOrd.fh"
#include "WrkSpc.fh"
#include "TwoDat.fh"
! Statement function
MulD2h(i,j) = ieor(i-1,j-1)+1

jSym = MulD2h(iSymp,iSymq)

if (NumCho(jSym) < 1) return

! save the value of pq1 because it belongs to a Common block
pq1_save = pq1
pq1 = ipq1

if (iSymp == iSymq) then
  Npq = nBas(iSymp)*(nBas(iSymp)+1)/2
else
  Npq = nBas(iSymp)*nBas(iSymq)
end if

if (iSymr == iSyms) then
  Nrs = nBas(iSymr)*(nBas(iSymr)+1)/2
else
  Nrs = nBas(iSymr)*nBas(iSyms)
end if

! Set up the batch procedure
! --------------------------
call GetMem('Maxmem','MAX ','REAL',KDUM,LWORK)

! Memory management
if (iSymp /= iSymr) then
  mNeed = 2*max(Npq,Nrs)+Nrs
else
  mNeed = 2*Npq
end if

if (mNeed > 0) then
  nVec = min(LWORK/mNeed,NumCho(jSym))
else
  ! ***QUIT*** bad initialization
  write(6,*) 'Gen_Int: bad initialization'
  rc = 99
  call Abend()
  nVec = -9999  ! dummy assignment - avoid compiler warnings
end if
if (nVec > 0) then
  mBatch = (NumCho(jSym)-1)/nVec+1
else
  ! ***QUIT*** insufficient memory
  write(6,*) 'Gen_Int: Insufficient memory for batch'
  write(6,*) 'LWORK= ',LWORK
  write(6,*) 'mNeed= ',mNeed
  write(6,*) 'NumCho= ',NumCho(jsym)
  write(6,*) 'jsym= ',jsym
  rc = rcRD05
  call Abend()
  mBatch = -9999  ! dummy assignment
end if

! Start the batch procedure for reading the vectors and computing the integrals

!call FZero(Xint,numpq*Nrs)
do i=1,numpq*Nrs
  Xint(i) = ZERO
end do

do iBatch=1,mBatch

  if (iBatch == mBatch) then
    NumV = NumCho(jSym)-nVec*(mBatch-1)
  else
    NumV = nVec
  end if

  iVec1 = nVec*(iBatch-1)+1

  if (iSymp /= iSymr) then

    LenMem1 = max(Npq,Nrs)*NumV
    LenMem2 = max(Npq,Nrs)*NumV
    LenMem3 = Nrs*NumV
    ! Allocate memory for reading the vectors and do the transposition
    call GetMem('MemC1','ALLO','REAL',kVec1,LenMem1)
    call GetMem('MemC2','ALLO','REAL',kVec2,LenMem2)
    call GetMem('MemC3','ALLO','REAL',kVec3,LenMem3)

  else

    LenMem1 = Npq*NumV  ! equal to LenMem2
    ! Allocate memory for reading the vectors and do the transposition
    call GetMem('MemC1','ALLO','REAL',kVec1,LenMem1)
    call GetMem('MemC2','ALLO','REAL',kVec2,LenMem1)

  end if

  !--- Copying out (and transpose) the elements of the 1st vector ---!
  !------------- L(pq,J) ---> L(qp,J)  ------------------------------!
  !------------------------------------------------------------------!
  if (iSymp /= iSymq) then  ! transposition needed
    call RdChoVec(Work(kVec1),Npq,NumV,iVec1,LuCVec(1))
    koff1 = 0
    koff2 = 0
    do jvec=1,NumV
      do jq=1,nBas(iSymq)
        do jp=1,nBas(iSymp)

          koff1 = kVec1+Npq*(jvec-1)+nBas(iSymp)*(jq-1)+(jp-1)
          koff2 = kVec2+Npq*(jvec-1)+nBas(iSymq)*(jp-1)+(jq-1)
          work(koff2) = work(koff1)

        end do
      end do
    end do
  else  ! no need to transpose "diagonal" vectors
    call RdChoVec(Work(kVec2),Npq,NumV,iVec1,LuCVec(1))
  end if
  kWqp = kVec2

  if (numpq == Npq) then

    if (iSymp /= iSymr) then !need to read the 2nd vector also

      if (iSymr /= iSyms) then
        call RdChoVec(Work(kVec1),Nrs,NumV,iVec1,LuCVec(2))
        koff1 = 0
        koff2 = 0
        do jvec=1,NumV
          do js=1,nBas(iSyms)
            do jr=1,nBas(iSymr)

              koff1 = kVec1+Nrs*(jvec-1)+nBas(iSymr)*(js-1)+(jr-1)
              koff2 = kVec3+Nrs*(jvec-1)+nBas(iSyms)*(jr-1)+(js-1)
              work(koff2) = work(koff1)

            end do
          end do
        end do
      else
        call RdChoVec(Work(kVec3),Nrs,NumV,iVec1,LuCVec(2))
      end if
      kWsr = kVec3


      ! Computing the integrals (II|JJ)
      ! -------------------------------
      ! (sr|{qp}) <- (sr|{qp}) + sum_I L(sr,#I)* L({qp},#I)
      !======================================================
      call DGEMM_('N','T',Nrs,numpq,NumV,ONE,Work(kWsr),Nrs,WORK(kWqp),numpq,ONE,Xint,Nrs)

    else ! isymp = isymr   (Npq=Nrs)

      ! Computing integrals of the type (II|II) and (IJ|IJ)

      call DGEMM_('N','T',Nrs,numpq,NumV,ONE,Work(kWqp),Nrs,WORK(kWqp),numpq,ONE,Xint,Nrs)

    end if

  else  ! numpq /= Npq

    !--- Copying out the elements of the 1st vector ---!
    !--------------------------------------------------!
    do J=1,NumV
      do jpq=1,numpq
        pq = pq1+jpq-1
        ! Address of the matrix element (pq,J) in the full matrix
        kOff1 = kVec2+Npq*(J-1)+(pq-1)
        ! Address of the matrix element (pq,J) in the sub-block matrix
        kOff2 = kVec1+numpq*(J-1)+(jpq-1)
        ! Copy out the elements of the sub-block matrix if not the full matrix
        Work(kOff2) = Work(kOff1)
      end do
    end do
    kXqp = kVec1
    kWsr = kVec2

    if (iSymp /= iSymr) then

      if (iSymr /= iSyms) then  !   L(rs,J) ---> L(sr,J)
        call RdChoVec(Work(kVec2),Nrs,NumV,iVec1,LuCVec(2))
        koff1 = 0
        koff2 = 0
        do jvec=1,NumV
          do js=1,nBas(iSyms)
            do jr=1,nBas(iSymr)

              koff1 = kVec2+Nrs*(jvec-1)+nBas(iSymr)*(js-1)+(jr-1)
              koff2 = kVec3+Nrs*(jvec-1)+nBas(iSyms)*(jr-1)+(js-1)
              work(koff2) = work(koff1)

            end do
          end do
        end do
      else
        call RdChoVec(Work(kVec3),Nrs,NumV,iVec1,LuCVec(2))
      end if
      kWsr = kVec3

    end if

    ! Computing the integrals
    ! -------------------------------
    ! (sr|{qp}) <- (sr|{qp}) + sum_I L(sr,#I)* L({qp},#I)
    !======================================================
    call DGEMM_('N','T',Nrs,numpq,NumV,ONE,Work(kWsr),Nrs,WORK(kXqp),numpq,ONE,Xint,Nrs)

  end if

  ! Free the memory
  if (iSymp /= iSymr) then
    call GetMem('MemC3','FREE','REAL',kVec3,LenMem3)
    call GetMem('MemC2','FREE','REAL',kVec2,LenMem2)
    call GetMem('MemC1','FREE','REAL',kVec1,LenMem1)
  else
    call GetMem('MemC2','FREE','REAL',kVec2,LenMem1)
    call GetMem('MemC1','FREE','REAL',kVec1,LenMem1)
  end if

end do  ! end of the batch procedure

rc = rc0000
pq1 = pq1_save

return

end subroutine GEN_INT
