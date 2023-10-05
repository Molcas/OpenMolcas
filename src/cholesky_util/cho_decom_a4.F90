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

subroutine Cho_Decom_A4(Diag,LstQSP,NumSP,iPass)
!
! Purpose: decompose qualified columns ("parallel" algorithm).

use Cholesky, only: Cho_Real_Par, INF_PASS, INF_PROGRESS, IPRINT, LQ_Tot, LQ, LuPri, LuSel, nnBstR, nQual, nSym, NumCho, NumChT, &
                    nVec_in_Buf, TDECOM
use Cholesky_procedures, only: Cho_P_GetLQ
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Diag(*)
integer(kind=iwp), intent(in) :: NumSP, LstQSP(NumSP), iPass
integer(kind=iwp) :: I, iEn, iK, iRed, iSt, iSym, iVec1, Jfi, Jin, jK, jVec, kI, kID, kK1, kK2, kK_1, kK_2, kQD, kV, l_IDKVec, &
                     l_KVec, l_LQ, l_Wrk1, l_xInt, LenLin, lK, MxQ, nkVec(8), nQual_Old(8), NumCho_Old(8), NumV(8)
real(kind=wp) :: C1, C2, W1, W2
integer(kind=iwp), allocatable :: IDKVec(:), iQScr(:)
real(kind=wp), allocatable :: KVec(:), KVScr(:), MQ(:), QDiag(:), Wrk1(:), xInt(:)
character(len=*), parameter :: SecNam = 'Cho_Decom_A4'

! Print header.
! -------------

LenLin = 0
if (iPrint >= Inf_Progress) then
  call Cho_Head(SecNam//': Decomposition of Qualified Diagonals','=',80,LUPRI)
  write(Lupri,'(/,A,I5,A,I4,A)') 'Integral pass number',iPass,' (',NumSP,' shell pair distributions calculated)'
  write(Lupri,'(A,8I8)') '#Cholesky vec.: ',(NumCho(iSym),iSym=1,nSym)
  write(Lupri,'(A,8I8)') '#vec. in buff.: ',(nVec_in_Buf(iSym),iSym=1,nSym)
  write(Lupri,'(A,8I8)') '#qualified    : ',(nQual(iSym),iSym=1,nSym)
  write(Lupri,'(A,8I8)') 'Current  dim. : ',(nnBstr(iSym,2),iSym=1,nSym)
  write(Lupri,'(A,8I8)') 'Original dim. : ',(nnBstr(iSym,1),iSym=1,nSym)
  write(Lupri,'(/,A,/,A)') '           #Vectors             Treated Diagonal', &
                           'Sym.     Sym.     Total     Index     Before      After   Conv. Neg.   New Max'
  LenLin = 79
  write(Lupri,'(80A)') ('-',I=1,LenLin)
  call XFlush(Lupri)
  NumCho_Old(1:nSym) = NumCho(1:nSym)
else if (iPrint >= Inf_Pass) then
  write(Lupri,'(/,A,I4)') 'Number of shell pair distributions calculated:',NumSP
  write(Lupri,'(A,8I8)') '#Cholesky vec.: ',(NumCho(iSym),iSym=1,nSym)
  write(Lupri,'(A,8I8)') '#vec. in buff.: ',(nVec_in_Buf(iSym),iSym=1,nSym)
  write(Lupri,'(A,8I8)') '#qualified    : ',(nQual(iSym),iSym=1,nSym)
  call XFlush(Lupri)
  NumCho_Old(1:nSym) = NumCho(1:nSym)
end if

! Allocations.
! ------------

call Cho_P_GetGV(numV,nSym)

l_KVec = nQual(1)**2
l_IDKVec = nQual(1)
l_LQ = nQual(1)*NumV(1)
do iSym=2,nSym
  l_KVec = l_KVec+nQual(iSym)**2
  l_IDKVec = l_IDKVec+nQual(iSym)
  l_LQ = l_LQ+nQual(iSym)*NumV(iSym)
end do
l_LQ = max(l_LQ,1) ! because there might not be any prev. vecs.
call mma_allocate(KVec,l_KVec,Label='KVec')
call mma_allocate(IDKVec,l_IDKVec,Label='IDKVec')
call mma_allocate(QDiag,l_IDKVec,Label='QDiag')

! Extract elements corresponding to qualified diagonals from
! previous Cholesky vectors (if any).
! ----------------------------------------------------------

call CWTime(C1,W1)
call mma_allocate(LQ_Tot,l_LQ,Label='LQ_Tot')

iEn = 0
iSt = 1
do iSym=1,nSym
  if (nQual(iSym)*NumV(iSym) > 0) then
    iEn = iEn+nQual(iSym)*NumV(iSym)
    LQ(iSym)%A(1:nQual(iSym),1:NumV(iSym)) => LQ_Tot(iSt:iEn)
    iSt = iEn+1
  else
    nullify(LQ(iSym)%A)
  end if
end do

call Cho_P_GetLQ(LQ_Tot,l_LQ,LstQSP,NumSP)
call CWTime(C2,W2)
tDecom(1,2) = tDecom(1,2)+C2-C1
tDecom(2,2) = tDecom(2,2)+W2-W1

! Extract qualified diagonal integral block.
! ------------------------------------------

call CWTime(C1,W1)

call mma_allocate(MQ,l_KVec,Label='MQ')
call Cho_P_GetMQ(MQ,size(MQ),LstQSP,NumSP)

! Decompose qualified diagonal block.
! The qualified diagonals are returned in QDiag.
! ----------------------------------------------

call Cho_Dec_Qual(Diag,LQ_Tot,MQ,KVec,IDKVec,nKVec,QDiag)

! Deallocate MQ.
! --------------

call mma_deallocate(MQ)

! Reorder the elements of the K-vectors according to IDK ordering.
! ----------------------------------------------------------------

MxQ = nQual(1)
do iSym=2,nSym
  MxQ = max(MxQ,nQual(iSym))
end do

call mma_allocate(KVScr,MxQ,Label='KVScr')

kK1 = 0
kK2 = kK1
kID = 0
do iSym=1,nSym
  do iK=1,nKVec(iSym)
    kK_1 = kK1+nQual(iSym)*(iK-1)
    KVScr(1:nQual(iSym)) = KVec(kK_1+1:kK_1+nQual(iSym))
    kK_2 = kK2+nKVec(iSym)*(iK-1)
    do jK=1,nKVec(iSym)
      lK = IDKVec(kID+jK)
      KVec(kK_2+jK) = KVScr(lK)
    end do
  end do
  kK1 = kK1+nQual(iSym)**2
  kK2 = kK2+nKVec(iSym)**2
  kID = kID+nQual(iSym)
end do

! Reorder QDiag to IDK ordering.
! ------------------------------

kID = 0
kQD = 0
do iSym=1,nSym
  KVScr(1:nQual(iSym)) = QDiag(kQD+1:kQD+nQual(iSym))
  do iK=1,nKVec(iSym)
    lK = IDKVec(kID+iK)
    QDiag(kQD+iK) = KVScr(lK)
  end do
  kQD = kQD+nQual(iSym)
  kID = kID+nQual(iSym)
end do

! Reorder elements of LQ vectors to IDK ordering.
! -----------------------------------------------

kID = 0
do iSym=1,nSym
  if (nQual(iSym) < 1) cycle
  do jVec=1,NumV(iSym)
    KVScr(1:nQual(iSym)) = LQ(iSym)%A(:,jVec)
    do iK=1,nKVec(iSym)
      lK = IDKVec(kID+iK)
      LQ(iSym)%A(iK,jVec) = KVScr(lK)
    end do
  end do
  kID = kID+nQual(iSym)
end do

call mma_deallocate(KVScr)

! Reset qualification index arrays to IDK ordering.
! Local as well as global are reordered.
! -------------------------------------------------

nQual_Old(1:nSym) = nQual(1:nSym)
call mma_allocate(iQScr,MxQ,Label='iQScr')

call Cho_P_ReoQual(iQScr,IDKVec,nKVec)

call mma_deallocate(iQScr)
nQual(1:nSym) = nKVec(1:nSym)

call CWTime(C2,W2)
tDecom(1,4) = tDecom(1,4)+C2-C1
tDecom(2,4) = tDecom(2,4)+W2-W1

! Compute vectors in each symmetry block.
! ---------------------------------------

kV = 1
kI = 1
kQD = 1
do iSym=1,nSym

  ! Cycle loop if nothing to do in this symmetry.
  ! ---------------------------------------------

  if (nQual(iSym) >= 1) then

    ! Set vector information.
    ! -----------------------

    call Cho_P_SetVecInf(nQual(iSym),iSym,iPass)

    ! Allocate memory for integrals/vectors.
    ! --------------------------------------

    l_xInt = max(nnBstR(iSym,2)*nQual(iSym),1)
    call mma_allocate(xInt,l_xInt,Label='xInt')

    if (nnBstR(iSym,2) > 0) then

      ! Read integral columns from disk, ordered according to IDK.
      ! ----------------------------------------------------------

      call CWTime(C1,W1)
      call Cho_RdQCol_Indx(xInt,IDKVec(kI),nnBstR(iSym,2),nQual(iSym),LuSel(iSym))
      call CWTime(C2,W2)
      tDecom(1,1) = tDecom(1,1)+C2-C1
      tDecom(2,1) = tDecom(2,1)+W2-W1

      ! Compute vectors.
      ! ----------------

      call mma_maxDBLE(l_Wrk1)
      call mma_allocate(Wrk1,l_Wrk1,Label='Wrk1')

      call Cho_CompVec(Diag,xInt,KVec(kV),QDiag(kQD),Wrk1,l_Wrk1,iSym,iPass)

      call mma_deallocate(Wrk1)

      ! Write vectors to disk and update vector counters.
      ! -------------------------------------------------

      call CWTime(C1,W1)
      iVec1 = NumCho(iSym)+1
      call Cho_PutVec(xInt,nnBstR(iSym,2),nQual(iSym),iVec1,iSym)
      call Cho_VecBuf_Copy(xInt,nQual(iSym),iSym)
      NumCho(iSym) = NumCho(iSym)+nQual(iSym)
      NumChT = NumChT+nQual(iSym)
      call CWTime(C2,W2)
      tDecom(1,2) = tDecom(1,2)+C2-C1
      tDecom(2,2) = tDecom(2,2)+W2-W1

    end if

    ! Transpose vectors on disk (parallel only).
    ! ------------------------------------------

    call CWTime(C1,W1)
    iRed = 2
    Jin = NumV(iSym)+1
    Jfi = NumV(iSym)+nQual(iSym)
    if (Cho_Real_Par) call Cho_VecTransp(xInt,Jin,Jfi,iSym,iRed,iPass)
    call CWTime(C2,W2)
    tDecom(1,2) = tDecom(1,2)+C2-C1
    tDecom(2,2) = tDecom(2,2)+W2-W1

    ! Deallocate memory for integrals/vectors.
    ! ----------------------------------------

    call mma_deallocate(xInt)

  end if

  ! Empty symmetry blocks jump here.
  ! --------------------------------

  kV = kV+nQual(iSym)**2
  kI = kI+nQual_Old(iSym)
  kQD = kQD+nQual_Old(iSym)

end do

! Deallocations.
! --------------

call mma_deallocate(LQ_Tot)
call mma_deallocate(QDiag)
call mma_deallocate(IDKVec)
call mma_deallocate(KVec)

! Print.
! ------

if (iPrint >= Inf_Progress) then
  NumCho_Old(1:nSym) = NumCho(1:nSym)-NumCho_Old(1:nSym)
  write(Lupri,'(80A)') ('-',I=1,LenLin)
  write(Lupri,'(A,8I8)') '#vec. gener.  : ',(NumCho_OLD(iSym),iSym=1,nSym)
else if (iPrint >= Inf_Pass) then
  NumCho_Old(1:nSym) = NumCho(1:nSym)-NumCho_Old(1:nSym)
  write(Lupri,'(A,8I8)') '#vec. gener.  : ',(NumCho_OLD(iSym),iSym=1,nSym)
end if

end subroutine Cho_Decom_A4
