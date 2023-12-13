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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_X_CompVec(irc,NVT,l_NVT,nBlock,l_nBlock,nV,l_nV1,l_nV2,iV1,l_iV11,l_iV12,ip_Z,l_Z1,l_Z2,Z,l_Z,Free_Z)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: Compute Cholesky vectors from Z vectors in core
!          and integrals computed on the fly:
!
!          L(uv,J) = 1/Z(J,J)
!                  * [ (uv|J) - sum[K=1,J-1] L(uv,K)*Z(J,K) ]
!
! NVT(i): total number of Cholesky vectors (all nodes), sym. i
! nBlock(i): number of vector blocks, sym. i
! nV(i,j): number of vectors in block i of sym. j
! iV1(i,j): index of first vector in block i of sym. j
! ip_Z(i,j): pointer to Z block i of sym. j
!
! If (Free_Z): Z array will be de-allocated here using ip_Z(1,1) as
! start point of the array - i.e. assuming that the Z blocks are
! stored as one coherent array. If this is not the case, simply put
! Free_Z=.False. Letting this routine deallocate Z maximizes the
! memory available for vector distribution.
!
! Vectors are distributed across nodes and stored according to
! reduced set 1 (all of them!).

use Index_Functions, only: iTri, nTri_Elem
#ifdef _DEBUGPRINT_
use Cholesky, only: iSP2F
#endif
use Cholesky, only: IndRSh, INF_PASS, iOff_Batch, IPRINT, iQuAB, iQuAB_here, LuPri, LuTmp, MaxQual, nDGM_call, nDim_Batch, &
                    nnBstRSh, nnShl, nQual, nSym, pTemp, TDECOM
use Cholesky_procedures, only: Cho_X_GetIP_InfVec
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: l_NVT, NVT(l_NVT), l_nBlock, nBlock(l_nBlock), l_nV1, l_nV2, nV(l_NV1,l_NV2), l_iV11, l_iV12, &
                                 iV1(l_IV11,l_iV12), l_Z1, l_Z2, ip_Z(l_Z1,l_Z2), l_Z
real(kind=wp), intent(inout) :: Z(l_Z)
logical(kind=iwp) :: Free_Z
integer(kind=iwp) :: iAdr(8), iCountSP, incZd, iSp, iSP_, iSP_1, iSP_2, iSym, J, J_inBlock, jBlock, K, K_inBlock, kBlock, kI, kL, &
                     kOffI, kOffZ, kZ, kZd, l_Int, l_Wrk, l_Zd, ldL, ldZ, Left, lTot, MaxQual_SAVE, n, nBatch, nSP, nSP_Max, &
                     nSP_this_batch
real(kind=wp) :: Byte, C0, C1, PDone, PMem, TotCPU, TotMem, TotWall, W0, W1, X0, X1, Y0, Y1
character(len=2) :: Unt
integer(kind=iwp), pointer :: InfVcT(:,:,:)
integer(kind=iwp), allocatable :: XCVLSP(:), XCVnBt(:), XCVTMP(:)
real(kind=wp), allocatable :: XCVInt(:), XCVZd(:)
real(kind=wp), parameter :: Mnt = 60.0_wp
character(len=*), parameter :: SecNam = 'Cho_X_CompVec'
integer(kind=iwp), external :: Cho_F2SP
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, nBlock_Max, nnB, nTot, nTot2
integer(kind=iwp), parameter :: myDebugInfo = 100
real(kind=wp), parameter :: Tol = 1.0e-14_wp
#endif

! Init return code
irc = 0

#ifdef _DEBUGPRINT_
! Check input
if ((l_NVT < nSym) .or. (l_nBlock < nSym) .or. (l_nV2 < nSym) .or. (l_iV12 < nSym) .or. (l_Z2 < nSym)) then
  irc = -2
  return
end if
nBlock_Max = nBlock(1)
do iSym=2,nSym
  nBlock_Max = max(nBlock_Max,nBlock(iSym))
end do
nnB = nTri_Elem(nBlock_Max)
if ((l_nV1 < nBlock_Max) .or. (l_iV11 < nBlock_Max) .or. (l_Z1 < nnB)) then
  irc = -1
  return
end if
#endif

! Parallel runs: open tmp vector files
call Cho_XCV_TmpFiles(irc,1)
if (irc /= 0) then
  write(LuPri,'(A,A,I8,A)') SecNam,': [1] Error in Cho_XCV_TmpFiles! (Return code:',irc,')'
  irc = 1
  return
end if

! Get pointer to InfVec array for all vectors
call Cho_X_GetIP_InfVec(InfVcT)

! Copy reduced set 1 to location 2.
! I.e. make rs1 the "current" reduced set.
call Cho_X_RSCopy(irc,1,2)
if (irc /= 0) then
  write(LuPri,'(A,A,I8,A)') SecNam,': Error in Cho_X_RSCopy! (Return code:',irc,')'
  irc = 1
  return
end if

! Allocate tmp array to keep track of which shell pairs are
! computed on this node.
call mma_allocate(XCVTMP,nnShl,Label='XCVTMP')

! Allocate and set list of parent shell pairs to compute
! NOTE: the tmp array allocated above is used here!!!
XCVTMP(:) = 0
nSP = 0
do iSym=1,nSym
  do J=1,NVT(iSym)
    iSP = Cho_F2SP(IndRSh(InfVcT(J,1,iSym))) ! Parent SP
    if ((iSP > 0) .and. (iSP <= nnShl)) then
      if (XCVTMP(iSP) == 0) nSP = nSP+1
      XCVTMP(iSP) = XCVTMP(iSP)+1
    else
      call Cho_Quit(SecNam//': Parent SP error!!',103)
    end if
  end do
end do
call mma_allocate(XCVLSP,nSP,Label='XCVLSP')
nSP = 0
do iSP=1,nnShl
  if (XCVTMP(iSP) > 0) then
    XCVLSP(nSP) = iSP
    nSP = nSP+1
  end if
end do
#ifdef _DEBUGPRINT_
if (nSP /= size(XCVLSP)) call Cho_Quit('SP counting error [1] in '//SecNam,103)
nTot = sum(NVT(1:nSym))
nTot2 = 0
do i=1,nSP
  nTot2 = nTot2+XCVTMP(XCVLSP(i))
end do
if (iPrint >= myDebugInfo) then
  call Cho_Head('Shell pair info from '//SecNam,'*',80,LuPri)
  write(LuPri,'(/,A,I8,/)') 'Number of shell pairs giving rise to vectors:',nSP
  do i=1,nSP
    iSP = XCVLSP(i)
    call Cho_InvPck(iSP2F(iSP),j,k,.true.)
    write(Lupri,'(A,I8,A,I4,A,I4,A,I8,A)') 'Shell pair',iSP,' (shells',j,',',k,') gives rise to',XCVTMP(iSP),' vectors'
  end do
  write(LuPri,'(A,I8,/,A,I8,/,A,I8)') 'Total number of vectors (from SP).................',nTot2, &
                                      'Total number of vectors (from NVT)................',nTot, &
                                      'Difference........................................',nTot2-nTot
end if
if (nTot /= nTot2) call Cho_Quit('SP counting error [2] in '//SecNam,103)
#endif

! Re-allocate and set qualification arrays
pTemp => iQuAB
MaxQual_SAVE = MaxQual
MaxQual = NVT(1)
do iSym=2,nSym
  MaxQual = max(MaxQual,NVT(iSym))
end do
call mma_allocate(iQuAB_here,MaxQual,nSym,Label='iQuAB_here')
iQuAB => iQuAB_here
nQual(1:nSym) = NVT(1:nSym)
iQuAB(1:nQual(iSym),1:nSym) = InfVcT(1:nQual(iSym),1,1:nSym) ! parent product for vector J

! Modify diagonal of Z matrix
! Z(J,J) <- 1/Z(J,J)
! If Z memory should not be de-allocated, save a backup of the
! diagonal.
if (Free_Z) then
  l_Zd = 1
  incZd = 0
else
  l_Zd = sum(NVT(1:nSym))
  incZd = 1
end if
call mma_allocate(XCVZd,l_Zd,Label='XCVZd')
kZd = 1
do iSym=1,nSym
  do jBlock=1,nBlock(iSym)
    kOffZ = ip_Z(iTri(jBlock,jBlock),iSym)-1
    do J_inBlock=1,nV(jBlock,iSym)
#     ifdef _DEBUGPRINT_
      ! Check for division by zero or negative diagonal
      ! This would be a bug....
      if ((abs(Z(kOffZ+iTri(J_inBlock,J_inBlock))) < Tol) .or. (Z(kOffZ+iTri(J_inBlock,J_inBlock)) < -Tol)) then
        write(LuPri,'(//,A,A)') SecNam,': Ooooops!'
        write(LuPri,'(A)') '....division by small or negative number....'
        write(LuPri,'(A,I8)') 'iSym=',iSym
        write(LuPri,'(A,I8)') 'jBlock=',jBlock
        write(LuPri,'(A,I8)') 'J=',iV1(jBlock,iSym)+J_inBlock-1
        write(LuPri,'(A,ES15.6)') 'Z(J,J)=',Z(kOffZ+iTri(J_inBlock,J_inBlock))
        write(LuPri,'(A,ES15.6)') 'Tolerance=',Tol
        call Cho_Quit('Division by small or negative number in '//SecNam,103)
      end if
#     endif
      XCVZd(kZd) = Z(kOffZ+iTri(J_inBlock,J_inBlock))
      Z(kOffZ+iTri(J_inBlock,J_inBlock)) = One/XCVZd(kZd)
      kZd = kZd+incZd
    end do
  end do
end do

! Distribute shell pairs across nodes according to dimension
! (helps load balance the linear algebra part below)
call Cho_XCV_Distrib_SP(XCVTMP,size(XCVTMP),iCountSP)

! Allocate batch dimension array
call mma_allocate(XCVnBt,max(iCountSP,1),Label='XCVnBt')

! Allocate offset array for batched SP loop
call mma_allocate(iOff_Batch,nSym,nnShl,Label='iOff_Batch')

! Split remaining memory in two parts.
! One for integrals/vectors, one for Seward.
call mma_maxDBLE(l_Wrk)
if (l_Wrk < 3) call Cho_Quit('Insufficient memory in '//SecNam,101)
l_Int = int(real(l_Wrk,kind=wp)/OneHalf)
call mma_allocate(XCVInt,l_Int,Label='XCVInt')
call mma_MaxDBLE(l_Wrk)
call xSetMem_Ints(l_Wrk)

! Print
if (iPrint >= Inf_Pass) then
  call Cho_Head('Generation of Cholesky vectors','=',80,LuPri)
  call Cho_Word2Byte(l_Int,8,Byte,Unt)
  write(LuPri,'(/,A,I12,A,F10.3,1X,A,A)') 'Memory available for integrals/vectors..',l_Int,' (',Byte,Unt,')'
  call Cho_Word2Byte(l_Wrk,8,Byte,Unt)
  write(LuPri,'(A,I12,A,F10.3,1X,A,A)') 'Memory available for Seward.............',l_Wrk,' (',Byte,Unt,')'
  write(LuPri,'(A,I12)') 'Number of shell pairs, total (nnShl)....',nnShl
  write(LuPri,'(A,I12)') 'Number of shell pairs, this node........',iCountSP
  write(LuPri,'(//,65X,A)') 'Time/min'
  write(LuPri,'(1X,A,5X,A,5X,A,2X,A,10X,A,16X,A,6X,A)') 'Batch','iSP1','iSP2','%Done','Memory','CPU','Wall'
  write(LuPri,'(A)') '----------------------------------------------------------------------------'
  call XFlush(LuPri)
  TotMem = Zero
  TotCPU = Zero
  TotWall = Zero
end if

! Compute Cholesky vectors in batched loop over shell pairs
iAdr(1:nSym) = 0 ! disk addresses
iSP_1 = 1
nBatch = 0
do while (iSP_1 <= iCountSP)
  if (iPrint >= Inf_Pass) call CWTime(X0,Y0)
  ! Set batch info
  Left = l_Int
  nSP_Max = iCountSP-iSP_1+1
  nSP_this_batch = 0
  do while ((Left > 0) .and. (nSP_this_batch < nSP_Max))
    iSP = iSP_1+nSP_this_batch
    n = sum(nnBstRSh(1:nSym,XCVTMP(iSP),2)*NVT(1:nSym))
    if (n <= Left) then
      Left = Left-n
      nSP_this_batch = nSP_this_batch+1
    else
      Left = 0 ! break while loop
    end if
  end do
  if (nSP_this_batch < 1) call Cho_Quit(SecNam//': Insufficient memory for shell pair batch',101)
  iSP_2 = iSP_1+nSP_this_batch-1
  do iSym=1,nSym
    nDim_Batch(iSym) = 0
    do iSP_=iSP_1,iSP_2
      iSP = XCVTMP(iSP_)
      iOff_Batch(iSym,iSP) = nDim_Batch(iSym)
      nDim_Batch(iSym) = nDim_Batch(iSym)+nnBstRSh(iSym,iSP,2)
    end do
  end do
  ! Calculate integrals (uv|J) for uv in shell pair batch and
  ! for all J (shell pairs that give rise to vectors are
  ! listed in ListSP).
  call Cho_XCV_GetInt(irc,XCVTMP(iSP_1),nSP_this_batch,XCVLSP,size(XCVLSP),NVT,l_NVT,XCVInt,size(XCVInt))
  if (irc /= 0) then
    write(LuPri,'(A,A,I8)') SecNam,': Cho_XCV_GetInt returned code',irc
    call Cho_Quit(SecNam//': Error in Cho_XCV_GetInt',104)
  end if
  ! Convert integrals into Cholesky vectors in each symmetry
  call CWTime(C0,W0)
  kOffI = 1
  do iSym=1,nSym
    ldL = max(nDim_Batch(iSym),1)
    do jBlock=1,nBlock(iSym)
      kOffZ = ip_Z(iTri(jBlock,jBlock),iSym)-1
      do J_inBlock=1,nV(jBlock,iSym)
        ! Convert integral column into Cholesky vector
        ! L(uv,J)=(uv|J)/Z(J,J)
        ! Note that the inverse was taken before the integral
        ! loop - i.e. Z(J,J) <- 1/Z(J,J)
        J = iV1(jBlock,iSym)+J_inBlock-1
        kL = kOffI+nDim_Batch(iSym)*(J-1)
        XCVInt(kL:kL+nDim_Batch(iSym)-1) = Z(kOffZ+iTri(J_inBlock,J_inBlock))*XCVInt(kL:kL+nDim_Batch(iSym)-1)
        ! Subtract from subsequent columns in current block
        ! using BLAS1
        ! (uv|K) <- (uv|K) - L(uv,J)*Z(K,J), K>J (in jBlock)
        do K_inBlock=J_inBlock+1,nV(jBlock,iSym)
          K = iV1(jBlock,iSym)+K_inBlock-1
          XCVInt(kOffI+nDim_Batch(iSym)*(K-1):kOffI+nDim_Batch(iSym)*K-1) = &
            XCVInt(kOffI+nDim_Batch(iSym)*(K-1):kOffI+nDim_Batch(iSym)*K-1)- &
            Z(kOffZ+iTri(K_inBlock,J_inBlock))*XCVInt(kL:kL+nDim_Batch(iSym)-1)
        end do
      end do
      ! Subtract from subsequent blocks using BLAS3
      ! (uv|K) <- (uv|K) - sum[J] L(uv,J)*Z(K,J),
      ! K>J (K in kBlock>jBlock containing J)
      kL = kOffI+nDim_Batch(iSym)*(iV1(jBlock,iSym)-1)
      do kBlock=jBlock+1,nBlock(iSym)
        kI = kOffI+nDim_Batch(iSym)*(iV1(kBlock,iSym)-1)
        kZ = ip_Z(iTri(kBlock,jBlock),iSym)
        ldZ = max(nV(kBlock,iSym),1)
        call dGeMM_('N','T',nDim_Batch(iSym),nV(kBlock,iSym),nV(jBlock,iSym),-One,XCVInt(kL),ldL,Z(kZ),ldZ,One,XCVInt(kI),ldL)
      end do
    end do
    nDGM_Call = nDGM_Call+nTri_Elem(nBlock(iSym)-1)
    kOffI = kOffI+nDim_Batch(iSym)*NVT(iSym)
  end do
  call CWTime(C1,W1)
  tDecom(1,3) = tDecom(1,3)+(C1-C0)
  tDecom(2,3) = tDecom(2,3)+(W1-W0)
  ! Write vectors to temp files
  call CWTime(C0,W0)
  kL = 1
  do iSym=1,nSym
    lTot = nDim_Batch(iSym)*NVT(iSym)
    if (lTot > 0) then
      call DDAFile(LuTmp(iSym),1,XCVInt(kL),lTot,iAdr(iSym))
      kL = kL+lTot
    end if
  end do
  call CWTime(C1,W1)
  tDecom(1,2) = tDecom(1,2)+(C1-C0)
  tDecom(2,2) = tDecom(2,2)+(W1-W0)
  ! Print
  if (iPrint >= Inf_Pass) then
    call CWTime(X1,Y1)
    PDone = 1.0e2_wp*real(iSP_2,kind=wp)/real(iCountSP,kind=wp)
    lTot = sum(nDim_Batch(1:nSym)*NVT(1:nSym))
    call Cho_Word2Byte(lTot,8,Byte,Unt)
    PMem = 1.0e2_wp*real(lTot,kind=wp)/real(l_Int,kind=wp)
    write(LuPri,'(I6,1X,I8,1X,I8,1X,F6.1,1X,F10.3,1X,A,A,F7.2,A,1X,F9.2,1X,F9.2)') nBatch+1,iSP_1,iSP_2,PDone,Byte,Unt,' (',PMem, &
                                                                                   '%)',(X1-X0)/Mnt,(Y1-Y0)/Mnt
    call XFlush(LuPri)
    TotMem = TotMem+real(lTot,kind=wp)
    TotCPU = TotCPU+(X1-X0)/Mnt
    TotWall = TotWall+(Y1-Y0)/Mnt
  end if
  ! Update counters and save batch dimension
  iSP_1 = iSP_1+nSP_this_batch
  nBatch = nBatch+1
  XCVnBt(nBatch) = nSP_this_batch
end do
if (iPrint >= Inf_Pass) then
  write(LuPri,'(A)') '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - '
  call Cho_RWord2Byte(TotMem,Byte,Unt)
  write(LuPri,'(32X,F10.3,1X,A,12X,F9.2,1X,F9.2)') Byte,Unt,TotCPU,TotWall
  write(LuPri,'(A)') '----------------------------------------------------------------------------'
  call XFlush(LuPri)
end if

! If not deallocate Z matrix then reconstruct diagonal of Z matrix
! Z(J,J) <- 1/Z(J,J)
if (.not. Free_Z) then
  kZd = 0
  do iSym=1,nSym
    do jBlock=1,nBlock(iSym)
      kOffZ = ip_Z(iTri(jBlock,jBlock),iSym)-1
      do J_inBlock=1,nV(jBlock,iSym)
        J = iV1(jBlock,iSym)+J_inBlock-1
        Z(kOffZ+iTri(J_inBlock,J_inBlock)) = XCVZd(kZd+J)
      end do
    end do
    kZd = kZd+NVT(iSym)
  end do
end if

! Deallocations
call xRlsMem_Ints()
call mma_deallocate(XCVInt)
call mma_deallocate(iOff_Batch)
call mma_deallocate(XCVZd)
nullify(iQuAB)
call mma_deallocate(iQuAB_here)
call mma_deallocate(XCVLSP)

! Reset qualification array pointers and MaxQual
iQuAB => pTemp
MaxQual = MaxQual_SAVE

! Parallel runs: distribute vectors across nodes (store on files)
! Serial runs: write vectors to permanent files
call CWTime(C0,W0)
call Cho_XCV_DistributeVectors(irc,XCVnBt,nBatch,XCVTMP,iCountSP,NVT,l_NVT)
if (irc /= 0) then
  write(LuPri,'(A,A,I8)') SecNam,': Cho_XCV_DistributeVectors returned code',irc
  call Cho_Quit(SecNam//': Error in Cho_XCV_DistributeVectors',104)
end if
call CWTime(C1,W1)
tDecom(1,2) = tDecom(1,2)+(C1-C0)
tDecom(2,2) = tDecom(2,2)+(W1-W0)

! Deallocations
call mma_deallocate(XCVnBt)
call mma_deallocate(XCVTMP)

! Parallel runs: close and erase tmp vector files
call Cho_XCV_TmpFiles(irc,3)
if (irc /= 0) then
  write(LuPri,'(A,A,I8,A)') SecNam,': [3] Error in Cho_XCV_TmpFiles! (Return code:',irc,')'
  call Cho_Quit(SecNam//': Error in Cho_XCV_TmpFiles',104)
end if

end subroutine Cho_X_CompVec
