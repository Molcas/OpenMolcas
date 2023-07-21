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
!               2012,2014, Victor P. Vysotskiy                         *
!***********************************************************************

subroutine Cho_XCV_DV_S(irc,SP_BatchDim,nSP_Batch,id_mySP,n_mySP)

use ChoSwp, only: nnBstRSh, iiBstRSh
use stdalloc

implicit none
integer irc
integer nSP_Batch
integer SP_BatchDim(nSP_Batch)
integer n_mySP
integer id_mySP(n_mySP)
#include "cholesky.fh"
#include "choprint.fh"
character*12 SecNam
parameter(SecNam='Cho_XCV_DV_S')
real*8 X0, X1, Y0, Y1
integer iSym, iSP, iSP1, iSP2, J1, nDim
integer ip_Mem, l_Mem
integer ip_V, ip_T
integer max_vector_dim, max_block_dim
integer nVec_per_batch, nVec_this_batch
integer iBatch, nBatch
integer kV, kT, kOffV, kOffT
integer lTot, iAdr, iAdr0
integer iSP_Batch, nSP_this_batch
integer i, j
real*8, allocatable :: DVSVEC(:)

! Init return code
irc = 0

! Find max vector dimension
max_vector_dim = nnBstR(1,2)
do iSym=2,nSym
  max_vector_dim = max(max_vector_dim,nnBstR(iSym,2))
end do
if (max_vector_dim < 1) then
  irc = -2
  return
end if

! Find max block dimension
max_block_dim = 0
iSP1 = 1
do iSP_Batch=1,nSP_Batch
  nSP_this_batch = SP_BatchDim(iSP_Batch)
  iSP2 = iSP1+nSP_this_batch-1
  do iSym=1,nSym
    nDim = 0
    do iSP=iSP1,iSP2
      nDim = nDim+nnBstRSh(iSym,id_mySP(iSP),2)
    end do
    max_block_dim = max(max_block_dim,nDim)
  end do
  iSP1 = iSP1+nSP_this_batch
end do
#ifdef _DEBUGPRINT_
if ((iSP1-1) /= n_mySP) call Cho_Quit(SecNam//': SP batch dimension error',103)
#endif

! Get largest memory block
call mma_maxDBLE(l_Mem)
lTot = max_vector_dim+max_block_dim
if (l_Mem < lTot) then
  irc = -1
  return
end if
call mma_allocate(DVSVEC,l_Mem,Label='DVSVEC')
ip_Mem = 1

! read, reorder, and write vectors, one symmetry at a time
do iSym=1,nSym
  if (iPrint >= Inf_Progress) then
    write(LuPri,'(/,A,I2,/,A)') 'Writing vectors, symmetry',iSym,'---------------------------'
    write(LuPri,'(3X,A,I8)') 'Total number of vectors:',NumCho(iSym)
    write(LuPri,'(3X,A,I8)') 'Vector dimension       :',nnBstR(iSym,2)
    write(LuPri,'(3X,A,I8)') 'Shell pair batches     :',nSP_Batch
    call Cho_Flush(LuPri)
  end if
  if ((NumCho(iSym) > 0) .and. (nnBstR(iSym,2) > 0)) then
    ! Set up batching
    lTot = max_block_dim+nnBstR(iSym,2)
    nVec_per_batch = min(l_Mem/lTot,NumCho(iSym))
    if (nVec_per_batch < 1) call Cho_Quit('Insufficient memory for batching in '//SecNam,101)
    nBatch = (NumCho(iSym)-1)/nVec_per_batch+1
    if (iPrint >= Inf_Progress) then
      write(LuPri,'(3X,A,I8)') 'Vector batches         :',nBatch
      call Cho_Flush(LuPri)
    end if
    ! Read and write vectors in batches
    do iBatch=1,nBatch
      ! Determine number of vectors in this batch
      if (iBatch == nBatch) then
        nVec_this_batch = NumCho(iSym)-nVec_per_batch*(nBatch-1)
      else
        nVec_this_batch = nVec_per_batch
      end if
      ! First vector in this batch
      J1 = nVec_per_batch*(iBatch-1)+1
      if (iPrint >= Inf_Progress) then
        write(LuPri,'(3X,A,I8,/,3X,A)') 'Vector batch number:',iBatch,'++++++++++++++++++++++++++++'
        write(LuPri,'(6X,A,I8)') 'Number of vectors in this batch:',nVec_this_batch
        write(LuPri,'(6X,A,I8,1X,I8)') 'First and last vector          :',J1,J1+nVec_this_batch-1
        call Cho_Flush(LuPri)
        call Cho_Timer(X0,Y0)
      end if
      ! Set memory pointers
      ip_T = ip_Mem
      ip_V = ip_T+max_block_dim*nVec_this_batch
      ! Read blocked vectors and copy into full vector array
      iAdr0 = 0
      iSP1 = 1
      do iSP_Batch=1,nSP_Batch
        nSP_this_batch = SP_BatchDim(iSP_Batch)
        iSP2 = iSP1+nSP_this_batch-1
        nDim = 0
        do iSP=iSP1,iSP2
          nDim = nDim+nnBstRSh(iSym,id_mySP(iSP),2)
        end do
        if (nDim > 0) then
          lToT = nDim*nVec_this_batch
          iAdr = iAdr0+nDim*(J1-1)
          call DDAFile(LuTmp(iSym),2,DVSVEC(ip_T),lTot,iAdr)
          kT = ip_T-1
          kV = ip_V-1
          do J=1,nVec_this_batch
            kOffT = kT+nDim*(J-1)
            kOffV = kV+nnBstR(iSym,2)*(J-1)+iiBstRSh(iSym,id_mySP(iSP1),2)
            do i=1,nDim
              DVSVEC(kOffV+i) = DVSVEC(kOffT+i)
            end do
          end do
          iAdr0 = iAdr0+nDim*NumCho(iSym)
        end if
        iSP1 = iSP1+nSP_this_batch
      end do
      if (iPrint >= Inf_Progress) then
        call Cho_Timer(X1,Y1)
        write(LuPri,'(6X,A,F12.2,1X,F12.2)') 'Time for read/reorder (sec)    :',(X1-X0),(Y1-Y0)
        call Cho_Flush(LuPri)
      end if
      ! Write full vectors to disk
      lTot = nnBstR(iSym,2)*nVec_this_batch
      iAdr = nnBstR(iSym,2)*(J1-1)
      call DDAFile(LuCho(iSym),1,DVSVEC(ip_V),lTot,iAdr)
      if (iPrint >= Inf_Progress) then
        call Cho_Timer(X0,Y0)
        write(LuPri,'(6X,A,F12.2,1X,F12.2)') 'Time for write (sec)           :',(X0-X1),(Y0-Y1)
        call Cho_Flush(LuPri)
      end if
    end do
  end if
end do

! Deallocation
call mma_deallocate(DVSVEC)

end subroutine Cho_XCV_DV_S
