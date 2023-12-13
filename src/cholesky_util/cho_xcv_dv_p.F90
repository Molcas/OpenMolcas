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

subroutine Cho_XCV_DV_P(irc,SP_BatchDim,nSP_Batch,id_mySP,n_mySP,NVT,l_NVT)

#ifdef _MOLCAS_MPP_
#ifndef _GA_
use Para_Info, only: nProcs
#endif
use Cholesky, only: iiBstRSh, INF_PROGRESS, IPRINT, LuCho, LuPri, LuTmp, nnBstR, nnBstRSh, nSym, NumCho
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp
#endif
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nSP_Batch, SP_BatchDim(nSP_Batch), n_mySP, id_mySP(n_mySP), l_NVT, NVT(l_NVT)
#ifdef _MOLCAS_MPP_
#include "mafdecls.fh"
integer(kind=iwp) :: g_a, iAdr, iAdr0, iBatch, iSP, iSP1, iSP2, iSP_, iSP_Batch, iSym, J0, J1, J2, kV, l_IDV, l_Mem, l_numV, l_V, &
                     lTot, max_vector_dim, my_nV, myEnd, myStart, nBatch, nDim, nSP_this_batch, nVec_per_batch, nVec_this_batch
real(kind=wp) :: X0, X1, Y0, Y1
logical(kind=iwp) :: ok
#ifndef _GA_
integer(kind=iwp) :: Jen, Jst
#else
integer(kind=iwp) :: i, j
#endif
integer(kind=iwp), allocatable :: GADIST(:), GAMNCH(:), GANUMV(:)
real(kind=wp), allocatable :: GAVEC(:)
integer(kind=iwp), parameter :: iOpt = 1
character(len=*), parameter :: SecNam = 'Cho_XCV_DV_P'
logical(kind=iwp), external :: ga_create, ga_destroy

! Init return code
irc = 0

! Find max vector dimension
max_vector_dim = nnBstR(1,2)
do iSym=1,nSym
  max_vector_dim = max(max_vector_dim,nnBstR(iSym,2))
end do
if (max_vector_dim < 1) then
  irc = -2
  return
end if

! Allocate my vector counter array
call mma_allocate(GAMNCH,nSym,Label='GAMNCH')
GAMNCH(:) = 0

! Allocate index arrays for GA use
l_numV = 1
call mma_allocate(GANUMV,l_numV,Label='GANUMV')

! Figure out how much 1/3 of total memory is
call mma_maxDBLE(l_Mem)
l_Mem = l_Mem/3
if (l_Mem < max_vector_dim) then
  ! OK, let us try 2/3 then...
  l_Mem = l_Mem*2
  if (l_Mem < max_vector_dim) then
    irc = -1
    return
  end if
end if

! distribute vectors, one symmetry at a time
do iSym=1,nSym
  if (iPrint >= Inf_Progress) then
    write(LuPri,'(/,A,I2,/,A)') 'Distributing vectors, symmetry',iSym,'--------------------------------'
    write(LuPri,'(3X,A,I8)') 'Total number of vectors:',NVT(iSym)
    write(LuPri,'(3X,A,I8)') 'Local number of vectors:',NumCho(iSym)
    write(LuPri,'(3X,A,I8)') 'Vector dimension       :',nnBstR(iSym,2)
    write(LuPri,'(3X,A,I8)') 'Shell pair batches     :',nSP_Batch
    call XFlush(LuPri)
  end if
  if ((NVT(iSym) > 0) .and. (nnBstR(iSym,2) > 0)) then
    ! Set up batching, ensuring that the number of vectors is
    ! the same on all nodes
    GANUMV(1) = min(l_Mem/nnBstR(iSym,2),NVT(iSym))
    call Cho_GAIGOp(GANUMV,1,'min')
    nVec_per_batch = GANUMV(1)
    if (nVec_per_batch < 1) call Cho_Quit('Insufficient memory for batching in '//SecNam,101)
    nBatch = (NVT(iSym)-1)/nVec_per_batch+1
    ! Allocate memory for vectors (will hold local as well as
    ! full vectors)
    l_V = nnBstR(iSym,2)*nVec_per_batch
    call mma_allocate(GAVEC,l_V,Label='GAVEC')
    ! Allocate vector ID array for distribution
    l_IDV = nVec_per_batch
    call mma_allocate(GADIST,l_IDV,Label='GADIST')
    ! Create global array with evenly distributed chunks
    ok = ga_create(mt_dbl,nnBstR(iSym,2),nVec_per_batch,'GA_XCV',0,0,g_a)
    if (.not. ok) call Cho_Quit(SecNam//': ga_create() failed!',101)
    if (iPrint >= Inf_Progress) then
      write(LuPri,'(3X,A,I8)') 'Vector batches         :',nBatch
      call XFlush(LuPri)
    end if
    ! Distribute in batches
    do iBatch=1,nBatch
      ! Sync: all must be at the same batch
      ! (so that we do not put while another is getting)
      call GASync()
      ! Determine number of vectors in this batch
      if (iBatch == nBatch) then
        nVec_this_batch = NVT(iSym)-nVec_per_batch*(nBatch-1)
      else
        nVec_this_batch = nVec_per_batch
      end if
      ! First and last vector in this batch
      J1 = nVec_per_batch*(iBatch-1)+1
      J2 = J1+nVec_this_batch-1
      if (iPrint >= Inf_Progress) then
        write(LuPri,'(3X,A,I8,/,3X,A)') 'Vector batch number:',iBatch,'++++++++++++++++++++++++++++'
        write(LuPri,'(6X,A,I8)') 'Number of vectors in this batch:',nVec_this_batch
        write(LuPri,'(6X,A,I8,1X,I8)') 'First and last vector          :',J1,J2
        call XFlush(LuPri)
        call CWTime(X0,Y0)
      end if
      ! Loop through SP blocks of vectors in this vector batch
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
          ! Read vector block
          lTot = nDim*nVec_this_batch
          iAdr = iAdr0+nDim*(J1-1)
          call DDAFile(LuTmp(iSym),2,GAVEC,lTot,iAdr)
          iAdr0 = iAdr0+nDim*NVT(iSym)
          ! Put vector block into global array
          kV = 1
          do iSP_=iSP1,iSP2
            iSP = id_mySP(iSP_)
            if (nnBstRSh(iSym,iSP,2) > 0) then
              myStart = iiBstRSh(iSym,iSP,2)+1
              myEnd = myStart+nnBstRSh(iSym,iSP,2)-1
              call ga_put(g_a,myStart,myEnd,1,nVec_this_batch,GAVEC(kV),nDim)
              kV = kV+nnBstRSh(iSym,iSP,2)
            end if
          end do
        end if
        iSP1 = iSP1+nSP_this_batch
      end do
      if (iPrint >= Inf_Progress) then
        call CWTime(X1,Y1)
        write(LuPri,'(6X,A,F12.2,1X,F12.2)') 'Time for read/ga_put (sec)     :',(X1-X0),(Y1-Y0)
        call XFlush(LuPri)
      end if
      ! Sync: all must be done putting data into g_a
      call GASync()
      ! Compute vector distribution for this batch
      my_nV = 0
      call Cho_P_Distrib_Vec(J1,J2,GADIST,my_nV)
      if (iPrint >= Inf_Progress) call CWTime(X0,Y0)
      ! Get vectors from global array
      J0 = J1-1
      kV = 1
#     ifdef _GA_
      do i=1,my_nV
        J = GADIST(i)-J0
        call ga_get(g_a,1,nnBstR(iSym,2),J,J,GAVEC(kV),nnBstR(iSym,2))
        kV = kV+nnBstR(iSym,2)
      end do
#     else
      if (my_nV > 0) then
        Jst = GADIST(1)-J0
        Jen = GADIST(my_nV)-J0
        call ga_get_striped(g_a,1,nnBstR(iSym,2),Jst,Jen,GAVEC(kV),nnBstR(iSym,2),nProcs)
        kV = kV+nnBstR(iSym,2)*my_nV
      end if
      ! VVP: First performs RMA and only then I/O
      call GASync()
#     endif
      ! Write vectors to disk
      lTot = nnBstR(iSym,2)*my_nV
      if (lTot > 0) then
        iAdr = nnBstR(iSym,2)*GAMNCH(iSym)
        call DDAFile(LuCho(iSym),iOpt,GAVEC,lTot,iAdr)
      end if
      if (iPrint >= Inf_Progress) then
        call CWTime(X1,Y1)
        write(LuPri,'(6X,A,F12.2,1X,F12.2)') 'Time for ga_get/write (sec)    :',(X1-X0),(Y1-Y0)
        call XFlush(LuPri)
      end if
      ! Update my vector counter
      GAMNCH(iSym) = GAMNCH(iSym)+my_nV
    end do
    ! Destroy global array
    ok = ga_destroy(g_a)
    ! Deallocations
    call mma_deallocate(GADIST)
    call mma_deallocate(GAVEC)
  end if
end do

! Deallocations
call mma_deallocate(GANUMV)
call mma_deallocate(GAMNCH)
#else

#include "macros.fh"
unused_var(SP_BatchDim)
unused_var(id_mySP)
unused_var(NVT)

irc = 999 ! should never be called in serial installation

return

#endif

end subroutine Cho_XCV_DV_P
