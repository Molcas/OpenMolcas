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

subroutine TRAMO(LBUF,OUTBUF,nOUTBUF,X1,nX1,X2,nX2,X3,nX3,VXPQ,nVXPQ,CMO,iDsk,mOVX)
! Purpose: two-electron transformation routine.
!          Transformation is made in core if all half-transformed
!          integrals fit. Otherwise sorted integrals
!          are written onto unit LUHALF.

#ifdef _HDF5_QCM_
use hdf5_utils, only: datadim, datadim_bound, datatag, file_id, hdf5_put_data, tagx
use motra_global, only: ihdf5
use stdalloc, only: mma_allocate, mma_deallocate
#endif
use motra_global, only: FnHalf, IAD13, iPrint, ISP, ISQ, ISR, ISS, LMOP, LMOQ, LMOR, LMOS, LTUVX, LuHalf, LuTwoMO, NBP, NBPQ, NBQ, &
                        NBR, NBRS, NBS, NOP, NOQ, NOR, NOS, NOVX
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: LBUF, nOUTBUF, nX1, nX2, nX3, nVXPQ, mOVX
real(kind=wp), intent(inout) :: OUTBUF(nOUTBUF), X1(nX1), VXPQ(nVXPQ)
real(kind=wp), intent(out) :: X2(nX2), X3(nX3)
real(kind=wp), intent(in) :: CMO(*)
integer(kind=iwp), intent(out) :: iDsk(3,mOVX)
#include "tratoc.fh"
integer(kind=iwp) :: I, IAD14, IAD14_, inBuf, IOPT, IOUT, IPQ, IPQMAX, IPQST, IPQUT, IPQX, IRC, IRSST, IST, ISTMOT, IVX, IX1, IX2, &
                     KBUF1, LOQ, LPKREC, LPQ, LVXPQ, NBYTES, NP, NPQ, NQ, NQM, NT, NTUVX, NUMAX, NV, NX, NXM
logical(kind=iwp) :: INCORE
#ifdef _HDF5_QCM_
integer(kind=iwp) :: iout_total
integer(kind=iwp), save :: total_number_2ints = 0
real(kind=wp), allocatable :: tmpbuf(:)
#endif

! Set some constants

iDsk(1,:) = -1
iDsk(2,:) = 0
iDsk(3,:) = 1

#ifdef _HDF5_QCM_
if (ihdf5 == 1) then
  !call mma_allocate(tmpbuf,nx1,label='tmpbuf')
  ! Stefan's dirty hack to increase the # of BF working in MOTRA with HDF5
  call mma_allocate(tmpbuf,(nop+1)*(noq+1)*(nor+1)*(nos+1),label='tmpbuf')
  tmpbuf(:) = 0
  write(u6,*) 'size of tmpbuf:',(nop+1)*(noq+1)*(nor+1)*(nos+1)
  iout_total = 0
end if
#endif
KBUF1 = 1
LTUVX = 0
IOUT = 0
IAD14 = 0
IPQUT = 0
IPQMAX = NBPQ
INCORE = .false.
if (NBPQ*NOVX > nVXPQ) then
  INCORE = .true.
  IPQMAX = nVXPQ/NOVX
  IPQMAX = (nVXPQ-IPQMAX)/NOVX
  IPQMAX = (nVXPQ-IPQMAX)/NOVX
  IPQMAX = (nVXPQ-IPQMAX)/NOVX
  write(u6,'(6X,A)') 'OUT OF CORE TRANSFORMATION'
  call DANAME_MF(LUHALF,FNHALF)
end if

! Start loop over sorted AO-integrals: npq pq-pairs in each buffer

if (IPRINT >= 20) write(u6,*) 'TRAMO 001'
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
    IPQUT = IPQUT+1

    ! Read in a block of integrals NPQ pq values

    if (LPQ == NPQ) then
      call RDORD(IRC,IOPT,ISP,ISQ,ISR,ISS,X1,LBUF,NPQ)
      call GADSum(X1,LBUF)
      IOPT = 2
      LPQ = 0
      IRSST = 1-NBRS
    end if
    LPQ = LPQ+1
    IRSST = IRSST+NBRS

    ! Start transformation of this pq pair

    if (ISR == ISS) then
      call SQUARE(X1(IRSST),X2,1,NBS,NBS)
      if (IPRINT >= 30) then
        write(u6,'(6X,A,2I4)') 'AO Integrals for the pair',NP,NQ
        call SQPRT(X2,NBR)
      end if
      if (NBR*NBS*NOS > 0) call DGEMM_('T','N',NBR,NOS,NBS,One,X2,NBS,CMO(LMOS),NBS,Zero,X3,NBR)
      if (NBR*NOR > 0) call DGEMM_Tri('T','N',NOR,NOR,NBR,One,X3,NBR,CMO(LMOR),NBR,ZERO,X2,NOR)
    else
      if (NBR*NBS*NOS > 0) call DGEMM_('T','N',NBR,NOS,NBS,One,X1(IRSST),NBS,CMO(LMOS),NBS,Zero,X3,NBR)
      if (NOS*NBR*NOR > 0) call DGEMM_('T','N',NOS,NOR,NBR,One,X3,NBR,CMO(LMOR),NBR,Zero,X2,NOS)
    end if

    ! Sort the matrix X2 into VXPQ (sort after PQ instead of VX)
    ! if INCORE also write sorted integrals on unit LUHALF.

    if (IPQUT > IPQMAX) then
      IPQUT = 1
      IST = 1
      do I=1,NOVX
        call PKR8(0,IPQMAX,NBYTES,VXPQ(IST),VXPQ(IST))
        LPKREC = (NBYTES+RtoB-1)/RtoB

        IAD14_ = IAD14
        call dDAFILE(LUHALF,1,VXPQ(IST),LPKREC,IAD14)
        call iDAFILE(LUHALF,1,iDsk(1,I),2,IAD14)
        iDsk(1,I) = IAD14_
        iDsk(2,I) = LPKREC
        iDsk(3,I) = iDsk(3,I)+IPQMAX

        IST = IST+IPQMAX
      end do
    end if

    IPQX = IPQUT
    do I=1,NOVX
      VXPQ(IPQX) = X2(I)
      IPQX = IPQX+IPQMAX
    end do
  end do !  NQ
end do   !  NP

LVXPQ = NBPQ*NOVX
if (IPRINT >= 20) write(u6,*) 'TRAMO 002'
if (IPRINT >= 30) then
  write(u6,'(A,I6)') ' HALF TRANSFORMED INTEGRALS:',LVXPQ
  write(u6,'(1X,10F11.6)') VXPQ(1:LVXPQ)
end if

! Empty last buffers

if (INCORE) then
  IST = 1
  do I=1,NOVX
    call PKR8(0,IPQMAX,NBYTES,VXPQ(IST),VXPQ(IST))
    LPKREC = (NBYTES+RtoB-1)/RtoB

    IAD14_ = IAD14
    call dDAFILE(LUHALF,1,VXPQ(IST),LPKREC,IAD14)
    call iDAFILE(LUHALF,1,iDsk(1,I),2,IAD14)
    iDsk(1,I) = IAD14_
    iDsk(2,I) = LPKREC
    iDsk(3,I) = iDsk(3,I)+IPQMAX

    IST = IST+IPQMAX
  end do
  IAD14 = 0
end if

! First half transformation is now done. VXPQ contains half trans-
! formed integrals for this symmetry block, if not INCORE
! otherwise integrals or one VX pair at a time will be read
! from unit LUHALF.

if (IPRINT >= 20) write(u6,*) 'TRAMO 004'
IVX = 0
do NV=1,NOR
  NXM = NV
  if (ISS /= ISR) NXM = NOS
  do NX=1,NXM
    IPQST = 1+NBPQ*IVX
    IVX = IVX+1

    ! Read one buffer of integrals back into core if INCORE

    if (INCORE) then

      ! Back chain buffers

      inBuf = iDsk(3,iVX)
      IPQ = inBuf-IPQMAX
      do
        IAD14 = iDsk(1,iVX)  ! Disk Address of the next buffer
        if (IAD14 < 0) exit
        LPKREC = iDsk(2,iVX)  ! Length of the  buffer
        call dDAFILE(LUHALF,2,VXPQ(inBuf),LPKREC,IAD14)
        call iDAFILE(LUHALF,2,iDsk(1,iVX),2,IAD14)
        call UPKR8(0,IPQMAX,NBYTES,VXPQ(inBuf),VXPQ(IPQ))
        IPQ = IPQ-IPQMAX
      end do

      IPQST = 1
    end if
    if (ISP /= ISR) then
      if (ISP == ISQ) then
        call SQUARE(VXPQ(IPQST),X2,1,NBQ,NBQ)
        if (NBP*NBQ*NOQ > 0) call DGEMM_('T','N',NBP,NOQ,NBQ,One,X2,NBQ,CMO(LMOQ),NBQ,Zero,X1,NBP)
        If (NOP*NBP > 0) call DGEMM_Tri('T','N',NOP,NOP,NBP,One,X1,NBP,CMO(LMOP),NBP,Zero,X2,NOP)
        IX2 = (NOP+NOP**2)/2
      else
        if (NBP*NBQ*NOQ > 0) call DGEMM_('T','N',NBP,NOQ,NBQ,One,VXPQ(IPQST),NBQ,CMO(LMOQ),NBQ,Zero,X1,NBP)
        if (NOQ*NBP*NOP > 0) call DGEMM_('T','N',NOQ,NOP,NBP,One,X1,NBP,CMO(LMOP),NBP,Zero,X2,NOQ)
        IX2 = NOP*NOQ
      end if
    else
      if (ISP == ISQ) then
        call SQUARE(VXPQ(IPQST),X2,1,NBP,NBP)
        if (NBP*NBQ*NOQ > 0) call DGEMM_('T','N',NBP,NOQ,NBQ,One,X2,NBQ,CMO(LMOQ),NBQ,Zero,X1,NBP)
      else
        if (NBP*NBQ*NOQ > 0) call DGEMM_('T','N',NBP,NOQ,NBQ,One,VXPQ(IPQST),NBQ,CMO(LMOQ),NBQ,Zero,X1,NBP)
      end if

      ! X1 now contains the matrix (pu/vx) for all p and all u
      ! for a fixed pair vx.
      ! If the first index t=v then the range of u is x,umax
      ! Otherwise the range is 1,umax. a loop over t is necessary her

      ISTMOT = LMOP+NBP*(NV-1)
      IX2 = 1
      do NT=NV,NOP
        NUMAX = NOQ
        if (ISP == ISQ) NUMAX = NT
        LOQ = NUMAX
        if (NT == NV) LOQ = NUMAX-NX+1
        IX1 = 1
        if (NT == NV) IX1 = 1+NBP*(NX-1)
        if (LOQ > 0) then
          if (NBP == 0) then
            X2(IX2:IX2+LOQ-1) = Zero
          else
            call DGEMM_('T','N',LOQ,1,NBP,One,X1(IX1),NBP,CMO(ISTMOT),NBP,Zero,X2(IX2),LOQ)
          end if
        end if
        IX2 = IX2+LOQ
        ISTMOT = ISTMOT+NBP
      end do
      IX2 = IX2-1
    end if

    ! Move integrals to output buffer and write them on LUTWOMO

    do NTUVX=1,IX2
      if (IOUT == nTraBuf) then
        IOUT = 0
        if (IPRINT >= 5) then
          write(u6,'(A)') 'SAVE TRANSFORMED INTEGRALS:'
          write(u6,'(A,3I8)') 'LU,nTraBuf,IDISK =',LUTWOMO,nTraBuf,IAD13
        end if
        call dDAFILE(LUTWOMO,1,OUTBUF(KBUF1),nTraBuf,IAD13)
#       ifdef _HDF5_QCM_
        !write(u6,*) 'SAVE TRANSFORMED INTEGRALS:',nTraBuf,KBUF1,iout_total,iout_total+nTraBuf
        if (ihdf5 == 1) then
          call dcopy_(nTraBuf,OUTBUF(KBUF1),1,tmpbuf(iout_total+1),1)
          iout_total = iout_total+nTraBuf
        end if
#       endif
        KBUF1 = nTraBuf+2-KBUF1
      end if
      IOUT = IOUT+1
      LTUVX = LTUVX+1
      OUTBUF(IOUT+KBUF1-1) = X2(NTUVX)
    end do
  end do
end do

! Empty last buffer (which is never empty)

if (IPRINT >= 5) then
  write(u6,'(A)') 'SAVE TRANSFORMED INTEGRALS:'
  write(u6,'(A,3I8)') 'LU,nTraBuf,IDISK =',LUTWOMO,nTraBuf,IAD13
end if
call dDAFILE(LUTWOMO,1,OUTBUF(KBUF1),nTraBuf,IAD13)
if (IPRINT >= 10) then
  write(u6,'(1X,A,4I2,A,I4)') 'TRANSFORMED INTEGRALS FOR SYMMETRY BLOCK',ISP,ISQ,ISR,ISS,' IOUT=',IOUT
  write(u6,'(1X,10F12.6)') OUTBUF(KBUF1:KBUF1+IOUT-1)
end if

#ifdef _HDF5_QCM_
if (ihdf5 == 1) then
  !> and the last buffer here as well...
  call dcopy_(iout,OUTBUF(KBUF1),1,tmpbuf(iout_total+1),1)
  iout_total = iout_total+iout
  !> put data to file
  write(datatag,'(a1,i3,a1,i3,a1,i3,a1,i3)') 'p',isp,'q',isq,'r',isr,'s',iss
  tagx = 16
  datadim(1) = 1
  datadim_bound = 1
  call hdf5_put_data(file_id(1),'XXXXXX',datadim,iout_total)
  datadim(1) = iout_total
  datadim_bound = 1
  call hdf5_put_data(file_id(1),'XXXXXX',datadim,tmpbuf)

  if (IPRINT >= 5) then
    total_number_2ints = total_number_2ints+iout_total
    write(u6,'(1X,a,i10)') ' total number of integrals (sym)  ...',iout_total
    write(u6,'(1X,a,i10)') ' total number of integrals so far ...',total_number_2ints
  end if

  call mma_deallocate(tmpbuf)
end if
#endif

if (INCORE) call DACLOS(LUHALF)

return

end subroutine TRAMO
