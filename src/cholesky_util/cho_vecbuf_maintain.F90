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

subroutine Cho_VecBuf_Maintain(irc,iRed,DoTime,DoStat)
!
! Purpose: maintain Cholesky vector buffer:
!
!          1) reorder vectors in buffer to current reduced set
!             storage (defined by location 2),
!
!          2) if possible, read in new vectors and store in current
!             reduced set storage (defined by location 2).
!
!          It is assumed that all vectors in the buffer are stored
!          according to the reduced set identified by iRed.
!
! DoTime: time as vector I/O.
! DoStat: update statistics info (#calls to system for I/O).
!
! Return code:  irc  = 0 : success
!               irc != 0 : failure
!
! Index arrays from Cholesky modified by this routine:
!
! NVEC_IN_BUF() -- #vectors stored in buffer in each symmetry

use Cholesky, only: CHVBUF, InfVec, ip_CHVBUF_SYM, iScr, l_CHVBUF_SYM, LuPri, nnBstR, nSym, nSys_call, NumCho, NumChT, &
                    nVec_in_Buf, TDECOM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: iRed
logical(kind=iwp), intent(in) :: DoTime, DoStat
integer(kind=iwp) :: iE, iMapC, iOff3, iRedC, iRS2, iS, iSym, iVec, iVec1, iVec2, jRed, jRS3, jVec, kVec, l_VRd, lCol, Left, lRow, &
                     mUsed, nDisk, nErr, nSys, nVec, nVRd
real(kind=wp) :: C1, C2, W1, W2
real(kind=wp), pointer :: V2(:,:), V3(:,:)
real(kind=wp), allocatable :: VRd(:)
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: LocDbg = _DBG_
character(len=*), parameter :: SecNam = 'Cho_VecBuf_Maintain'

! Set return code.
! ----------------

irc = 0

! Return if there is no buffer to maintain.
! -----------------------------------------

if (.not. allocated(CHVBUF)) then
  if (LocDbg) write(Lupri,*) SecNam,': returning: no buffer to maintain!'
  return
end if

! Debug print.
! ------------

if (LocDbg) then
  write(Lupri,*)
  write(Lupri,*) '>>>>> Enter ',SecNam,' <<<<<'
  write(Lupri,*) 'iRed = ',iRed
  write(Lupri,*) 'l_ChVBuf  = ',size(CHVBUF),'   ip_ChVBuf = ',1
  write(Lupri,'(A,8I16)') 'l_ChVBuf_Sym : ',(l_ChVBuf_Sym(iSym),iSym=1,nSym)
  write(Lupri,'(A,8I16)') 'ip_ChVBuf_Sym: ',(ip_ChVBuf_Sym(iSym),iSym=1,nSym)
  write(Lupri,'(A,8I16)') 'nVec_in_Buf  : ',(nVec_in_Buf(iSym),iSym=1,nSym)
end if

! If there are no vectors yet, return.
! ------------------------------------

if (NumChT < 1) then
  if (LocDbg) then
    write(Lupri,*) SecNam,': returning: no vectors!'
    write(Lupri,*) SecNam,': NumChT = ',NumChT
  end if
  return
end if

! Check that iScr array has been allocated.
! -----------------------------------------

if (.not. allocated(iScr)) then
  write(Lupri,*) SecNam,': iScr array not allocated!'
  irc = 102
  return
end if

! Start timing.
! -------------

if (DoTime) call CWTime(C1,W1)

! Set index arrays for reduced set iRed at location 3.
! ----------------------------------------------------

if (iRed < 1) then
  nErr = 0
  do iSym=1,nSym
    if (nVec_in_Buf(iSym) /= 0) then
      write(Lupri,*) SecNam,': sym. block ',iSym,':'
      write(Lupri,*) '   iRed = ',iRed,'   #vectors in buffer: ',nVec_in_Buf(iSym),' (should be 0)'
      nErr = nErr+1
    end if
  end do
  if (nErr /= 0) then
    irc = 103
    return
  end if
  iRedC = 1
else
  iRedC = iRed
end if

call Cho_X_SetRed(irc,3,iRedC)
if (irc /= 0) then
  write(Lupri,*) SecNam,': Cho_X_SetRed returned ',irc
  irc = 104
  return
end if

! Reordering.
! ===========

do iSym=1,nSym
  if ((nnBstR(iSym,2) > 0) .and. (nVec_in_Buf(iSym) > 0)) then

    ! Check reduced set dimensions.
    ! -----------------------------

    if (nnBstR(iSym,2) > nnBstR(iSym,3)) then
      write(Lupri,*) SecNam,': dimension of reduced set at 2 is larger than that at 3'
      write(Lupri,*) 'Symmetry: ',iSym,'  ID of 3: ',iRedC
      write(Lupri,*) 'Dimension of reduced set 2: ',nnBstR(iSym,2)
      write(Lupri,*) 'Dimension of reduced set 3: ',nnBstR(iSym,3)
      write(Lupri,*) 'Unable to continue...'
      irc = 104
      return
    end if

    ! Define mapping from reduced set at location 2 to that at
    ! location 3.
    ! --------------------------------------------------------

    call Cho_RS2RS(iScr,size(iScr),2,3,iRedC,iSym)

    ! Reorder vectors.
    ! ----------------

    lRow = nnBstR(iSym,2)
    lCol = nVec_in_Buf(iSym)
    iS = ip_ChVBuf_Sym(iSym)
    iE = iS-1+lRow*lCol
    V2(1:lRow,1:lCol) => CHVBUF(iS:iE)

    lRow = nnBstR(iSym,3)
    iE = iS-1+lRow*lCol
    V3(1:lRow,1:lCol) => CHVBUF(iS:iE)

    do iVec=1,nVec_in_Buf(iSym)
      do iRS2=1,nnBstR(iSym,2)
        jRS3 = iScr(iRS2)
#       ifdef _DEBUGPRINT_
        if ((iRS2 < 1) .or. (iRS2 > size(iScr))) then
          write(LuPri,*) 'iRS2=',iRS2
          write(LuPri,*) 'SIZE(iScr)=',size(iScr)
          call Cho_Quit('RS-2-RS map error in '//SecNam,104)
        end if
        if ((jRS3 < 1) .or. (jRS3 > nnBstR(iSym,3))) then
          write(LuPri,*) 'jRS3=',JRS3
          write(LuPri,*) 'nnBstR(iSym,3)=',nnBstR(iSym,3)
          call Cho_Quit('RS-2-RS map error in '//SecNam,104)
        end if
#       endif
        V2(iRS2,iVec) = V3(jRS3,iVec)
      end do
    end do

  end if
end do
nullify(V2)
nullify(V3)

! Read in more vectors.
! =====================

nSys = 0 ! #calls to reading routine (counter)

call mma_maxDBLE(l_VRd)
call mma_allocate(VRd,l_VRd,Label='VRd')
do iSym=1,nSym
  nDisk = NumCho(iSym)-nVec_in_Buf(iSym)
# ifdef _DEBUGPRINT_
  if (nDisk < 0) call Cho_Quit('nDisk < 0 in '//SecNam,103)
# endif
  iMapC = -1
  if ((nnBstR(iSym,2) > 0) .and. (nDisk > 0)) then

    ! Compute how many more vectors can be stored in buffer taking
    ! into account the number of vectors on disk.
    ! ------------------------------------------------------------

    Left = l_ChVBuf_Sym(iSym)-nnBstR(iSym,2)*nVec_in_Buf(iSym)
    if (Left >= 0) then
      nVec = min(Left/nnBstR(iSym,2),nDisk)
    else
      call Cho_Quit('Left < 0 in '//SecNam,103)
      nVec = 0
    end if
    iVec1 = nVec_in_Buf(iSym)+1
    iVec2 = iVec1+nVec-1

    ! Read and reorder vectors.
    ! -------------------------

    iVec = iVec1
    do while (iVec <= iVec2)

      ! Read vectors.
      ! -------------

      nVRd = 0
      mUsed = 0
      call Cho_VecRd(VRd,l_VRd,iVec,iVec2,iSym,nVRd,iRedC,mUsed)
      if (nVRd < 1) call Cho_Quit('Insufficient memory for read in '//SecNam,101)
      nSys = nSys+1

      ! Reorder the vectors and store in buffer in current
      ! reduced set.
      ! --------------------------------------------------

      lRow = nnBstR(iSym,2)
      lCol = nVec_in_Buf(iSym)+nVRd
      iS = ip_ChVBuf_Sym(iSym)
      iE = iS-1+lRow*lCol

      V2(1:lRow,1:lCol) => CHVBUF(iS:iE)

      iOff3 = 0
      do kVec=1,nVRd

        jVec = iVec+kVec-1
        jRed = InfVec(jVec,2,iSym)
        if (jRed /= iRedC) then
          call Cho_X_SetRed(irc,3,jRed)
          if (irc /= 0) then
            write(Lupri,*) SecNam,': Cho_X_SetRed [2] returned ',irc
            irc = 104
            return
          end if
          iRedC = jRed
        end if

        if (jRed /= iMapC) then
          call Cho_RS2RS(iScr,size(iScr),2,3,jRed,iSym)
          iMapC = jRed
        end if

        do iRS2=1,nnBstR(iSym,2)
          jRS3 = iScr(iRS2)
#         ifdef _DEBUGPRINT_
          if ((jRS3 < 1) .or. (jRS3 > nnBstR(iSym,3))) call Cho_Quit('RS-2-RS map error [2] in '//SecNam,104)
#         endif
          V2(iRS2,jVec) = VRd(iOff3+jRS3)
        end do

        iOff3 = iOff3+nnBstR(iSym,3)

      end do

      ! Update counters.
      ! ----------------

      iVec = iVec+nVRd
      nVec_in_Buf(iSym) = nVec_in_Buf(iSym)+nVRd

    end do

  end if
end do
nullify(V2)
call mma_deallocate(VRd)

! Update global timing.
! ---------------------

if (DoStat) nSys_Call = nSys_Call+nSys

! Update global timing.
! ---------------------

if (DoTime) then
  call CWTime(C2,W2)
  tDecom(1,2) = tDecom(1,2)+C2-C1
  tDecom(2,2) = tDecom(2,2)+W2-W1
end if

! Debug print.
! ------------

if (LocDbg) then
  write(Lupri,*) 'After updating: '
  write(Lupri,'(A,8I8)') 'nVec_in_Buf  : ',(nVec_in_Buf(iSym),iSym=1,nSym)
  write(Lupri,*) '>>>>> Exit  ',SecNam,' <<<<<'
end if

end subroutine Cho_VecBuf_Maintain
