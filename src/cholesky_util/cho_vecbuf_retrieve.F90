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
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_VecBuf_Retrieve(Vec,lVec,jVec1,iVec2,iSym,jNum,iRedC,mUsed)
!
! Thomas Bondo Pedersen, June 2006.
!
! Purpose: copy as many vectors as possible from buffer to array
!          Vec, starting at vector jVec1 and copying at most until
!          vector iVec2. On exit, jNum is the number of vectors
!          copied. On entry as well as on exit, iRedC identifies the
!          reduced set stored at location 3 (use "-1" if none or
!          unknown). On exit, mUsed is the actual amount of memory
!          used (in array Vec).
!
! NOTE: if no vectors can be copied, jNum=0 and mUsed=0 are returned
!       but execution is NOT stopped here!!!
!
! NOTE: it is assumed that the vectors are stored in their
!       respective reduced sets (thus, should only be used with
!       RUN_MODE = RUN_EXTERNAL).

use Cholesky, only: CHVBUF, InfVec, ip_CHVBUF_SYM, l_CHVBFI_SYM, l_CHVBUF_SYM, LuPri, nDimRS, nnBstR, nVec_in_Buf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lVec, jVec1, iVec2, iSym
real(kind=wp), intent(out) :: Vec(lVec)
integer(kind=iwp), intent(out) :: jNum, mUsed
integer(kind=iwp), intent(inout) :: iRedC
logical(kind=iwp) :: Full
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
integer(kind=iwp) :: iLoc, irc, iV2, iVec, jAdr, jRed, jVec, kB, kOffV, lTot, nErr, nTst
real(kind=wp) :: xNrm
logical(kind=iwp), parameter :: LocDbg = _DBG_
character(len=*), parameter :: SecNam = 'Cho_VecBuf_Retrieve'
real(kind=wp), external :: ddot_

! Initialize.
! -----------

jNum = 0
mUsed = 0

! Check that a buffer has been allocated and contains vectors in the
! requested range.
! ------------------------------------------------------------------

if (l_ChvBuf_Sym(iSym) < 1) then
  if (LocDbg) write(Lupri,*) SecNam,': returning immediately. No buffer allocated.'
  return
end if
if (l_ChvBfI_Sym(iSym) > 0) call Cho_VecBuf_Check()
if (jVec1 > nVec_in_Buf(iSym)) then
  if (LocDbg) write(Lupri,*) SecNam,': returning immediately. jVec1 = ',jVec1,'  >  nVec_in_Buf = ',nVec_in_Buf(iSym),' (sym. ', &
                             iSym,')'
  return
end if

! Count how many vectors can be copied.
! -------------------------------------

lTot = 0
Full = lTot >= lVec
jVec = jVec1-1
iV2 = min(nVec_in_Buf(iSym),iVec2)
if (.not. allocated(nDimRS)) then
  iLoc = 3
  do while ((jVec < iV2) .and. (.not. Full))
    jVec = jVec+1
    jRed = InfVec(jVec,2,iSym)
    if (jRed /= iRedC) then
      irc = 0
      call Cho_X_SetRed(irc,iLoc,jRed)
      if (irc /= 0) then
        write(Lupri,*) SecNam,': Cho_X_SetRed returned ',irc
        call Cho_Quit('Error in '//SecNam,104)
      end if
      iRedC = jRed
    end if
    lTot = lTot+nnBstR(iSym,iLoc)
    if (lTot > lVec) then
      jVec = jVec-1
      lTot = lTot-nnBstR(iSym,iLoc)
      Full = .true.
    else
      jNum = jNum+1
    end if
  end do
else
  do while ((jVec < iV2) .and. (.not. Full))
    jVec = jVec+1
    jRed = InfVec(jVec,2,iSym)
    lTot = lTot+nDimRS(iSym,jRed)
    if (lTot > lVec) then
      jVec = jVec-1
      lTot = lTot-nDimRS(iSym,jRed)
      Full = .true.
    else
      jNum = jNum+1
    end if
  end do
end if

! Copy vectors (if any).
! ----------------------

if (lTot > 0) then
  kB = ip_ChVBuf_Sym(iSym)
  if (jVec1 > 1) then
    if (.not. allocated(nDimRS)) then
      iLoc = 3
      do jVec=1,jVec1-1
        jRed = InfVec(jVec,2,iSym)
        if (iRedC /= jRed) then
          irc = 0
          call Cho_X_SetRed(irc,iLoc,jRed)
          if (irc /= 0) then
            write(Lupri,*) SecNam,': Cho_X_SetRed returned ',irc
            call Cho_Quit('Error [2] in '//SecNam,104)
          end if
          iRedC = jRed
        end if
        kB = kB+nnBstR(iSym,iLoc)
      end do
    else
      do jVec=1,jVec1-1
        jRed = InfVec(jVec,2,iSym)
        kB = kB+nDimRS(iSym,jRed)
      end do
    end if
  end if
  Vec(1:lTot) = CHVBUF(kB:kB+lTot-1)
  ! Check copy operation (may fail if molcas is compiled for
  ! 64 bit but linked to a 32 bit blas library)
  ! Note: check is not done unless it is enabled when the buffer
  ! is initialized.
  ! This only happens if _DEBUGPRINT_ is defined, so
  ! the following section is normally not executed.
  ! For debugging without turning on _DEBUGPRINT_ compilation:
  !   Call Cho_VecBuf_EnableIntegrityCheck(irc)
  ! somewhere in the code and for every subsequent call to this
  ! routine, the integrity is checked.
  if (l_ChVBfI_Sym(iSym) > 0) then
    nErr = 0
    kB = 1
    do iVec=1,jNum
      jVec = jVec1+iVec-1
      jRed = InfVec(jVec,2,iSym)
      call Cho_VecBuf_CompareNormAndSum(nDimRS(iSym,jRed),1,Vec(kB),jVec,iSym,irc)
      if (irc /= 0) then
        nErr = nErr+1
        write(LuPri,'(A,I9,A,I2,A)') 'Buffer copy failed for vector',jVec,' (sym.',iSym,')'
      end if
      kB = kB+nDimRS(iSym,jRed)
    end do
    if (nErr > 0) then
      call XFlush(LuPri)
      write(LuPri,'(A,I9,A)') 'Cho_VecBuf_Retrieve: buffer copy failed for',nErr,' vectors. Going to check buffer integrity...'
      call XFlush(LuPri)
      call Cho_VecBuf_Check()
      write(LuPri,'(A)') 'Buffer integrity checked: OK --- error occurs in the copy operation.'
#     ifdef _I8_
      write(LuPri,'(A)') 'This appears to be a 64-bit version of MOLCAS. Did you link to a 32-bit version of the BLAS library?'
#     endif
      call Cho_Quit('Cho_VecBuf_Retrieve: buffer copy failed',104)
    end if
  end if
end if

! Set memory used.
! ----------------

mUsed = lTot

! Debug: print.
! -------------

if (LocDbg) then
  write(Lupri,*)
  write(Lupri,*) SecNam,':'
  if (jNum < 1) then
    write(Lupri,*) 'No vectors copied!'
  else
    write(Lupri,*) 'Vectors ',jVec1,' to ',jVec1+jNum-1,' of symmetry ',iSym,' copied from buffer.'
    if (allocated(nDimRS)) then
      kOffV = 1
      do iVec=1,jNum
        jVec = jVec1+iVec-1
        jRed = InfVec(jVec,2,iSym)
        jAdr = InfVec(jVec,3,iSym)
        xNrm = sqrt(dDot_(nDimRS(iSym,jRed),Vec(kOffV),1,Vec(kOffV),1))
        write(Lupri,*) 'Vector:',jVec,' disk address: ',jAdr,' norm: ',xNrm
        kOffV = kOffV+nDimRS(iSym,jRed)
      end do
      nTst = kOffV-1
      if (nTst /= mUsed) call Cho_Quit('Vector dimension error in '//SecNam,104)
    end if
  end if
  call XFlush(Lupri)
end if

end subroutine Cho_VecBuf_Retrieve
