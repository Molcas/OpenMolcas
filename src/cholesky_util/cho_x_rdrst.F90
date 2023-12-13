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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_X_RdRst(ifail)
!
! T.B. Pedersen, july 2004.
!
! Purpose: read decomposition info and store in common
!          block. If ifail != 0 on exit, some error occurred and,
!          most likely, some of the restart info is not
!          defined/initialized.

use Index_Functions, only: nTri_Elem
use Cholesky, only: Cho_AdrVec, Damp, InfRed, InfRed_Hidden, InfVec, InfVec_Hidden, InfVec_N2, LuRst, MaxRed, MaxVec, nBas, nnShl, &
                    nShell, nSym, NumCho, Span, ThrCom, ThrDiag, ThrNeg, TooNeg, WarNeg, XCho_AdrVec, XnPass, XDamp, XScDiag, &
                    XSpan, XThrCom, XThrDiag, XThrNeg, XTooNeg, XWarNeg
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ifail
integer(kind=iwp), parameter :: lScr = 8
integer(kind=iwp) :: iAdr, iOpt, iSym, j, jScr(lScr), nRd, nSP_UpLim
real(kind=wp) :: dScr(lScr)
character(len=*), parameter :: SecNam = 'Cho_X_RdRst'

! Set return code.
! ----------------

ifail = 0

! Read molecular info.
! --------------------

iAdr = 0

iOpt = 2
nRd = 4
call iDAFile(LuRst,iOpt,jScr,nRd,iAdr)
nShell = jScr(2)
nnShl = jScr(3)
if (jScr(2) < 1) then
  write(u6,'(A,A,I10)') SecNam,': #shells from restart file:',jScr(2)
  ifail = 1
  call Finish_this()
  return
end if
nSP_UpLim = nTri_Elem(nShell)
if ((jScr(3) < 1) .or. (jScr(3) > nSP_UpLim)) then
  write(u6,'(A,A,I10)') SecNam,': #shell pairs from restart file:',jScr(3)
  ifail = 1
  call Finish_this()
  return
end if
if (jScr(1) /= nSym) then
  write(u6,'(A,A,I10)') SecNam,': #irreps from restart file:',jScr(1)
  ifail = 1
  call Finish_this()
  return
else
  iOpt = 2
  call iDAFile(LuRst,iOpt,jScr,nSym,iAdr)
  do iSym=1,nSym
    if (jScr(iSym) /= nBas(iSym)) then
      write(u6,'(A,A,I2,A,I10)') SecNam,': #basis functions in sym.',iSym,' from restart file:',jScr(iSym)
      ifail = 2
      call Finish_this()
      return
    end if
  end do
end if

! Read decomposition configuration info.
! --------------------------------------

iOpt = 2
nRd = 2
call iDAFile(LuRst,iOpt,jScr,nRd,iAdr)
if (jScr(1) == 0) then
  XScDiag = .false.
else if (jScr(1) == 1) then
  XScDiag = .true.
else
  write(u6,'(A,A,I10)') SECNAM,': integer flag for screening not recognized:',jScr(1)
  ifail = 2
  call Finish_this()
  return
end if
if ((jScr(2) > 0) .and. (jScr(2) < 3)) then
  XCho_AdrVec = jScr(2)
else
  write(u6,'(A,A,I10)') SECNAM,': vector file address mode not recognized:',jScr(2)
  ifail = 3
  call Finish_this()
  return
end if
if (XCho_AdrVec /= Cho_AdrVec) then
  write(u6,'(A,A,I10)') SECNAM,': vector file address mode from restart file:',XCho_AdrVec
  write(u6,'(A,A,I10)') SECNAM,': vector file address mode from runfile     :',Cho_AdrVec
  ifail = 3
  call Finish_this()
  return
end if

iOpt = 2
nRd = 8
call dDAFile(LuRst,iOpt,dScr,nRd,iAdr)
XThrCom = dScr(1)
XThrDiag = dScr(2)
XDamp(1) = dScr(3)
XDamp(2) = dScr(4)
XSpan = dScr(5)
XThrNeg = dScr(6)
XWarNeg = dScr(7)
XTooNeg = dScr(8)
ThrCom = XThrCom
ThrDiag = XThrDiag
Damp(1) = XDamp(1)
Damp(2) = XDamp(2)
Span = XSpan
ThrNeg = XThrNeg
WarNeg = XWarNeg
TooNeg = XTooNeg

! Allocate InfVec array.
! ----------------------

call mma_allocate(InfVec_Hidden,MaxVec,InfVec_N2,nSym,Label='InfVec_Hidden')
InfVec => InfVec_Hidden

! Allocate and initialize (read) InfRed array.
! --------------------------------------------

iOpt = 2
nRd = 1
call iDAFile(LuRst,iOpt,jScr,nRd,iAdr)
MaxRed = jScr(1)
XnPass = MaxRed
if (MaxRed < 1) then
  write(u6,'(A,A,I10)') SecNam,': #reduced sets from restart file:',MaxRed
  ifail = 4
  call Finish_this()
  return
else
  call mma_allocate(InfRed_Hidden,MaxRed,Label='InfRed_Hidden')
  InfRed => InfRed_Hidden
  iOpt = 2
  call iDAFile(LuRst,iOpt,InfRed,size(InfRed),iAdr)
  if (InfRed(1) /= 0) then
    write(u6,'(A,A,I10)') SecNam,': disk address of 1st reduced set:',InfRed(1)
    ifail = 5
    call Finish_this()
    return
  end if
end if

! Read InfVec array.
! ------------------

do iSym=1,nSym
  iOpt = 2
  nRd = 1
  call iDAFile(LuRst,iOpt,jScr,nRd,iAdr)
  if (jScr(1) /= NumCho(iSym)) then
    write(u6,'(A,A,I2,A,I10)') SecNam,': #Cholesky vectors (sym.',iSym,'): ',NumCho(iSym)
    write(u6,'(A,A,I10)') SecNam,': ....and from restart file: ',jScr(iSym)
    ifail = 6
    call Finish_this()
    return
  else
    if (NumCho(iSym) < 1) then
      InfVec(:,:,iSym) = 0
    else
      InfVec(:,:,iSym) = 0
      do j=1,size(InfVec,2)
        iOpt = 2
        call iDAFile(LuRst,iOpt,InfVec(:,j,iSym),NumCho(iSym),iAdr)
      end do
    end if
  end if
end do

! Return.
! -------

call Finish_this()

contains

! failures jump to this point
subroutine Finish_this()

  if (ifail /= 0) write(u6,'(A,A)') SecNam,': refusing to read more restart info!'

end subroutine Finish_this

end subroutine Cho_X_RdRst
