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

subroutine Cho_SubScr_Dia(ChoVec,nVec,iSym,iLoc,DSPNorm)
!
! Purpose: compute diagonal from Cholesky vectors and find norm
!          for each shell pair. Needed for screening of vector
!          subtraction. Which norm is used is determined by
!          string DSPNorm:
!
!          DSPNorm = 'Max' : max. element
!          DSPNorm = 'Fro' : Frobenius norm
!
!          Any other norm is taken to be 'Max'.

use ChoSubScr, only: DSPNm, DSubScr
use ChoSwp, only: iiBstRSh, nnBstRSh
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: ChoVec(*)
integer(kind=iwp) :: nVec, iSym, iLoc
character(len=*) :: DSPNorm
#include "cholesky.fh"
integer(kind=iwp) :: iAB, iAB1, iAB2, iSP, iVec, kOff, lstr
character(len=3) :: myDSPNorm
character(len=*), parameter :: SecNam = 'Cho_SubScr_Dia'

#ifdef _DEBUGPRINT_
if ((iLoc < 1) .or. (iLoc > 3)) call Cho_Quit('iLoc error in '//SecNam,104)
if ((iSym < 1) .or. (iSym > nSym)) call Cho_Quit('iSym error in '//SecNam,104)
if (.not. Cho_SScreen) call Cho_Quit('Cho_SScreen is .False. in '//SecNam,104)
#endif

! Initialize and check for early return.
! --------------------------------------

call FZero(DSubScr,nnBstR(iSym,iLoc))
call FZero(DSPNm,nnShl)
if ((nVec < 1) .or. (nnBstR(iSym,iLoc) < 1)) return

! Compute diagonal.
! -----------------

do iVec=1,nVec
  kOff = nnBstR(iSym,iLoc)*(iVec-1)
  do iAB=1,nnBstR(iSym,iLoc)
    DSubScr(iAB) = DSubScr(iAB)+ChoVec(kOff+iAB)*ChoVec(kOff+iAB)
  end do
end do

! Find diagonal norm in each shell pair.
! --------------------------------------

lstr = len(DSPNorm)
if (lstr < 3) then
  myDSPNorm = 'MAX'
else
  myDSPNorm = DSPNorm(1:3)
  call UpCase(myDSPNorm)
end if

#ifdef _DEBUGPRINT_
if (lstr < 1) then
  write(Lupri,*) SecNam,': input norm: (null string)'
else if (lstr < 3) then
  write(Lupri,*) SecNam,': input norm: ',DSPNorm,' (incomplete)'
else
  write(Lupri,*) SecNam,': input norm: ',DSPNorm
end if
write(Lupri,*) SecNam,': norm used : ',myDSPNorm
#endif

if (myDSPNorm == 'MAX') then
  do iSP=1,nnShl
    iAB1 = iiBstRSh(iSym,iSP,iLoc)+1
    iAB2 = iAB1+nnBstRSh(iSym,iSP,iLoc)-1
    do iAB=iAB1,iAB2
      DSPNm(iSP) = max(DSPNm(iSP),DSubScr(iAB))
    end do
  end do
else if (myDSPNorm == 'FRO') then
  do iSP=1,nnShl
    iAB1 = iiBstRSh(iSym,iSP,iLoc)+1
    iAB2 = iAB1+nnBstRSh(iSym,iSP,iLoc)-1
    do iAB=iAB1,iAB2
      DSPNm(iSP) = DSPNm(iSP)+DSubScr(iAB)*DSubScr(iAB)
    end do
    DSPNm(iSP) = sqrt(DSPNm(iSP))
  end do
else
  write(Lupri,*) SecNam,': WARNING: unkown norm: ',DSPNorm
  write(Lupri,*) SecNam,': WARNING: using max element...'
  do iSP=1,nnShl
    iAB1 = iiBstRSh(iSym,iSP,iLoc)+1
    iAB2 = iAB1+nnBstRSh(iSym,iSP,iLoc)-1
    do iAB=iAB1,iAB2
      DSPNm(iSP) = max(DSPNm(iSP),DSubScr(iAB))
    end do
  end do
end if

end subroutine Cho_SubScr_Dia
