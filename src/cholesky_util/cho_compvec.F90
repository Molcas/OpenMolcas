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

subroutine Cho_CompVec(Diag,xInt,VecK,QDiag,Wrk,lWrk,iSym,iPass)
!
! Purpose: compute vectors from decomposition of the qualified
!          diagonal integral block. VecK contains the nQual(iSym)
!          vectors from that decomposition. xInt contains the
!          integral columns. The vectors are returned in xInt.
!          Wrk (dimension lWrk) is work space allocated in calling
!          routine. QDiag contains the qualified diagonals.

use ChoSwp, only: IndRed

implicit real*8(a-h,o-z)
real*8 Diag(*), xInt(*), VecK(*), QDiag(*), Wrk(lWrk)
#include "cholesky.fh"
#include "choprint.fh"
character*11 SecNam
parameter(SecNam='Cho_CompVec')
logical LocDbg
parameter(LocDbg=.false.)
integer Cho_P_IndxParentDiag
external Cho_P_IndxParentDiag

! Subtract previous vectors.
! --------------------------

call Cho_Subtr(xInt,Wrk,lWrk,iSym)

! Debug: check diagonal elements in updated integrals.
! ----------------------------------------------------

if (Cho_DiaChk .or. LocDbg) then
  Tol = Tol_DiaChk
  nErr = 0
  call Cho_P_ChkInt(xINT,Diag,iSym,nErr,Tol,.true.)
  if (nErr /= 0) then
    write(Lupri,*) SecNam,': ',nErr,' diagonal errors found!'
    write(Lupri,*) '          Integral pass: ',iPass
    write(Lupri,*) '          #Tests: ',nQual(iSym)
    !write(Lupri,*) '          Printing integrals:'
    !call Cho_Output(xINT,1,nnBstR(iSym,2),1,nQual(iSym),nnBstR(iSym,2),nQual(iSym),1,Lupri)
    call Cho_Quit('Diagonal errors in '//SecNam,104)
  else
    write(Lupri,*) SecNam,': comparison of qual. integrals and current diagonal: no errors !'
  end if
end if

! Set max. diagonal for screening.
! --------------------------------

QDmax = QDiag(1)

! Compute vectors.
! ----------------

do i=1,nQual(iSym)

  kOff0 = nnBstR(iSym,2)*(i-1)
  kK0 = nQual(iSym)*(i-1)

  ! Compute vector.
  ! ---------------

  xC = QDiag(i)
  Fac = 1.0d0/sqrt(abs(xC))
  kOff = kOff0+1
  call dScal_(nnBstR(iSym,2),Fac,xInt(kOff),1)

  ! Zero elements corresponding to zero diagonals.
  ! ----------------------------------------------

  do jAB=1,nnBstR(iSym,2)
    jAB1 = IndRed(iiBstR(iSym,2)+jAB,2)
    if (Diag(jAB1) == 0.0d0) xInt(kOff0+jAB) = 0.0d0
  end do

  ! Update diagonal.
  ! ----------------

  do jAB=1,nnBstR(iSym,2)
    jAB1 = IndRed(iiBstR(iSym,2)+jAB,2)
    kOff = kOff0+jAB
    Diag(jAB1) = Diag(jAB1)-xInt(kOff)**2
  end do

  ! Update Qdiag.
  ! -------------

  do j=i,nQual(iSym)
    QDiag(j) = QDiag(j)-VecK(kK0+j)**2
  end do
  OlDiag = QDiag(i)
  QDiag(i) = 0.0d0

  ! Get index (1st reduced set) of the parent diagonal for this
  ! vector (global index!).
  ! -----------------------------------------------------------

  iABG = Cho_P_IndxParentDiag(i,iSym)

  ! Zero treated diagonal element.
  ! ------------------------------

  call Cho_P_ZeroDiag(Diag,iSym,iABG)

  ! Check for and zero negative diagonals, count converged, and
  ! find max. diagonal element.
  ! -----------------------------------------------------------

  call Cho_ChkDia_A4(Diag,QDmax,iSym,nNeg,nNegT,nConv,xM,yM,zM)

  ! nnZTot is the total number of zeroed negative diagonals,
  ! updated for statistics purposes.
  ! --------------------------------------------------------

  nnZTot = nnZTot+nNeg

  ! Subtract this vector from remaining integral columns,
  ! M([gd],{ab}_j) -= K({ab}_j,i)*L([gd],NumCho+i)
  ! -----------------------------------------------------

  kOff = kOff0+1
  do j=i+1,nQual(iSym)
    Fac = -VecK(kK0+j)
    kInt = nnBstR(iSym,2)*(j-1)+1
    call dAXPY_(nnBstR(iSym,2),Fac,xInt(kOff),1,xInt(kInt),1)
  end do

  ! Print.
  ! ------

  if (iPrint >= Inf_Progress) then
    iVec = NumCho(iSym)+i
    iVecT = NumChT+i
    do j=i+1,nQual(iSym)
      xM = max(xM,QDiag(j))
    end do
    write(Lupri,'(I3,3(1X,I9),2(1X,D11.3),2(1X,I4),1X,D11.3)') iSym,iVec,iVecT,iABG,xC,OlDiag,nConv,nNeg,xM
  end if

end do

if (iPrint >= Inf_Progress) call Cho_Flush(Lupri)

end subroutine Cho_CompVec
