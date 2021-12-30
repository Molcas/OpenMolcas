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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ComputeFuncER(ERFun,CMO,nBas,nOcc,nFro,nSym,Timing)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: compute Edmiston-Ruedenberg functional.
!
! =====================================================
!    WORKS *ONLY* WITH CHOLESKY DECOMPOSED INTEGRALS
! =====================================================

implicit none
real*8 ERFun
real*8 CMO(*)
integer nSym
integer nBas(nSym), nOcc(nSym), nFro(nSym)
logical Timing
#include "WrkSpc.fh"
character*13 SecNam
parameter(SecNam='ComputeFuncER')
character*80 Txt
real*8 FracMem
integer irc, ipERFun, lERFun, iSym, i, kOff, nFroT
integer nOccT(8)

! Initializations.
! ----------------

irc = 0

FracMem = 0.0d0 ! no buffer allocated
call Cho_X_Init(irc,FracMem)
if (irc /= 0) then
  write(Txt,'(A,I4)') 'Cho_X_Init returned',irc
  call SysAbendMsg(SecNam,'Cholesky initialization failure!',Txt)
end if

! Check dimensions.
! -----------------

call ERChk_Localisation(irc,nBas,nOcc,nFro,nSym)
if (irc /= 0) then
  write(Txt,'(A,I4)') 'ERChk_Localisation returned',irc
  call SysAbendMsg(SecNam,'Cholesky initialization mismatch!',Txt)
end if

! Compute ER functional.
! ----------------------

nOccT(1) = nOcc(1)+nFro(1)
do iSym=2,nSym
  nOccT(iSym) = nOcc(iSym)+nFro(iSym)
end do

lERFun = nOccT(1)
nFroT = nFro(1)
do iSym=2,nSym
  lERFun = lERFun+nOccT(iSym)
  nFroT = nFroT+nFro(iSym)
end do

call GetMem('ERFun','Allo','Real',ipERFun,lERFun)
ERFun = 0.0d0
call EvalERFun(ERFun,Work(ipERFun),CMO,nOccT,nSym,Timing)
if (nFroT > 0) then
  kOff = ipERFun-1
  do iSym=1,nSym
    do i=1,nFro(iSym)
      ERFun = ERFun-Work(kOff+i)
    end do
    kOff = kOff+nOccT(iSym)
  end do
end if
call GetMem('ERFun','Free','Real',ipERFun,lERFun)

! Finalizations.
! --------------

call Cho_X_Final(irc)
if (irc /= 0) then
  write(Txt,'(A,I4)') 'Cho_X_Final returned',irc
  call SysAbendMsg(SecNam,'Cholesky finalization failure!',Txt)
end if

end subroutine ComputeFuncER
