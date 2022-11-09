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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               1993, Per-Olof Widmark                                 *
!***********************************************************************

subroutine GetOrd(rc,Square,nSym,nBas,nSkip)
!***********************************************************************
!                                                                      *
!    Purpose: Read the table of content from the OrdInt file           *
!                                                                      *
!    Calling parameters:                                               *
!    nSym   : contains on return the number of irred. rep.             *
!    nBas   : contains on return the number of basis functions per     *
!             irred. rep.                                              *
!    nSkip  : contains on return the skipping flag per irred. rep.     *
!    rc     : return code                                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher and P.O. Widmark                                 *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use TwoDat, only: AuxTwo, isBas, isDAdr, isMxDa, isOrd, isPkPa, isPkTh, isSkip, isSym, nBatch, rcTwo, TocTwo
use Pack_mod, only: isPack, PkThrs
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: rc, nSym, nBas(0:7), nSkip(0:7)
logical(kind=iwp), intent(out) :: Square
#include "Molcas.fh"
integer(kind=iwp) :: iBatch, ijPair, iPack, iSyBlk, iSym, iTab, jSym, klPair, kSym, lSym, mxDAdr, nPairs, ntBas
logical(kind=iwp) :: DoCholesky
character(len=*), parameter :: TheName = 'GetOrd'

!----------------------------------------------------------------------*
! Start the procedure                                                  *
!----------------------------------------------------------------------*
rc = rcTwo%good

!----------------------------------------------------------------------*
! Cholesky cheating: set output variables regardless of OrdInt file    *
! Note that all elements of nSkip are zero (no skipping)               *
! Square is set false (arbitrarily)                                    *
! Skip data relating to this file (Toc etc.) entirely....              *
!----------------------------------------------------------------------*
call decideoncholesky(DoCholesky)
if (DoCholesky) then
  call Get_iScalar('nSym',nSym)
  call Get_iArray('nBas',nBas,nSym)
  nSkip(:) = 0
  Square = .false.
  return
end if
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
if (.not. AuxTwo%Opn) then
  rc = rcTwo%TC01
  call SysAbendMsg(TheName,'The ORDINT file has not been opened',' ')
end if
!----------------------------------------------------------------------*
! Check the ordering parameter                                         *
!----------------------------------------------------------------------*
if ((TocTwo(isOrd) < 0) .or. (TocTwo(isOrd) > 1)) then
  rc = rcTwo%TC02
  call SysWarnMsg(TheName,'The file carries an invalid ordering parameter',' ')
  call SysValueMsg('TocTwo(isOrd)',TocTwo(isOrd))
end if
Square = TocTwo(isOrd) == 1
!----------------------------------------------------------------------*
! Check the number of symmetry operations                              *
!----------------------------------------------------------------------*
nSym = TocTwo(isSym)
if ((nSym /= 1) .and. (nSym /= 2) .and. (nSym /= 4) .and. (nSym /= 8)) then
  rc = rcTwo%TC03
  call SysWarnMsg(TheName,'The file carries an invalid number of irreducible representations',' ')
  call SysValueMsg('nSym',nSym)
end if
nPairs = nTri_Elem(nSym)
iBatch = 0
do iSym=1,nSym
  do jSym=1,iSym
    do kSym=1,nSym
      do lSym=1,kSym
        if (Mul(iSym,jSym) == Mul(kSym,lSym)) then
          ijPair = jSym+nTri_Elem(iSym-1)
          klPair = lSym+nTri_Elem(kSym-1)
          iSyBlk = (ijPair-1)*nPairs+klPair
          iBatch = iBatch+1
          nBatch(iSyBlk) = iBatch
        end if
      end do
    end do
  end do
end do
!----------------------------------------------------------------------*
! Check number of basis function                                       *
!----------------------------------------------------------------------*
ntBas = 0
do iSym=0,(nSym-1)
  nBas(iSym) = TocTwo(isBas+iSym)
  ntBas = ntBas+nBas(iSym)
  if (nBas(iSym) < 0) then
    call SysWarnMsg(TheName,'Invalid number of basis functions',' ')
    call SysValueWarnMsg('iSym',iSym)
    call SysCondMsg('nBas(iSym) < 0',nBas(iSym),'<',0)
  end if

  if (nBas(iSym) > mxBas) then
    call SysWarnMsg(TheName,'Invalid number of basis functions',' ')
    call SysValueWarnMsg('iSym',iSym)
    call SysCondMsg('nBas(iSym) > mxBas',nBas(iSym),'>',mxBas)
  end if
end do
if (ntBas <= 0) then
  call SysWarnMsg(TheName,'Invalid number of basis functions',' ')
  call SysCondMsg('ntBas <= 0',ntBas,'<=',0)
end if
if (ntBas > mxOrb) then
  call SysWarnMsg(TheName,'Invalid number of basis functions',' ')
  call SysCondMsg('ntBas > mxOrb',ntBas,'>',mxOrb)
end if
!----------------------------------------------------------------------*
! Check the skip parameters                                            *
!----------------------------------------------------------------------*
do iSym=0,(nSym-1)
  nSkip(iSym) = TocTwo(isSkip+iSym)
  if (nSkip(iSym) < 0) call SysAbendMsg(TheName,'The table of skiping parameters is spoiled',' ')
end do
!----------------------------------------------------------------------*
! Check table of disk addresses                                         *
!----------------------------------------------------------------------*
mxDAdr = TocTwo(isMxDa)
if (mxDAdr < 0) then
  call SysWarnMsg(TheName,'The file carries an invalid disk address',' ')
  call SysCondMsg('mxDAdr < 0',mxDAdr,'<',0)
end if
do iTab=0,175
  if ((TocTwo(isDAdr+iTab) < 0) .or. (TocTwo(isDAdr+iTab) > mxDAdr)) then
    call SysWarnMsg(TheName,'The table of disk addresses is spoiled',' ')
    call SysValueWarnMsg('iTab',iTab)
    call SysCondMsg('TocTwo(isDAdr+iTab) > mxDAdr',TocTwo(isDAdr+iTab),'>',mxDAdr)
  end if
end do
!----------------------------------------------------------------------*
! Generate and check packing table                                     *
!----------------------------------------------------------------------*
call Int2Real(TocTwo(isPkTh),PkThrs)
if (PkThrs < 0) call SysAbendMsg(TheName,'The accuracy threshold for unpacking is spoiled',' ')
iPack = TocTwo(isPkPa)
if ((iPack < 0) .or. (iPack > 1)) then
  call SysWarnMsg(TheName,'The packing flag is spoiled',' ')
  call SysValueMsg('iPack',iPack)
end if
isPack = TocTwo(isPkPa) == 0

!----------------------------------------------------------------------*
! exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine GetOrd
