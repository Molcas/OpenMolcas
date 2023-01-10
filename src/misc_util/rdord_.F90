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

subroutine RdOrd_(rc,iOpt,iSym,jSym,kSym,lSym,Buf,lBuf,nMat)
!***********************************************************************
!                                                                      *
!    Purpose: Read a buffer of ordered integrals of length lBuf        *
!             The subroutine returns the number of submatrices         *
!             it has picked up.                                        *
!                                                                      *
!    Calling parameters:                                               *
!    iSym   : irred. representation of first symmetry label            *
!    jSym   : irred. representation of second symmetry label           *
!    kSym   : irred. representation of third symmetry label            *
!    lSym   : irred. representation of fourth symmetry label           *
!    Buf    : containts on output the integrals                        *
!    lBuf   : length of the integral buffer                            *
!    nMat   : number of submatrices read in                            *
!    iOpt   : option code (iOpt=1:start reading at first integral)     *
!                         (iOpt=2:continue reading)                    *
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
use TwoDat, only: AuxTwo, isBas, isOrd, isPkPa, isSkip, isSym, nBatch, RAMD, rcTwo, TocTwo
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(out) :: rc, nMat
integer(kind=iwp), intent(in) :: iOpt, iSym, jSym, kSym, lSym, lBuf
real(kind=wp), intent(_OUT_) :: Buf(*)
integer(kind=iwp) :: iB, iBatch, ijB, ijS, iSkip, iSyBlk, jB, jSkip, kB, klB, klS, kSkip, lB, Leftpq, lSkip, nInts, nPairs, nSym
logical(kind=iwp) :: Square

!----------------------------------------------------------------------*
! Start the procedure                                                  *
!----------------------------------------------------------------------*
rc = rcTwo%good
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
if (.not. AuxTwo%Opn) then
  rc = rcTwo%RD08
  write(u6,*) 'RdOrd: ORDINT not opened yet!'
  call Abend()
end if
!----------------------------------------------------------------------*
! Check if the packing table has been loaded                           *
!----------------------------------------------------------------------*
if ((TocTwo(isPkPa) < 0) .or. (TocTwo(isPkPa) > 1)) then
  rc = rcTwo%RD09
  write(u6,*) 'RdOrd: the packing flags are spoiled'
  call Abend()
end if
!----------------------------------------------------------------------*
! Check the symmetry labels                                            *
!----------------------------------------------------------------------*
Square = TocTwo(isOrd) == 1
if (Mul(iSym,jSym) /= Mul(kSym,lSym)) then
  rc = rcTwo%RD01
  write(u6,*) 'RdOrd: Wrong symmetry labels, direct product is not total symmetric'
  call Abend()
end if
if ((iSym < jSym) .or. (kSym < lSym)) then
  rc = rcTwo%RD02
  write(u6,*) 'RdOrd: invalid order of symmetry labels'
  call Abend()
end if
ijS = jSym+nTri_Elem(iSym-1)
klS = lSym+nTri_Elem(kSym-1)
if ((ijS < klS) .and. (.not. Square)) then
  rc = rcTwo%RD03
  write(u6,*) 'RdOrd: invalid combination of symmetry labels'
  call Abend()
end if
nSym = TocTwo(isSym)
nPairs = nTri_Elem(nSym)
iSyBlk = (ijS-1)*nPairs+klS
iBatch = nBatch(iSyBlk)
!----------------------------------------------------------------------*
! Check the skipping flags                                             *
!----------------------------------------------------------------------*
iSkip = TocTwo(isSkip+iSym-1)
jSkip = TocTwo(isSkip+jSym-1)
kSkip = TocTwo(isSkip+kSym-1)
lSkip = TocTwo(isSkip+lSym-1)
if ((iSkip+jSkip+kSkip+lSkip) /= 0) then
  rc = rcTwo%RD07
  write(u6,*) 'RdOrd: Requested symmetry block has not been computed'
  call Abend()
end if
!----------------------------------------------------------------------*
! Check options                                                        *
!----------------------------------------------------------------------*
if ((iOpt /= 1) .and. (iOpt /= 2)) then
  rc = rcTwo%RD06
  write(u6,*) 'RdOrd: Invalid option'
  write(u6,*) 'iOpt=',iOpt
  call Abend()
end if
!----------------------------------------------------------------------*
! Check the buffer size                                                *
!----------------------------------------------------------------------*
if (lBuf <= 0) then
  rc = rcTwo%RD04
  write(u6,*) 'RdOrd: invalid buffer size'
  write(u6,*) 'lbuf=',lBuf
  call Abend()
end if
!----------------------------------------------------------------------*
! Compute matrix dimensions                                            *
!----------------------------------------------------------------------*
iB = TocTwo(isBas+iSym-1)
jB = TocTwo(isBas+jSym-1)
kB = TocTwo(isBas+kSym-1)
lB = TocTwo(isBas+lSym-1)
if (iSym == jSym) then
  ijB = nTri_Elem(iB)
else
  ijB = iB*jB
end if
if (kSym == lSym) then
  klB = nTri_Elem(kB)
else
  klB = kB*lB
end if
!----------------------------------------------------------------------*
! Compute submatrix dimensions                                         *
!----------------------------------------------------------------------*
if (lBuf <= 0) then
  rc = rcTwo%RD04
  write(u6,*) 'RdOrd: invalid buffer size'
  write(u6,*) 'lbuf=',lBuf
  call Abend()
end if
if (klB <= 0) then
  nMat = 0
  return
else
  nMat = (lBuf-1)/klB
end if
if (nMat > ijB) nMat = ijB
if (nMat == 0) then
  rc = rcTwo%RD05
  write(u6,*) 'RdOrd: too small buffer'
  write(u6,*) 'Buffer size is lBuf  =',lBuf
  write(u6,*) 'Size of submatrix klB=',klB
  write(u6,*) 'Call parameters to rdord are:'
  write(u6,*) 'iOpt=',iOpt
  write(u6,*) 'iSym=',iSym
  write(u6,*) 'jSym=',jSym
  write(u6,*) 'kSym=',kSym
  write(u6,*) 'lSym=',lSym
  write(u6,*) 'lBuf=',lBuf
  write(u6,*) 'nMat=',nMat
  write(u6,*) 'Symmetry block iSyBlk=',iSyBlk
  write(u6,*) 'Batch nr       iBatch=',iBatch
  write(u6,*) 'iB=TocTwo(isBas+iSym-1), etc:'
  write(u6,*) 'iB=',iB
  write(u6,*) 'jB=',jB
  write(u6,*) 'kB=',kB
  write(u6,*) 'lB=',lB
  call Abend()
end if
!----------------------------------------------------------------------*
! Check that the number of submatrices do not run beyond               *
! the last integral of a symmetry allowed batch                        *
!----------------------------------------------------------------------*
Leftpq = AuxTwo%Npq
if (iOpt == 1) then
  Leftpq = ijB
else
  nMat = min(Leftpq,nMat)
end if
Leftpq = Leftpq-nMat
AuxTwo%Npq = Leftpq
nInts = nMat*klB
!----------------------------------------------------------------------*
! Transfer integrals                                                   *
!----------------------------------------------------------------------*
if (RAMD%act) then
  call ORDIN2(iOpt,Buf,nInts,iBatch)
else
  call ORDIN1(iOpt,Buf,nInts,iBatch)
end if

!----------------------------------------------------------------------*
! exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine RdOrd_
