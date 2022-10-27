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
!    Global data declarations (Include files) :                        *
!    TwoDat : table of contents and auxiliary information              *
!    TowRc  : Table of return code                                     *
!                                                                      *
!    Local data declarations:                                          *
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

use Symmetry_Info, only: Mul

implicit integer(A-Z)
#include "Molcas.fh"
#include "TwoDat.fh"
real*8 Buf(*)
logical Square

!----------------------------------------------------------------------*
! Start the procedure                                                  *
!----------------------------------------------------------------------*
rc = rc0000
!----------------------------------------------------------------------*
! Pick up the file definitions                                         *
!----------------------------------------------------------------------*
open = AuxTwo(isStat)
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
if (open /= 1) then
  rc = rcRD08
  write(6,*) 'RdOrd: ORDINT not opened yet!'
  call Abend()
end if
!----------------------------------------------------------------------*
! Check if the packing table has been loaded                           *
!----------------------------------------------------------------------*
if ((TocTwo(isPkPa) < 0) .or. (TocTwo(isPkPa) > 1) .or. (TocTwo(isPkAs) < 0) .or. (TocTwo(isPkAs) > 1)) then
  rc = rcRD09
  write(6,*) 'RdOrd: the packing flags are spoiled'
  call Abend()
end if
!----------------------------------------------------------------------*
! Check the symmetry labels                                            *
!----------------------------------------------------------------------*
Square = TocTwo(isOrd) == 1
if (Mul(iSym,jSym) /= Mul(kSym,lSym)) then
  rc = rcRD01
  write(6,*) 'RdOrd: Wrong symmetry labels, direct product is not total symmetric'
  call Abend()
end if
if ((iSym < jSym) .or. (kSym < lSym)) then
  rc = rcRD02
  write(6,*) 'RdOrd: invalid order of symmetry labels'
  call Abend()
end if
ijS = jSym+iSym*(iSym-1)/2
klS = lSym+kSym*(kSym-1)/2
if ((ijS < klS) .and. (.not. Square)) then
  rc = rcRD03
  write(6,*) 'RdOrd: invalid combination of symmetry labels'
  call Abend()
end if
nSym = TocTwo(isSym)
nPairs = nSym*(nSym+1)/2
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
  rc = rcRD07
  write(6,*) 'RdOrd: Requested symmetry block has not been computed'
  call Abend()
end if
!----------------------------------------------------------------------*
! Check options                                                        *
!----------------------------------------------------------------------*
if ((iOpt /= 1) .and. (iOpt /= 2)) then
  rc = rcRD06
  write(6,*) 'RdOrd: Invalid option'
  write(6,*) 'iOpt=',iOpt
  call Abend()
end if
!----------------------------------------------------------------------*
! Check the buffer size                                                *
!----------------------------------------------------------------------*
if (lBuf <= 0) then
  rc = rcRD04
  write(6,*) 'RdOrd: invalid buffer size'
  write(6,*) 'lbuf=',lBuf
  call Abend()
end if
!----------------------------------------------------------------------*
! Compute matrix dimensions                                            *
!----------------------------------------------------------------------*
iB = TocTwo(isBas+iSym-1)
jB = TocTwo(isBas+jSym-1)
kB = TocTwo(isBas+kSym-1)
lB = TocTwo(isBas+lSym-1)
ijB = iB*jB
if (iSym == jSym) ijB = jB*(iB+1)/2
klB = kB*lB
if (kSym == lSym) klB = lB*(kB+1)/2
!----------------------------------------------------------------------*
! Compute submatrix dimensions                                         *
!----------------------------------------------------------------------*
if (lBuf <= 0) then
  rc = rcRD04
  write(6,*) 'RdOrd: invalid buffer size'
  write(6,*) 'lbuf=',lBuf
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
  rc = rcRD05
  write(6,*) 'RdOrd: too small buffer'
  write(6,*) 'Buffer size is lBuf  =',lBuf
  write(6,*) 'Size of submatrix klB=',klB
  write(6,*) 'Call parameters to rdord are:'
  write(6,*) 'iOpt=',iOpt
  write(6,*) 'iSym=',iSym
  write(6,*) 'jSym=',jSym
  write(6,*) 'kSym=',kSym
  write(6,*) 'lSym=',lSym
  write(6,*) 'lBuf=',lBuf
  write(6,*) 'nMat=',nMat
  write(6,*) 'Symmetry block iSyBlk=',iSyBlk
  write(6,*) 'Batch nr       iBatch=',iBatch
  write(6,*) 'iB=TocTwo(isBas+iSym-1), etc:'
  write(6,*) 'iB=',iB
  write(6,*) 'jB=',jB
  write(6,*) 'kB=',kB
  write(6,*) 'lB=',lB
  call Abend()
end if
!----------------------------------------------------------------------*
! Check that the number of submatrices do not run beyond               *
! the last integral of a symmetry allowed batch                        *
!----------------------------------------------------------------------*
Leftpq = AuxTwo(isNpq)
if (iOpt == 1) then
  Leftpq = ijB
else
  nMat = min(Leftpq,nMat)
end if
Leftpq = Leftpq-nMat
AuxTwo(isNpq) = Leftpq
nInts = nMat*klB
!----------------------------------------------------------------------*
! Transfer integrals                                                   *
!----------------------------------------------------------------------*
if (RAMD) then
  call ORDIN2(iOpt,Buf,nInts,iBatch)
else
  call ORDIN1(iOpt,Buf,nInts,iBatch)
end if

!----------------------------------------------------------------------*
! exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine RdOrd_
