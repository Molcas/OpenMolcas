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
! Copyright (C) 1991, Markus P. Fuelscher                              *
!               1991, Per Ake Malmqvist                                *
!***********************************************************************

subroutine MkSrt1()
!***********************************************************************
!                                                                      *
!     Purpose: Determine the splitting of the integrals into           *
!              submatrices (slices) and hence determine the            *
!              number of bins. Also construct the offsets              *
!              which allow to reduce an integral symmetry block        *
!              number into a relative bin number with respect to       *
!              the bin number of the first integral in an given        *
!              symmetry block.                                         *
!                                                                      *
!     Called from: Sort0                                               *
!                                                                      *
!     Calls to : none                                                  *
!                                                                      *
!     Calling parameters: none                                         *
!                                                                      *
!     Global data declarations (Include files) :                       *
!     TwoDef  : definitions of the record structure                    *
!     Srt0    : common block containing information pertinent to       *
!               the calculation of 2el integral sequence numbers       *
!     Srt1    : common block containing information the number of      *
!               bins and partitioning of symmetry blocks               *
!     Srt2    : common block containing information pertinent to       *
!               the bin sorting algorithm                              *
!                                                                      *
!     Local data declarations: none                                    *
!                                                                      *
!*** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****

use srt2
implicit integer(A-Z)

#include "srt0.fh"
#include "srt1.fh"

#include "SysDef.fh"
#include "print.fh"
#include "warnings.fh"

iRout = 80
iPrint = nPrint(iRout)
if (iPrint > 10) write(6,*) ' >>> Enter MKSRT1 <<<'
!----------------------------------------------------------------------*
!                                                                      *
!     grab memory                                                      *
!                                                                      *
!     The available memory should not be smaller than the number       *
!     of unique symmetry blocks times the size of a bin.               *
!     If nSyOp=1 mxBin=  1 and If Square=.true. mxBin=  1              *
!        nSyOp=2 mxBin=  4                      mxBin=  5              *
!        nSyOp=4 mxBin= 19                      mxBin= 28              *
!        nSyOp=8 mxBin=106                      mxBin=176              *
!                                                                      *
!----------------------------------------------------------------------*

lSrtA = 1
if (nSyOp == 2) then
  lSrtA = 4
  if (Square) lSrtA = 5
else if (nSyOp == 4) then
  lSrtA = 19
  if (Square) lSrtA = 28
else if (nSyOp == 8) then
  lSrtA = 106
  if (Square) lSrtA = 176
end if
lSrtA = lSrtA*lBin
lSrtA = ((1+RtoI)*lSrtA)/RtoI
call mma_MaxDBLE(MaxMem)
lSrtA = max(lSrtA,MaxMem/2)

!----------------------------------------------------------------------*
!     determine the partitioning of the 2el integrals into             *
!     submatrices by dividing the symmetry blocks into slices that     *
!     fit into the available memory.                                   *
!----------------------------------------------------------------------*

do iSyBlk=1,mSyBlk
  nSln(iSyBlk) = 0
  lSll(iSyBlk) = 0
end do

#ifndef _I8_
lim_32 = 2**30
#endif
nSym = nSyOp
do iSymi=1,nSym
  ib = nBs(iSymi)
  iSkip = nSkip(iSymi)
  if ((ib == 0) .or. (iSkip == 1)) Go To 100
  do jSymj=1,iSymi
    iSymj = 1+ieor(iSymi-1,jSymj-1)
    iSyblj = jSymj+iSymi*(iSymi-1)/2
    jb = nBs(jSymj)
    ibj = ib*jb
    if (iSymi == jSymj) ibj = ib*(ib+1)/2
    jSkip = nSkip(jSymj)
    if ((jb == 0) .or. (jSkip == 1)) Go To 200
    kSymMx = iSymi
    if (Square) kSymMx = nSym
    do kSymk=1,kSymMx
      kb = nBs(kSymk)
      kSkip = nSkip(kSymk)
      if ((kb == 0) .or. (kSkip == 1)) Go To 300
      lSymMx = jSymj
      if ((kSymk /= iSymi) .or. Square) lSymMx = kSymk
      do lSyml=1,lSymMx
        kSyml = 1+ieor(kSymk-1,lSyml-1)
        kSybll = lSyml+kSymk*(kSymk-1)/2
        if (ieor(iSymj-1,kSyml-1) /= 0) Go To 400

        !write(6,*) 'i,j,k,l=',iSymi,jSymj,kSymk,lSyml
        lb = nBs(lSyml)
        lSkip = nSkip(lSyml)
        if ((lb == 0) .or. (lSkip == 1)) Go To 400
        kbl = kb*lb
#       ifndef _I8_
        ix = lim_32/kbl
#       endif
        if (kSymk == lSyml) kbl = kb*(kb+1)/2
        iSyBlk = kSybll+mxSyP*(iSyblj-1)
        nij = 0
        do ij=1,ibj
#         ifdef _I8_
          if ((ij*kbl) < lSrtA) nij = ij
#         else
          if (ij <= ix) then
            if ((ij*kbl) < lSrtA) nij = ij
          end if
#         endif
        end do
        if (nij == 0) then
          write(6,*)
          write(6,'(2X,A,I3.3,A)') '*** Error in MKSRT1 ***'
          write(6,'(2X,A)') 'Insufficient memory'
          write(6,'(2X,A)') 'Increase the value of the MOLCAS_MEM environement variable'
          write(6,*)
          call Quit(_RC_MEMORY_ERROR_)
        end if
        nSlice = 1+(ibj-1)/nij
        do ij=nij,1,-1
          if ((ij*nSlice) >= ibj) nij = ij
        end do
        nSln(iSyBlk) = nSlice
        lSll(iSyBlk) = nij*kbl
        !write(6,*) 'iSyBlk=',iSyBlk
        !write(6,*) 'nSln(iSyBlk)=',nSln(iSyBlk)
        !write(6,*) 'lSll(iSyBlk)=',lSll(iSyBlk)
400     continue
      end do
300   continue
    end do
200 continue
  end do
100 continue
end do

!----------------------------------------------------------------------*
!     Determine the final amount of memory and bins that will be used. *
!     Check again for consistency                                      *
!----------------------------------------------------------------------*

nBin = 0
MxSrtA2 = 0
do iSyBlk=1,mSyBlk
  nSlice = nSln(iSyBlk)
  lSlice = lSll(iSyBlk)
  nBin = nBin+nSlice
  MxSrtA2 = max(MxSrtA2,lSlice)
end do
MxSrtA1 = ((1+RtoI)*nBin*lBin)/RtoI

if (iPrint > 5) then
  write(6,*)
  write(6,'(A,I12,A,I4,A)') '  SEWARD will use a sorting area of',MxSrtA1,' Words(Real*8) in the first phase (=',nBin,' bins).'
  write(6,'(A,I12,A)') '  SEWARD will use a sorting area of',MxSrtA2,' Words(Real*8) in the second phase.'
  write(6,*)
end if

if (nBin > mxBin) then
  write(6,*)
  write(6,'(2X,A,I3.3,A)') '*** Error in MKSRT1 ***'
  write(6,'(2X,A)') 'Insufficient memory'
  write(6,'(2X,A)') 'nBin exceeds limits (nBin>mxBin)'
  write(6,*)
  write(6,*) 'nBin=',nBin
  write(6,*) 'mxBin=',mxBin
  write(6,*)
  write(6,*) 'Increase MOLCAS_MEM and try again!'
  write(6,*)
  call Quit(_RC_MEMORY_ERROR_)
end if

if (MxSrtA1 > lSrtA) then
  write(6,*)
  write(6,'(2X,A,I3.3,A)') '*** Error in MKSRT1 ***'
  write(6,'(2X,A)') 'Insufficient memory'
  write(6,'(2X,A)') 'MxSrtA1>lSrtA'
  write(6,*)
  write(6,*) 'MxSrtA1=',MxSrtA1
  write(6,*) 'lSrtA=',lSrtA
  write(6,*)
  write(6,*) 'Increase MOLCAS_MEM and try again!'
  write(6,*)
  call Quit(_RC_MEMORY_ERROR_)
end if

!----------------------------------------------------------------------*
!     compute offsets                                                  *
!----------------------------------------------------------------------*

iOff = 1
do iSyBlk=1,mSyBlk
  nSlice = nSln(iSyBlk)
  IstBin(iSyBlk) = iOff
  iOff = iOff+nSlice
end do

!----------------------------------------------------------------------*
!     initialize counter for integrals per slice                       *
!----------------------------------------------------------------------*

call ICopy(3*nBin,[0],0,mInt,1)

return

end subroutine MkSrt1
