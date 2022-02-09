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
!              the bin number of the first integral in a given         *
!              symmetry block.                                         *
!                                                                      *
!     Called from: Sort0                                               *
!                                                                      *
!     Calls to : none                                                  *
!                                                                      *
!     Calling parameters: none                                         *
!                                                                      *
!*** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****

use sort_data, only: iDIBin, iDVBin, iStBin, lBin, lSll, mInt, mSyBlk, mxSyP, n_Int, nBin, nBs, nSkip, nSln, nSyOp, nRec, Square
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6, RtoI

implicit none
#include "print.fh"
#include "warnings.h"
integer(kind=iwp) :: ib, ibj, ij, iOff, iPrint, iRout, iSkip, iSyblj, iSyBlk, iSymi, iSymj, jb, jSkip, jSymj, kb, kbl, kSkip, &
                     kSybll, kSymk, kSyml, kSymMx, lb, lSkip, lSlice, lSrtA, lSyml, lSymMx, MaxMem, mnBin, MxSrtA1, MxSrtA2, nij, &
                     nSlice, nSym
#ifndef _I8_
integer(kind=iwp) :: ix
integer(kind=iwp), parameter :: lim_32 = 2**30
#endif
logical(kind=iwp), external:: Reduce_Prt

iRout = 80
iPrint = nPrint(iRout)
if (iPrint > 10) write(u6,*) ' >>> Enter MKSRT1 <<<'
!----------------------------------------------------------------------*
!                                                                      *
!     grab memory                                                      *
!                                                                      *
!     The available memory should not be smaller than the number       *
!     of unique symmetry blocks times the size of a bin.               *
!     If nSyOp=1 mnBin=  1 and If Square=.true. mnBin=  1              *
!        nSyOp=2 mnBin=  4                      mnBin=  5              *
!        nSyOp=4 mnBin= 19                      mnBin= 28              *
!        nSyOp=8 mnBin=106                      mnBin=176              *
!                                                                      *
!----------------------------------------------------------------------*

select case (nSyOp)
  case (2)
    mnBin = 4
    if (Square) mnBin = 5
  case (4)
    mnBin = 19
    if (Square) mnBin = 28
  case (8)
    mnBin = 106
    if (Square) mnBin = 176
  case default
    mnBin = 1
end select
lSrtA = mnBin*lBin
lSrtA = ((1+RtoI)*lSrtA)/RtoI
call mma_MaxDBLE(MaxMem)
lSrtA = max(lSrtA,MaxMem/2)

!----------------------------------------------------------------------*
!     determine the partitioning of the 2el integrals into             *
!     submatrices by dividing the symmetry blocks into slices that     *
!     fit into the available memory.                                   *
!----------------------------------------------------------------------*

mSyBlk = mxSyP**2
call mma_allocate(iStBin,mSyBlk,label='iStBin')
call mma_allocate(nSln,mSyBlk,label='nSln')
call mma_allocate(lSll,mSyBlk,label='lSll')

nSln(:) = 0
lSll(:) = 0

nSym = nSyOp
do iSymi=1,nSym
  ib = nBs(iSymi)
  iSkip = nSkip(iSymi)
  if ((ib == 0) .or. (iSkip == 1)) cycle
  do jSymj=1,iSymi
    iSymj = 1+ieor(iSymi-1,jSymj-1)
    iSyblj = jSymj+iSymi*(iSymi-1)/2
    jb = nBs(jSymj)
    ibj = ib*jb
    if (iSymi == jSymj) ibj = ib*(ib+1)/2
    jSkip = nSkip(jSymj)
    if ((jb == 0) .or. (jSkip == 1)) cycle
    kSymMx = iSymi
    if (Square) kSymMx = nSym
    do kSymk=1,kSymMx
      kb = nBs(kSymk)
      kSkip = nSkip(kSymk)
      if ((kb == 0) .or. (kSkip == 1)) cycle
      lSymMx = jSymj
      if ((kSymk /= iSymi) .or. Square) lSymMx = kSymk
      do lSyml=1,lSymMx
        kSyml = 1+ieor(kSymk-1,lSyml-1)
        kSybll = lSyml+kSymk*(kSymk-1)/2
        if (ieor(iSymj-1,kSyml-1) /= 0) cycle

        !write(u6,*) 'i,j,k,l=',iSymi,jSymj,kSymk,lSyml
        lb = nBs(lSyml)
        lSkip = nSkip(lSyml)
        if ((lb == 0) .or. (lSkip == 1)) cycle
        kbl = kb*lb
#       ifndef _I8_
        ix = lim_32/kbl
#       endif
        if (kSymk == lSyml) kbl = kb*(kb+1)/2
        iSyBlk = kSybll+mxSyP*(iSyblj-1)
        nij = 0
        do ij=1,ibj
#         ifndef _I8_
          if (ij <= ix) then
#         endif
            if ((ij*kbl) < lSrtA) nij = ij
#         ifndef _I8_
          end if
#         endif
        end do
        if (nij == 0) then
          write(u6,*)
          write(u6,'(2X,A,I3.3,A)') '*** Error in MKSRT1 ***'
          write(u6,'(2X,A)') 'Insufficient memory'
          write(u6,'(2X,A)') 'Increase the value of the MOLCAS_MEM environement variable'
          write(u6,*)
          call Quit(_RC_MEMORY_ERROR_)
        end if
        nSlice = 1+(ibj-1)/nij
        do ij=nij,1,-1
          if ((ij*nSlice) >= ibj) nij = ij
        end do
        nSln(iSyBlk) = nSlice
        lSll(iSyBlk) = nij*kbl
        !write(u6,*) 'iSyBlk=',iSyBlk
        !write(u6,*) 'nSln(iSyBlk)=',nSln(iSyBlk)
        !write(u6,*) 'lSll(iSyBlk)=',lSll(iSyBlk)
      end do
    end do
  end do
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

if ((.not. Reduce_Prt()) .and. (iPrint > 5)) then
  write(u6,*)
  write(u6,'(A,I12,A,I4,A)') '  SEWARD will use a sorting area of',MxSrtA1,' Words(Real*8) in the first phase (=',nBin,' bins).'
  write(u6,'(A,I12,A)') '  SEWARD will use a sorting area of',MxSrtA2,' Words(Real*8) in the second phase.'
  write(u6,*)
end if

call mma_allocate(iDIBin,3,max(mnBin,nBin),label='iDIBin')
call mma_allocate(iDVBin,4,max(mnBin,nBin),label='iDVBin')
call mma_allocate(mInt,3,max(mnBin,nBin),label='mInt')
call mma_allocate(nRec,max(mnBin,nBin),label='nRec')
call mma_allocate(n_Int,max(mnBin,nBin),label='n_Int')
iDIBin(:,:) = 0
iDVBin(:,:) = 0
mInt(:,:) = 0
nRec(:) = 0
n_Int(:) = 0

if (MxSrtA1 > lSrtA) then
  write(u6,*)
  write(u6,'(2X,A,I3.3,A)') '*** Error in MKSRT1 ***'
  write(u6,'(2X,A)') 'Insufficient memory'
  write(u6,'(2X,A)') 'MxSrtA1>lSrtA'
  write(u6,*)
  write(u6,*) 'MxSrtA1=',MxSrtA1
  write(u6,*) 'lSrtA=',lSrtA
  write(u6,*)
  write(u6,*) 'Increase MOLCAS_MEM and try again!'
  write(u6,*)
  call Quit(_RC_MEMORY_ERROR_)
end if

!----------------------------------------------------------------------*
!     compute offsets                                                  *
!----------------------------------------------------------------------*

iOff = 1
do iSyBlk=1,mSyBlk
  nSlice = nSln(iSyBlk)
  iStBin(iSyBlk) = iOff
  iOff = iOff+nSlice
end do

return

end subroutine MkSrt1
