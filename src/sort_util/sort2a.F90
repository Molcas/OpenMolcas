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

subroutine SORT2A(iBin,lSrtA,SrtArr,IOStk,lStk,nStk)
!***********************************************************************
!                                                                      *
!     Purpose: Reload all integral belonging to a given slice          *
!              and sort them.                                          *
!                                                                      *
!     Called from: Sort2                                               *
!                                                                      *
!     Calls to : UPkI4,UPkR8,DaFile                                    *
!                                                                      *
!     Calling parameters:                                              *
!     iBin   : Slice number                                            *
!     SrtArr : Work space to keep the 2el integrals                    *
!                                                                      *
!     Global data declarations (Include files) :                       *
!     TwoDef  : definitions of the record structure                    *
!     Srt0    : common block containing information pertinent to       *
!               the calculation of 2el integral sequence numbers       *
!     Srt1    : common block containing information the number of      *
!               bins and partitioning of symmetry blocks               *
!     Srt2    : common block containing information pertinent to       *
!               the bin sorting algorithm                              *
!     Srt3    : dynamic stack to control inout and output of           *
!               integral buffers                                       *
!                                                                      *
!     local data declarations:                                         *
!     PkVBin : I/O buffer contains packed integral values              *
!     PkIBin : I/O buffer contains packed ordering numbers             *
!     ValBin : contains unpacked integral values                       *
!     IndBin : contains unpacked ordering numbers                      *
!                                                                      *
!*** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****

use srt2
implicit real*8(A-H,O-Z)

#include "TwoDat.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "SysDef.fh"
#include "print.fh"
#include "PkCtl.fh"
#include "warnings.fh"

dimension SrtArr(lSrtA)
dimension PkVBin(lStRec)
integer PkIBin(lStRec)
integer IOStk(lStk)

!----------------------------------------------------------------------*
!     as the packed integral labels add an extra 1-2 Byte              *
!     disk space per integral we have to adjust the record             *
!     length of LuTmp to the different machines.                       *
!----------------------------------------------------------------------*

idiv = ItoB/2
if (Pack) idiv = idiv/2
mStRec = (lStRec/idiv)
mDaRec = (lDaRec/idiv)

!----------------------------------------------------------------------*
!     pick up the print level                                          *
!----------------------------------------------------------------------*

iRout = 85
iPrint = nPrint(iRout)
if (iPrint >= 10) then
  write(6,*) ' >>> Enter SORT2A <<<'
  write(6,*) ' iBin  ',iBin
  write(6,*) ' lSrtA ',lSrtA
end if

iZero = lSrtA-mInt(1,iBin)
iInt = mInt(2,iBin)*RtoB
iInd = mInt(3,iBin)
iP_Storage = (iZero+iInt+RtoB)/RtoB
iI_Storage = (iInd+iInt+RtoB)/RtoB
if (iP_Storage <= iI_Storage) then
  iDVBin(4,iBin) = 0  ! Dense mode.
  !write(6,*) 'Mode: Dense'
else
  iDVBin(4,iBin) = 1  ! Sparse mode.
  !write(6,*) 'Mode: Sparse'
end if

#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'Processing slice                   :',iBin
write(6,*) 'Actual number of non-zero integrals:',mInt(1,iBin)
write(6,*) 'Effective number of integrals      :',mInt(2,iBin)
write(6,*) 'Effective number of indicies       :',mInt(3,iBin)
write(6,*) 'Total number of integrals          :',lSrtA
write(6,*) 'Packed storage                     :',iP_Storage
write(6,*) 'Indexed storage                    :',iI_Storage
#endif

!----------------------------------------------------------------------*
!     Start reading packed buffers                                     *
!----------------------------------------------------------------------*

iDaTmp = iDIBin(2,iBin)
iDaTwo = iDVBin(2,iBin)
do while (iDaTmp >= 0)
  nStk = nStk+1
  if (nStk > lStk) then
    write(6,*)
    write(6,'(2X,A,I3.3,A)') '*** Error in SORT2A ***'
    write(6,'(2X,A)') 'nStk exceeds limits (nStk>lStk)'
    write(6,'(2X,A,I8)') 'nStk =',nStk
    write(6,'(2X,A,I8)') 'lStk =',lStk
    write(6,'(2X,A,I8)') 'iBin =',iBin
    write(6,*)
    write(6,*) 'Action: rerun with a larger MOLCAS_MEM'
    call Quit(_RC_MEMORY_ERROR_)
  end if
  IOStk(nStk) = iDaTwo
  iOpt = 2
  if (iPrint >= 10) then
    write(6,*) ' read records: iDaTmp,iDaTwo ',iDaTmp,iDaTwo
  end if
  call iDAFILE(LuTmp,iOpt,PkIBin,mStRec,iDaTmp)
  call dDAFILE(LuTwo,iOpt,PkVBin,lStRec,iDaTwo)
  ist1 = lTop
  ist2 = lTop
  do iSec=1,nSect
    nInts1 = PkIBin(ist1-1)
    nInts2 = int(PkVBin(ist2-1))
    if (nInts1 /= nInts2) then
      write(6,*)
      write(6,'(2X,A,I3.3,A)') '*** Error in SORT2A ***'
      write(6,'(2X,A)') 'An inconsistency has been deteced'
      write(6,'(2X,A)') 'nInts1#nInts2'
      write(6,*)
      call xFlush(6)
      call Abend
    end if
    nInts = nInts1
    if (nInts > lBin) then
      write(6,*)
      write(6,'(2X,A,I3.3,A)') '*** Error in SORT2A ***'
      write(6,'(2X,A)') 'An inconsistency has been deteced'
      write(6,'(2X,A)') 'nInts>lBin'
      write(6,*)
      call xFlush(6)
      call Abend
    end if
    if (nInts > 0) then

      !----------------------------------------------------------------*
      ! Unpack buffers                                                 *
      !----------------------------------------------------------------*

      call UPKI4(nInts,lIBin,PkIBin(ist1+1),IndBin)
      iOpt = 0 ! Always tight mode
      call UPKR8(iOpt,nInts,lVBin,PkVBin(ist2+1),ValBin)

      !----------------------------------------------------------------*
      ! Sort 2el integrals                                             *
      !----------------------------------------------------------------*

      do indx=1,nInts
        SrtArr(IndBin(indx)) = ValBin(indx)
      end do
      ist1 = ist1+mDaRec
      ist2 = ist2+lDaRec
    end if
  end do

  !--------------------------------------------------------------------*
  !   Get the disk adress of the next record                           *
  !--------------------------------------------------------------------*

  iDaTmp = PkIBin(1)
  iDaTwo = int(PkVBin(1))
end do
if (iPrint >= 99) call dVcPrt('sorted ERIs',' ',SrtArr,lSrtA)

return

end subroutine SORT2A
