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
! Copyright (C) 1993,1996,1998, Markus P. Fuelscher                    *
!               1993, Per Ake Malmqvist                                *
!               1998, Roland Lindh                                     *
!***********************************************************************

subroutine SORT3(MaxDax)
!***********************************************************************
!                                                                      *
!     Purpose: Up to this point the two electron integral file         *
!              has been used like a random stack of records. To        *
!              save hunting for records when reading the records       *
!              bring them into sequential order.                       *
!                                                                      *
!     Called from: Seward_main                                         *
!                                                                      *
!     Calls to : DaFile,DCopy,MkOrd,ClsOrd,ErrOrd                      *
!                                                                      *
!     Calling parameters:                                              *
!     MaxDax  : Higest disk adress of the final 2el integral file      *
!                                                                      *
!     local data declarations:                                         *
!     Buf    : I/O buffer contains packed integral values              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher and P.-AA. Malmqvist                             *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!     - modified to use a virtual disk                                 *
!       M. P. Fuelscher, University of Lund, Sweden, 1996              *
!     - modified to get rid of all disk address calculations           *
!       M. P. Fuelscher, University of Lund, Sweden, 1998              *
!     - modified to eliminate copy statement                           *
!       R. Lindh, University of Lund, Sweden, 1998                     *
!     - modified to eliminate free records                             *
!       R. Lindh, University of Lund, Sweden, 1998                     *
!                                                                      *
!***********************************************************************

use TwoDat, only: lStRec
use sort_data, only: iDaTw0, iDIBin, iDVBin, iStBin, lSll, LuTmp, LuTwo, mInt, MxOrd, n_Int, nBin, nRec, nSln
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: MaxDax
integer(kind=iwp) :: i, iB1, iB2, iBin, iDisk, iDummy, iOpt, iOrd, iRc, iRd, iTmp, iWr, j, j1, j2
real(kind=wp) :: Buf(2*lStRec)
integer(kind=iwp), allocatable :: SrtKey(:), SrtAdr(:)
#ifdef _DEBUGPRINT_
#include "print.fh"
integer(kind=iwp) :: iPrint, iRout
#endif

!----------------------------------------------------------------------*
!     pick up the print level                                          *
!----------------------------------------------------------------------*

#ifdef _DEBUGPRINT_
iRout = 88
iPrint = nPrint(iRout)
if (iPrint > 5) write(u6,*) ' >>> Enter SORT3 <<<'
#endif

!----------------------------------------------------------------------*
!     Scan once the two-electron integral file a pick up the sort      *
!     key as well as disk addresses                                    *
!----------------------------------------------------------------------*

call mma_allocate(SrtKey,MxOrd,Label='SrtKey')
call mma_allocate(SrtAdr,MxOrd,Label='SrtAdr')
iOpt = 2
iDisk = iDaTw0
MaxDax = 0
do iOrd=1,MxOrd
  SrtAdr(iOrd) = iDisk
  MaxDax = max(iDisk,MaxDax)
  call dDAFILE(LuTwo,iOpt,Buf,lStRec,iDisk)
  SrtKey(iOrd) = int(Buf(2),kind=iwp)
end do
MaxDax = iDisk
#ifdef _DEBUGPRINT_
if (iPrint >= 10) then
  call iVcPrt('Sort keys',' ',SrtKey,MxOrd)
  call iVcPrt('Disk addresses',' ',SRtAdr,MxOrd)
end if
#endif

!----------------------------------------------------------------------*
!     Sort records in ascending order of the sort key                  *
!----------------------------------------------------------------------*

iWr = 1
iRd = 2

! Loop over all records

do i=1,MxOrd
  j1 = i
  j2 = SrtKey(i)
  if (j2 /= i) then
    iB1 = 1
    iB2 = lStRec+1
    iDisk = SrtAdr(j1)
    call dDAFILE(LuTwo,iRd,Buf(iB1),lStRec,iDisk)
    do while (j2 /= i)
      iDisk = SrtAdr(j2)
      call dDAFILE(LuTwo,iRd,Buf(iB2),lStRec,iDisk)
      iDisk = SrtAdr(j2)
      call dDAFILE(LuTwo,iWr,Buf(iB1),lStRec,iDisk)
      iTmp = iB2
      iB2 = iB1
      iB1 = iTmp
      j1 = j2
      j2 = SrtKey(j2)
      SrtKey(j1) = j1
    end do
    iDisk = SrtAdr(j2)
    call dDAFILE(LuTwo,iWr,Buf(iB1),lStRec,iDisk)
    SrtKey(j2) = j2
  end if
end do
#ifdef _DEBUGPRINT_
if (iPrint >= 10) then
  call iVcPrt('Sort keys',' ',SrtKey,MxOrd)
end if
#endif

!----------------------------------------------------------------------*
!     Update the disk start adressed of each slice                     *
!----------------------------------------------------------------------*

j = 1
do iBin=1,nBin
  iDVBin(2,iBin) = SrtAdr(j)
  j = j+nRec(iBin)
end do
call mma_deallocate(SrtKey)
call mma_deallocate(SrtAdr)

!----------------------------------------------------------------------*
!     Write the final table of content to disk                         *
!----------------------------------------------------------------------*

call MkOrd(iDummy)

!----------------------------------------------------------------------*
!     Close the ordered 2el integral file                              *
!----------------------------------------------------------------------*

iRc = -1
call ClsOrd(iRc)
if (iRc /= 0) then
  write(u6,*) 'SORT3: Error closing ORDINT'
  call Abend()
end if
call DaClos(LuTmp)

!----------------------------------------------------------------------*
!     Release memory                                                   *
!----------------------------------------------------------------------*

call mma_deallocate(iDIBin)
call mma_deallocate(iDVBin)
call mma_deallocate(mInt)
call mma_deallocate(nRec)
call mma_deallocate(n_Int)
call mma_deallocate(iStBin)
call mma_deallocate(nSln)
call mma_deallocate(lSll)

!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

return

end subroutine SORT3
