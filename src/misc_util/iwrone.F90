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
! Copyright (C) 1993, Per-Olof Widmark                                 *
!               1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine iWrOne(rc,Option,InLab,Comp,data,SymLab)
!***********************************************************************
!                                                                      *
!     Purpose: write data to one-electron integral file                *
!                                                                      *
!     Calling parameters:                                              *
!     rc      : Return code                                            *
!     Option  : Switch to set options                                  *
!     InLab   : Identifier for the data to write                       *
!     Comp    : Composite identifier to select components              *
!     Data    : contains on input the data to store on disk            *
!     SymLab  : symmetry label of the provided data                    *
!                                                                      *
!     Global data declarations (Include file):                         *
!     Parm    : Table of paramaters                                    *
!     rcParm  : Table of return codes                                  *
!     Switch  : Table of options                                       *
!     Common  : Common block containing TocOne                         *
!     Data    : Data definitions                                       *
!                                                                      *
!     Local data declarations:                                         *
!     Label   : character*8, used to covert incoming names             *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark and M.P. Fuelscher                                  *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)
#include "OneDat.fh"
character*(*) InLab
dimension data(*)
character*8 Label
dimension LabTmp(2)
logical debug, close
! Statement function
MulTab(i,j) = ieor(i-1,j-1)+1

!----------------------------------------------------------------------*
! Start procedure:                                                     *
! Pick up the file definitions                                         *
!----------------------------------------------------------------------*
rc = rc0000
LuOne = AuxOne(pLu)
open = AuxOne(pOpen)
!----------------------------------------------------------------------*
! Check the file status                                                *
!----------------------------------------------------------------------*
close = .false.
if (open /= 1) then

  ! Well, I'll open and close it for you under the default name

  LuOne = 77
  LuOne = isFreeUnit(LuOne)
  Label = 'ONEINT  '
  !write(6,*) 'WrOne: opening OneInt'
  iRC = -1
  iOpt = 0
  call OpnOne(iRC,iOpt,Label,LuOne)
  if (iRC /= 0) then
    write(6,*) 'WrOne: Error opening file'
    call Abend()
  end if
  close = .true.
end if
!----------------------------------------------------------------------*
! Truncate the label to 8 characters and convert it to upper case      *
!----------------------------------------------------------------------*
Label = InLab
call UpCase(Label)
Length = len(Label)/ItoB
LabTmp(:Length) = transfer(Label,LabTmp,Length)
!----------------------------------------------------------------------*
! Print debugging information                                          *
!----------------------------------------------------------------------*
debug = .false.
if (iand(option,1024) /= 0) debug = .true.
if (debug) then
  call DmpOne()
  write(6,*) '<<< Entering WrOne >>>'
  write(6,'(a,z8)') ' rc on entry:     ',rc
  write(6,'(a,a)') ' Label on entry:  ',Label
  write(6,'(a,z8)') ' Comp on entry:   ',Comp
  write(6,'(a,z8)') ' SymLab on entry: ',SymLab
  write(6,'(a,z8)') ' Option on entry: ',Option
end if
!----------------------------------------------------------------------*
! Store operators as individual records on disk                        *
! Note: If the incoming operator has already been stored               *
! previously (label, component and symmetry labels are identical)      *
! it will replace the existing one.                                    *
!----------------------------------------------------------------------*
k = 0
do i=MxOp,1,-1
# ifdef _I8_
  if ((TocOne(pOp+LenOp*(i-1)+oLabel) == LabTmp(1)) .and. (TocOne(pOp+LenOp*(i-1)+oComp) == Comp) .and. &
      (TocOne(pOp+LenOp*(i-1)+oSymLb) == SymLab)) k = i
# else
  if ((TocOne(pOp+LenOp*(i-1)+oLabel) == LabTmp(1)) .and. (TocOne(pOp+LenOp*(i-1)+oLabel+1) == LabTmp(2)) .and. &
      (TocOne(pOp+LenOp*(i-1)+oComp) == Comp) .and. (TocOne(pOp+LenOp*(i-1)+oSymLb) == SymLab)) k = i
# endif
end do
iDisk = TocOne(pOp+LenOp*(k-1)+oAddr)
if (k == 0) then
  do i=MxOp,1,-1
    if (TocOne(pOp+LenOp*(i-1)+oLabel) == NaN) k = i
  end do
  iDisk = TocOne(pNext)
end if
if (k == 0) then
  rc = rcWR11
  write(6,*) 'WrOne: The total number of operators exceeds the limit'
  write(6,*) 'k == 0'
  call Abend()
end if
Length = 0
do i=1,nSym
  do j=1,i
    ij = MulTab(i,j)-1
    if (iand(2**ij,SymLab) /= 0) then
      if (i == j) then
        Length = Length+nBas(i)*(nBas(i)+1)/2
      else
        Length = Length+nBas(i)*nBas(j)
      end if
    end if
  end do
end do
Length = RtoI*(Length+nAuxDt)
TocOne(pOp+LenOp*(k-1)+oLabel) = LabTmp(1)
#ifndef _I8_
TocOne(pOp+LenOp*(k-1)+oLabel+1) = LabTmp(2)
#endif
TocOne(pOp+LenOp*(k-1)+oComp) = Comp
TocOne(pOp+LenOp*(k-1)+oSymLb) = SymLab
TocOne(pOp+LenOp*(k-1)+oAddr) = iDisk
call iDaFile(LuOne,1,data,Length,iDisk)
TocOne(pNext) = max(TocOne(pNext),iDisk)
!----------------------------------------------------------------------*
! Finally copy the TocOne back to disk                                 *
!----------------------------------------------------------------------*
iDisk = 0
call iDaFile(LuOne,1,TocOne,lToc,iDisk)

if (close) then
  iRC = -1
  iOpt = 0
  call ClsOne(iRC,iOpt)
  if (iRC /= 0) then
    write(6,*) 'WrOne: Error closing file'
    call Abend()
  end if
end if

!----------------------------------------------------------------------*
! Terminate procedure                                                  *
!----------------------------------------------------------------------*
return

end subroutine iWrOne
