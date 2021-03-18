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
! Copyright (C) 2017, Valera Veryazov                                  *
!***********************************************************************

subroutine XMatReader(iZMUnit,LuWr,nAtoms,nXAtoms,nBasis,nAskAtoms,nxbas,xb_label,xb_bas,iErr)
! nAtoms : Total number of real atoms. Include X dummy atoms: NAT(i)= 0
! nXAtoms: Total number of ghost (Z) atoms:                   NAT(i)=-1
! nBasis : Nummer of atom types requiring Basis Set
!  ***  nAskAtoms == -1  =>  Seward ZMAT input  => Use "End of"
!  ***  nAskAtoms /= -1  =>  GateWay ZMAT input => Use nAskAtoms

use ZMatConv_Mod, only: BasReq, iZmat, MaxAtoms, NAT, Symbols, Zmat
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iZMUnit, LuWr, nAskAtoms
integer(kind=iwp), intent(out) :: nAtoms, nXAtoms, nBasis, nxbas, iErr
character(len=*), intent(inout) :: xb_label(*), xb_bas(*)
integer(kind=iwp) :: i, IreadHere, iXU, NA, NAtom, Nwords
real(kind=wp) :: Dist
character(len=80) :: Line
character(len=24) :: Words(7)
character(len=3) :: Command

xb_label(1) = ' '
xb_bas(1) = ' '
nxbas = 1

iErr = 0
nAtoms = 0
nXAtoms = 0
nBasis = 0
do i=1,100 ! MaxNat
  BasReq(i) = .false.
end do
do i=1,MaxAtoms
  Symbols(i) = '     '
  NAT(i) = 0
  iZmat(i,:) = 0
  Zmat(i,:) = Zero
end do

! Read Line (or COMMAND)
10 if ((nAtoms+nXAtoms) == nAskAtoms) goto 100
read(iZMUnit,'(A)',Err=9906,end=9999) Line
if (Line(1:1) == '*') goto 10
if (Line == ' ') goto 100
Command = Line(1:3)
call UpCase(Command)
if (Command == 'END') goto 100
iErr = 0
NA = 0
Dist = Zero

! Here we read number or a file.
read(Line,*,err=666,end=666) NA
IreadHere = 1
goto 667
666 continue
Ireadhere = 0
iXU = iZMUnit+1
call molcas_open(iXU,line)
read(iXU,*) NA
667 continue
if (Ireadhere == 1) then
  read(iZMUnit,'(A)',Err=9906,end=9999) Line
else
  read(iXU,'(A)',Err=9906,end=9999) Line
end if
do i=1,NA
  if (Ireadhere == 1) then
    read(iZMUnit,'(A)',Err=9906,end=9999) Line
  else
    read(iXU,'(A)',Err=9906,end=9999) Line
  end if
  call Pick_Words(Line,4,Nwords,Words)
  if (Nwords < 4) goto 9993
  call FoundAtomicNumber(LuWr,Words(1),NAtom,iErr)
  if (iErr /= 0) goto 9998
  if (NAtom >= 0) nAtoms = nAtoms+1
  if (NAtom == -1) nXAtoms = nXAtoms+1
  NAT(nAtoms+nXAtoms) = NAtom
  Symbols(nAtoms+nXAtoms) = trim(Words(1))
  if (NAtom > 0) BasReq(NAtom) = .true.

  call Get_dNumber(Words(2),Dist,iErr)
  Zmat(nAtoms+nXAtoms,1) = Dist
  call Get_dNumber(Words(3),Dist,iErr)
  Zmat(nAtoms+nXAtoms,2) = Dist
  call Get_dNumber(Words(4),Dist,iErr)
  Zmat(nAtoms+nXAtoms,3) = Dist

end do
if (Ireadhere == 0) close(iXU)
! Pre-check Basis Set consistency  BasReq: Atom requiring Basis Set
100 nBasis = 0
do i=1,100
  if (BasReq(i)) nBasis = nBasis+1
end do

goto 9999

9906 iErr = 1
write(LuWr,*) ' [XMatReader]: Unable to read x-matrix file !'
goto 9999

9993 iErr = 1
write(LuWr,*) ' [XMatReader]: X-Matrix incomplete in line'
write(LuWr,*) '               ',Line
goto 9999

9998 iErr = 1
write(LuWr,*) ' [XMatReader]: Error in line'
write(LuWr,*) '               ',Line
goto 9999

9999 return

end subroutine XMatReader
