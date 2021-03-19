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

use ZMatConv_Mod, only: BasReq, NAT, Symbols, Zmat
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iZMUnit, LuWr, nAskAtoms
integer(kind=iwp), intent(out) :: nAtoms, nXAtoms, nBasis, nxbas, iErr
character(len=*), intent(inout) :: xb_label(*), xb_bas(*)
integer(kind=iwp) :: i, istatus, iXU, NA, NAtom, Nwords
logical(kind=iwp) :: IreadHere, Skip
real(kind=wp) :: Dist
character(len=80) :: Line
character(len=24) :: Words(7)
character(len=3) :: Command
integer(kind=iwp), external :: isFreeUnit

xb_label(1) = ' '
xb_bas(1) = ' '
nxbas = 1

iErr = 0
nAtoms = 0
nXAtoms = 0
nBasis = 0

! Read Line (or COMMAND)
Skip = .false.
do
  if ((nAtoms+nXAtoms) == nAskAtoms) then
    Skip = .true.
    exit
  end if
  read(iZMUnit,'(A)',iostat=istatus) Line
  if (istatus > 0) call error()
  if (istatus /= 0) return
  if (Line(1:1) /= '*') exit
end do
if (Line == ' ') Skip = .true.
Command = Line(1:3)
call UpCase(Command)
if (Command == 'END') Skip = .true.
if (.not. Skip) then
  NA = 0
  Dist = Zero

  ! Here we read number or a file.
  read(Line,*,iostat=istatus) NA
  if (istatus == 0) then
    IreadHere = .true.
  else
    IreadHere = .false.
    iXU = isFreeUnit(iZMUnit+1)
    call molcas_open(iXU,line)
    read(iXU,*) NA
  end if
  if (IreadHere) then
    read(iZMUnit,'(A)',iostat=istatus) Line
  else
    read(iXU,'(A)',iostat=istatus) Line
  end if
  if (istatus > 0) call error()
  if (istatus /= 0) return
  do i=1,NA
    if (IreadHere) then
      read(iZMUnit,'(A)',iostat=istatus) Line
    else
      read(iXU,'(A)',iostat=istatus) Line
    end if
    if (istatus > 0) call error()
    if (istatus /= 0) return
    call Pick_Words(Line,4,Nwords,Words)
    if (Nwords < 4) then
      iErr = 1
      write(LuWr,*) ' [XMatReader]: X-Matrix incomplete in line'
      write(LuWr,*) '               ',Line
      return
    end if
    call FoundAtomicNumber(LuWr,Words(1),NAtom,iErr)
    if (iErr /= 0) then
      iErr = 1
      write(LuWr,*) ' [XMatReader]: Error in line'
      write(LuWr,*) '               ',Line
      return
    end if
    if (NAtom >= 0) nAtoms = nAtoms+1
    if (NAtom == -1) nXAtoms = nXAtoms+1
    if (nAtoms+nXAtoms > size(NAT)) then
      iErr = 1
      write(LuWr,*) ' [XMatReader]: Too many atoms'
      write(LuWr,*) '               ',Line
      return
    end if
    NAT(nAtoms+nXAtoms) = NAtom
    Symbols(nAtoms+nXAtoms) = trim(Words(1))
    if (NAtom > 0) BasReq(NAtom) = .true.

    call Get_dNumber(Words(2),Dist,iErr)
    Zmat(1,nAtoms+nXAtoms) = Dist
    call Get_dNumber(Words(3),Dist,iErr)
    Zmat(2,nAtoms+nXAtoms) = Dist
    call Get_dNumber(Words(4),Dist,iErr)
    Zmat(3,nAtoms+nXAtoms) = Dist

  end do
  if (.not. IreadHere) close(iXU)
end if
! Pre-check Basis Set consistency  BasReq: Atom requiring Basis Set
nBasis = 0
do i=1,size(BasReq)
  if (BasReq(i)) nBasis = nBasis+1
end do

return

contains

subroutine error()
  iErr = 1
  write(LuWr,*) ' [XMatReader]: Unable to read x-matrix file !'
end subroutine error

end subroutine XMatReader
