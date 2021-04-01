!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine ZMatReader(iZMUnit,LuWr,nAtoms,nXAtoms,nBasis,nAskAtoms,iErr)
! nAtoms : Total number of real atoms. Include X dummy atoms: NAT(i)= 0
! nXAtoms: Total number of ghost (Z) atoms:                   NAT(i)=-1
! nBasis : Nummer of atom types requiring Basis Set
!  ***  nAskAtoms == -1  =>  Seward ZMAT input  => Use "End of"
!  ***  nAskAtoms /= -1  =>  GateWay ZMAT input => Use nAskAtoms

use ZMatConv_Mod, only: BasReq, iZmat, NAT, Symbols, Zmat
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iZMUnit, LuWR, nAskAtoms
integer(kind=iwp), intent(out) :: nAtoms, nXAtoms, nBasis, iErr
integer(kind=iwp) :: i, NA, NAtom, NB, NT, Nwords
integer(kind=iwp) :: istatus
real(kind=wp) :: Beta, Dist, Theta
character(len=80) :: Line
character(len=24) :: Words(7)
character(len=3) :: Command

iErr = 0
nAtoms = 0
nXAtoms = 0
nBasis = 0

do
  ! Read Line (or COMMAND)
  if ((nAtoms+nXAtoms) == nAskAtoms) exit
  read(iZMUnit,'(A)',iostat=istatus) Line
  if (istatus < 0) write(LuWr,*) ' [ZMatReader]: Unable to read z-matrix file !'
  if (istatus /= 0) return
  if (Line(1:1) == '*') cycle
  if (Line == ' ') exit
  Command = Line(1:3)
  call UpCase(Command)
  if (Command == 'END') exit
  iErr = 0
  NA = 0
  NB = 0
  NT = 0
  Dist = Zero
  Beta = Zero
  Theta = Zero

  ! Read Symbol           [ Symb ]
  call Pick_Words(Line,7,Nwords,Words)
  if (Nwords < 1) then
    call error(5)
    return
  end if
  call FoundAtomicNumber(LuWr,Words(1),NAtom,iErr)
  if (iErr /= 0) then
    call error(6)
    return
  end if
  if (NAtom >= 0) nAtoms = nAtoms+1
  if (NAtom == -1) nXAtoms = nXAtoms+1
  if (nAtoms+nXAtoms > size(NAT)) then
    call error(7)
    return
  end if
  NAT(nAtoms+nXAtoms) = NAtom
  Symbols(nAtoms+nXAtoms) = trim(Words(1))
  if (NAtom > 0) BasReq(NAtom) = .true.
  if ((nAtoms+nXAtoms) == 1) cycle ! Read Only the First Atom

  ! Read Distance         [ Symb   NA Dist ]
  call Pick_Words(Line,7,Nwords,Words)
  if (Nwords < 3) then
    call error(5)
    return
  end if
  call Get_iNumber(Words(2),NA,iErr)
  if (iErr /= 0) then
    call error(6)
    return
  end if
  if (NA >= (nAtoms+nXAtoms)) then
    call error(1)
    return
  end if
  call Get_dNumber(Words(3),Dist,iErr)
  if (iErr /= 0) then
    call error(6)
    return
  end if
  if (Dist <= Zero) then
    call error(2)
    return
  end if
  iZmat(1,nAtoms+nXAtoms) = NA
  Zmat(1,nAtoms+nXAtoms) = Dist
  if ((nAtoms+nXAtoms) == 2) cycle ! Read Only the Second Atom

  ! Read Planar angle     [ Symb   NA Dist   NB Beta ]
  call Pick_Words(Line,7,Nwords,Words)
  if (Nwords < 5) then
    call error(5)
    return
  end if
  call Get_iNumber(Words(4),NB,iErr)
  if (iErr /= 0) then
    call error(6)
    return
  end if
  if (NB >= (nAtoms+nXAtoms)) then
    call error(1)
    return
  end if
  call Get_dNumber(Words(5),Beta,iErr)
  if (iErr /= 0) then
    call error(6)
    return
  end if
  if ((Beta <= Zero) .or. (Beta >= 180.0_wp)) then
    call error(3)
    return
  end if
  iZmat(2,nAtoms+nXAtoms) = NB
  Zmat(2,nAtoms+nXAtoms) = Beta
  if (NA == NB) then
    call error(4)
    return
  end if
  if ((nAtoms+nXAtoms) == 3) cycle ! Read Only the Second Atom

  ! Read Dihedral angle   [ Symb   NA Dist   NB Beta   NT Theta]
  call Pick_Words(Line,7,Nwords,Words)
  if (Nwords < 7) then
    call error(5)
    return
  end if
  call Get_iNumber(Words(6),NT,iErr)
  if (iErr /= 0) then
    call error(6)
    return
  end if
  if (NT >= (nAtoms+nXAtoms)) then
    call error(1)
    return
  end if
  call Get_dNumber(Words(7),Theta,iErr)
  if (iErr /= 0) then
    call error(6)
    return
  end if
  iZmat(3,nAtoms+nXAtoms) = NT
  Zmat(3,nAtoms+nXAtoms) = Theta
  if ((NA == NB) .or. (NB == NT) .or. (NA == NT)) then
    call error(4)
    return
  end if
end do

! Pre-check Basis Set consistency  BasReq: Atom requiring Basis Set
nBasis = 0
do i=1,size(BasReq)
  if (BasReq(i)) nBasis = nBasis+1
end do

return

contains

subroutine error(code)
  integer(kind=iwp), intent(in) :: code
  iErr = 1
  select case (code)
    case (1)
      write(LuWr,*) ' [ZMatReader]: Wrong index in line'
    case (2)
      write(LuWr,*) ' [ZMatReader]: Wrong distance in line'
    case (3)
      write(LuWr,*) ' [ZMatReader]: Wrong planar angle in line'
    case (4)
      write(LuWr,*) ' [ZMatReader]: Multiple index in line'
    case (5)
      write(LuWr,*) ' [ZMatReader]: Z-Matrix incomplete in line'
    case (6)
      write(LuWr,*) ' [ZMatReader]: Error in line'
    case (7)
      write(LuWr,*) ' [ZMatReader]: Too many atoms'
    case default
  end select
  write(LuWr,*) '               ',Line
end subroutine error

end subroutine ZMatReader
