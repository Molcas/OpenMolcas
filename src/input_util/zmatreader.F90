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

implicit integer(i-n)
implicit real*8(a-h,o-z)
#include "g_zmatconv.fh"
character*80 Line, Blank
character*3 Command
character*24 Words(7)
! nAtoms : Total number of real atoms. Include X dummy atoms: NAT(i)= 0
! nXAtoms: Total number of ghost (Z) atoms:                   NAT(i)=-1
! nBasis : Nummer of atom types requiring Basis Set
!  ***  nAskAtoms.EQ.-1  =>  Seward ZMAT input  => Use "End of"
!  ***  nAskAtoms.NE.-1  =>  GateWay ZMAT input => Use nAskAtoms

iErr = 0
Blank = ' '
nAtoms = 0
nXAtoms = 0
nBasis = 0
do i=1,100 ! MaxNat
  BasReq(i) = .false.
end do
do i=1,MaxAtoms
  Symbols(i) = '     '
  NAT(i) = 0
  iZmat(i,1) = 0
  iZmat(i,2) = 0
  iZmat(i,3) = 0
  Zmat(i,1) = 0.0d0
  Zmat(i,2) = 0.0d0
  Zmat(i,3) = 0.0d0
end do

! Read Line (or COMMAND)
10 if ((nAtoms+nXAtoms) == nAskAtoms) goto 100
read(iZMUnit,'(A)',Err=9906,end=9999) Line
if (Line(1:1) == '*') goto 10
if (Line == Blank) goto 100
Command = Line(1:3)
call UpCase(Command)
if (Command == 'END') goto 100
iErr = 0
NA = 0
NB = 0
NT = 0
Dist = 0.0d0
Beta = 0.0d0
Theta = 0.0d0

! Read Symbol           [ Symb ]
call Pick_Words(Line,7,Nwords,Words)
if (Nwords < 1) goto 9993
call FoundAtomicNumber(LuWr,Words(1),NAtom,iErr)
if (iErr /= 0) goto 9998
if (NAtom >= 0) nAtoms = nAtoms+1
if (NAtom == -1) nXAtoms = nXAtoms+1
NAT(nAtoms+nXAtoms) = NAtom
Symbols(nAtoms+nXAtoms) = trim(Words(1))
if (NAtom > 0) BasReq(NAtom) = .true.
if ((nAtoms+nXAtoms) == 1) goto 10 ! Raed Only the First Atom

! Read Distance         [ Symb   NA Dist ]
call Pick_Words(Line,7,Nwords,Words)
if (Nwords < 3) goto 9993
call Get_iNumber(Words(2),NA,iErr)
if (iErr /= 0) goto 9998
if (NA >= (nAtoms+nXAtoms)) goto 9997
call Get_dNumber(Words(3),Dist,iErr)
if (iErr /= 0) goto 9998
if (Dist <= 0.0d0) goto 9996
iZmat(nAtoms+nXAtoms,1) = NA
Zmat(nAtoms+nXAtoms,1) = Dist
if ((nAtoms+nXAtoms) == 2) goto 10 ! Raed Only the Second Atom

! Read Planar angle     [ Symb   NA Dist   NB Beta ]
call Pick_Words(Line,7,Nwords,Words)
if (Nwords < 5) goto 9993
call Get_iNumber(Words(4),NB,iErr)
if (iErr /= 0) goto 9998
if (NB >= (nAtoms+nXAtoms)) goto 9997
call Get_dNumber(Words(5),Beta,iErr)
if (iErr /= 0) goto 9998
if (Beta <= 0.0d0 .or. Beta >= 180.0d0) goto 9995
iZmat(nAtoms+nXAtoms,2) = NB
Zmat(nAtoms+nXAtoms,2) = Beta
if (NA == NB) goto 9994
if ((nAtoms+nXAtoms) == 3) goto 10 ! Raed Only the Second Atom

! Read Dihedral angle   [ Symb   NA Dist   NB Beta   NT Theta]
call Pick_Words(Line,7,Nwords,Words)
if (Nwords < 7) goto 9993
call Get_iNumber(Words(6),NT,iErr)
if (iErr /= 0) goto 9998
if (NT >= (nAtoms+nXAtoms)) goto 9997
call Get_dNumber(Words(7),Theta,iErr)
if (iErr /= 0) goto 9998
iZmat(nAtoms+nXAtoms,3) = NT
Zmat(nAtoms+nXAtoms,3) = Theta
if (NA == NB .or. NB == NT .or. NA == NT) goto 9994
goto 10

! Pre-check Basis Set consistency  BasReq: Atom requiring Basis Set
100 nBasis = 0
do i=1,100
  if (BasReq(i)) nBasis = nBasis+1
end do

goto 9999

9906 iErr = 1
write(LuWr,*) ' [ZMatReader]: Unable to read z-matrix file !'
goto 9999

9997 iErr = 1
write(LuWr,*) ' [ZMatReader]: Wrong index in line'
write(LuWr,*) '               ',Line
goto 9999

9996 iErr = 1
write(LuWr,*) ' [ZMatReader]: Wrong distance in line'
write(LuWr,*) '               ',Line
goto 9999

9995 iErr = 1
write(LuWr,*) ' [ZMatReader]: Wrong planar angle in line'
write(LuWr,*) '               ',Line
goto 9999

9994 iErr = 1
write(LuWr,*) ' [ZMatReader]: Multiple index in line'
write(LuWr,*) '               ',Line
goto 9999

9993 iErr = 1
write(LuWr,*) ' [ZMatReader]: Z-Matrix incomplete in line'
write(LuWr,*) '               ',Line
goto 9999

9998 iErr = 1
write(LuWr,*) ' [ZMatReader]: Error in line'
write(LuWr,*) '               ',Line
goto 9999

9999 return

end subroutine ZMatReader
