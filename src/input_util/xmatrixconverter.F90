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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************

subroutine XMatrixConverter(LuRd,LuWr,mxAtom,STDINP,lSTDINP,iglobal,nxbas,xb_label,xb_bas,iErr)
!***********************************************************************
! Author: Valera Veryazov                                              *
!                                                                      *
! This is an adaptation of GG Program ZMatrixConverter                 *
!***********************************************************************

use ZMatConv_Mod, only: BasAva, Base, BasReq, Coords, MaxAtoms, NAT, Num_Elem, Symbols, Zmat
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LuRd, LuWr, mxAtom, iglobal
character(len=180), intent(out) :: STDINP(mxAtom*2)
integer(kind=iwp), intent(out) :: lSTDINP, iErr
integer(kind=iwp), intent(inout) :: nxbas
character(len=*), intent(inout) :: xb_label(*), xb_bas(*)
integer(kind=iwp) :: i, i4, i5, ineedfix, iSTDINP, j, k, nAtoms, NATprev, nBase, nBasis, nXAtoms
real(kind=wp) :: r
logical(kind=iwp) :: IfTest
character(len=180) :: aDebug
character(len=12) :: Angstring
character(len=5) :: ll
character(len=4) :: lll

#ifdef _DEBUGPRINT_
IfTest = .true.
#else
IfTest = .false.
#endif

!  ***  H-Fm (Atomic numbers 1-100)
!  ***  X dummy atoms (NA = 0 )
!  ***  Z ghost atoms (NA =-1 )
!  ***  nAskAtoms.EQ.-1  =>  Seward ZMAT input
!  ***  nAskAtoms.NE.-1  =>  GateWay ZMAT input

! nAtoms : nr. of atoms passed to SEWARD (includes X dummy atoms).
! nXAtoms: nr. of ghost Z atoms (not passed to SEWARD but resumed
!          by OutZMat in SLAPAF).
! nBase  : number of BasisSets found in input.
! Base(i): BasisSet for atom with Atomic Number -i-.
! BasAva(i) & BasReq(i): Logical to check BasisSet-consistency.
! Coords(_,i): X, Y, Z, coordinates (in Angstrom) for atom -i-.

nAtoms = 0
nXAtoms = 0
nBase = 0
lSTDINP = 0
do i=1,Num_Elem
  Base(i) = ' '
  BasAva(i) = .false.
  BasReq(i) = .false.
end do
do i=1,MaxAtoms
  Coords(i,:) = Zero
end do
nBasis = 0
iErr = 0
Angstring = '  / Angstrom'

! Reading input
iErr = 0
call BasisReader(LuWr,nBase,iglobal,nxbas,xb_label,xb_bas,iErr)
if (IfTest) then
  write(LuWr,*)
  write(LuWr,*) '------------------------------------------------'
  write(LuWr,*) 'XMatrixConverter - From BasisReader :'
  write(LuWr,*) '                   nBase=',nBase
  do i=1,Num_Elem
    if (BasAva(i)) write(LuWr,'(I23,3X,A)') i,Base(i)
  end do
  write(LuWr,*)
end if
if (iErr /= 0) goto 9905
iErr = 0
call XMatReader(LuRd,LuWr,nAtoms,nXAtoms,nBasis,-1,nxbas,xb_label,xb_bas,iErr)
if (IfTest) then
  do i=1,nAtoms+nXAtoms
    write(LuWr,'(1X,A,I3,3(F10.6))') Symbols(i),NAT(i),(Zmat(i,j),j=1,3)
  end do
  write(LuWr,*)
end if
if (iErr /= 0) goto 9906

! Some checks
if (nBase == 0) then
  write(LuWr,*) 'ERROR: No basis set specified !'
  goto 9999
end if
if (nAtoms == 0) then
  write(LuWr,*) 'ERROR: No atom coordinates specified !'
  goto 9999
end if
if (nBase < nBasis) then
  write(LuWr,*) 'ERROR: Wrong number of basis sets !'
  write(LuWr,*) '       Available=',nBase,'  Required=',nBasis
  goto 9999
end if
call BasisConsistency(LuWr,iErr)
if (iErr /= 0) then
  write(LuWr,*) 'ERROR: Basis set inconsistency !'
  goto 9999
end if
call Put_iScalar('N ZMAT',0)
do i=1,nAtoms+nXAtoms
  Coords(i,1) = Zmat(i,1)
  Coords(i,2) = Zmat(i,2)
  Coords(i,3) = Zmat(i,3)
end do

if (IfTest) then
  write(LuWr,*)
  write(LuWr,*) '------------------------------------------------'
  write(LuWr,*) 'ZMatrixConverter - XYZCoords (Angstroms) :'
  do i=1,nAtoms+nXAtoms
    write(LuWr,99) i,NAT(i),(Coords(i,j),j=1,3)
  end do
  write(LuWr,*)
end if

! Check for superposed atoms
do i=1,nAtoms+nXAtoms
  if (NAT(i) > 0) then
    do j=i+1,nAtoms+nXAtoms
      if (NAT(j) > 0) then
        r = Zero
        do k=1,3
          r = r+(Coords(i,k)-Coords(j,k))**2
        end do
        if (r < 1.0e-4_wp) goto 9907
      end if
    end do
  end if
end do

! Writing
!2000  Continue
if (NAT(1) == -1) then
  NATprev = -1
else
  NATprev = -9999
end if
iSTDINP = 1
do i=1,nAtoms+nXAtoms
  if (NAT(i) == -1) goto 2100
  if (NAT(i) /= NATprev) then
    if (i /= 1 .and. NATprev /= -1) then
      write(STDINP(iSTDINP),'(A)') 'End of basis'
      iSTDINP = iSTDINP+1
    end if
    write(STDINP(iSTDINP),'(A)') 'Basis set'
    iSTDINP = iSTDINP+1
    if (NAT(i) > 0) then
      write(STDINP(iSTDINP),'(A)') Base(NAT(i))
    else
      write(STDINP(iSTDINP),'(A)') 'X..... / InLine '
      iSTDINP = iSTDINP+1
      write(STDINP(iSTDINP),'(A)') '0.   0'
      iSTDINP = iSTDINP+1
      write(STDINP(iSTDINP),'(A)') '0    0'
    end if
    iSTDINP = iSTDINP+1
    NATprev = NAT(i)
  end if
  ll = Symbols(i)
  write(lll,'(i4)') i
  ineedfix = 1
  do i5=1,5
    if (index('0123456789',ll(i5:i5)) /= 0) ineedfix = 0
  end do
  if (ineedfix == 1) then
    i5 = index(ll,' ')
    do i4=1,4
      if (lll(i4:i4) /= ' ') then
        ll(i5:i5) = lll(i4:i4)
        i5 = i5+1
      end if
    end do
  end if
  write(STDINP(iSTDINP),'(A5,3F24.18,A)') ll,(Coords(i,j),j=1,3),Angstring
  iSTDINP = iSTDINP+1
2100 continue
end do
write(STDINP(iSTDINP),'(A)') 'End of basis'
lSTDINP = iSTDINP
if (IfTest) then
  write(LuWr,*)
  write(LuWr,*) '------------------------------------------------'
  write(LuWr,*) 'XMatrixConverter - The input passed to SEWARD : '
  do i=1,iSTDINP
    aDebug = STDINP(i)
    write(LuWr,*) aDebug
  end do
  write(LuWr,*)
end if
goto 9999

99 format(I3,1X,I3,1X,3(F12.6))

9905 write(LuWr,*) ' ERROR: Wrong input in Bases Set definition !'
goto 9999
9906 write(LuWr,*) ' ERROR: Wrong input in Z-Matrix definition !'
goto 9999
9907 write(LuWr,*) ' ERROR: Superimposed atoms: ',i,j,'  r=',sqrt(r)
goto 9999

9999 continue

end subroutine XMatrixConverter
