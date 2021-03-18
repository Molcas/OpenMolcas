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
! Copyright (C) 2006, Giovanni Ghigo                                   *
!***********************************************************************

subroutine ZMatrixConverter_GW(LuRd,LuWr,LuOut,nAskAtoms,iErr)
!***********************************************************************
! Author: Giovanni Ghigo                                               *
!         Torino (Italy)  October-November 2006                        *
!                                                                      *
! This is an adaptation of Subroutine ZMatrixConverter for GateWay     *
!***********************************************************************

implicit real*8(a-h,o-z)
implicit integer(i-n)
#include "constants.fh"
#include "g_zmatconv.fh"

nAtoms = 0
do i=1,Num_Elem
  Base(i) = ' '
  BasAva(i) = .false.
  BasReq(i) = .false.
end do
do i=1,MaxAtoms
  Coords(i,1) = 0.0d0
  Coords(i,2) = 0.0d0
  Coords(i,3) = 0.0d0
end do
nBasis = 0
iErr = 0

! Reading input
iErr = 0
call ZMatReader(LuRd,LuWr,nAtoms,nXAtoms,nBasis,nAskAtoms,iErr)
if (iErr /= 0) goto 9906
write(LuOut,*) nAtoms
write(LuOut,'(A)') 'Angstrom'

! Some checks
if (nAtoms == 0) then
  write(LuWr,*) 'ERROR: No atom coordinates specified !'
  goto 9999
end if

call Put_iScalar('N ZMAT',nAtoms+nXAtoms)
call Put_cArray('Symbol ZMAT',Symbols(1),(nAtoms+nXAtoms)*5)
call Put_iArray('Index ZMAT',iZmat,MaxAtoms*3)
call Put_iArray('NAT ZMAT',NAT,nAtoms+nXAtoms)

! Calculate coordinates
torad = CONST_PI_/180.0d0
! Atom #1
if (nAtoms+nXAtoms == 1) goto 2000
! Atom #2
Coords(2,3) = Zmat(2,1)  ! Z(2)=R
if (nAtoms+nXAtoms == 2) goto 2000
! Atom #3
if (iZmat(3,1) == 1) then
  Coords(3,1) = Zmat(3,1)*sin(Zmat(3,2)*torad) ! X(2)=R sin(A)
  Coords(3,3) = Zmat(3,1)*cos(Zmat(3,2)*torad) ! Z(3)=R cos(A)
else
  Coords(3,1) = Zmat(3,1)*sin(Zmat(3,2)*torad)
  Coords(3,3) = Coords(2,3)-Zmat(3,1)*cos(Zmat(3,2)*torad)
end if
if (nAtoms+nXAtoms == 3) goto 2000
! Atom #4 ->
do iAtom=4,nAtoms+nXAtoms
  call ZMatConv(LuWr,iAtom,iErr)
end do
if (iErr /= 0) goto 9999

! Check for superposed atoms
2000 do i=1,nAtoms+nXAtoms
  if (NAT(i) > 0) then
    do j=i+1,nAtoms+nXAtoms
      if (NAT(j) > 0) then
        r = 0.0d0
        do k=1,3
          r = r+(Coords(i,k)-Coords(j,k))**2
        end do
        if (r < 0.0001d0) goto 9907
      end if
    end do
  end if
end do

! Writing

do i=1,nAtoms+nXAtoms
  if (NAT(i) > 0) write(LuOut,999) Symbols(i),(Coords(i,k),k=1,3)
end do
999 format(A5,1X,3(F12.6))
goto 9999

9906 write(LuWr,*) ' ERROR: Wrong input in Z-Matrix definition !'
goto 9999
9907 write(LuWr,*) ' ERROR: Superimposed atoms: ',i,j,'  r=',sqrt(r)
goto 9999

9999 return

end subroutine ZMatrixConverter_GW
