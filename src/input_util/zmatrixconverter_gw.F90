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

#include "compiler_features.h"
#ifdef _HAVE_EXTRA_

subroutine ZMatrixConverter_GW(LuRd,LuWr,LuOut,nAskAtoms,iErr)
!***********************************************************************
! Author: Giovanni Ghigo                                               *
!         Torino (Italy)  October-November 2006                        *
!                                                                      *
! This is an adaptation of Subroutine ZMatrixConverter for GateWay     *
!***********************************************************************

use ZMatConv_Mod, only: BasReq, Coords, iZmat, NAT, Symbols, Zmat
use isotopes, only: MaxAtomNum
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, deg2rad
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LuRd, LuWr, LuOut, nAskAtoms
integer(kind=iwp), intent(out) :: iErr
integer(kind=iwp) :: i, iAtom, j, k, nAtoms, nBasis, nXAtoms
real(kind=wp) :: r

nAtoms = 0
nBasis = 0
iErr = 0

call mma_allocate(BasReq,MaxAtomNum,label='BasReq')
call mma_allocate(NAT,nAskAtoms,label='NAT')
call mma_allocate(Symbols,nAskAtoms,label='Symbols')
call mma_allocate(iZmat,3,nAskAtoms,label='iZmat')
call mma_allocate(Zmat,3,nAskAtoms,label='Zmat')
BasReq(:) = .false.
NAT(:) = 0
Symbols(:) = ''
iZmat(:,:) = 0
Zmat(:,:) = Zero

! Reading input
call ZMatReader(LuRd,LuWr,nAtoms,nXAtoms,nBasis,nAskAtoms,iErr)
if (iErr /= 0) then
  write(LuWr,*) ' ERROR: Wrong input in Z-Matrix definition !'
  return
end if
write(LuOut,*) nAtoms
write(LuOut,'(A)') 'Angstrom'

! Some checks
if (nAtoms == 0) then
  iErr = 1
  write(LuWr,*) 'ERROR: No atom coordinates specified !'
  return
end if

call Put_iScalar('N ZMAT',nAtoms+nXAtoms)
call Put_cArray('Symbol ZMAT',Symbols(1),(nAtoms+nXAtoms)*len(Symbols))
call Put_iArray('Index ZMAT',iZmat,(nAtoms+nXAtoms)*3)
call Put_iArray('NAT ZMAT',NAT,nAtoms+nXAtoms)

call mma_allocate(Coords,3,nAtoms+nXAtoms,label='Coords')
Coords(:,:) = Zero

! Calculate coordinates
! Atom #1
if (nAtoms+nXAtoms > 1) then
  ! Atom #2
  Coords(3,2) = Zmat(1,2)  ! Z(2)=R
end if
if (nAtoms+nXAtoms > 2) then
  ! Atom #3
  if (iZmat(1,3) == 1) then
    Coords(1,3) = Zmat(1,3)*sin(Zmat(2,3)*deg2rad) ! X(2)=R sin(A)
    Coords(3,3) = Zmat(1,3)*cos(Zmat(2,3)*deg2rad) ! Z(3)=R cos(A)
  else
    Coords(1,3) = Zmat(1,3)*sin(Zmat(2,3)*deg2rad)
    Coords(3,3) = Coords(3,2)-Zmat(1,3)*cos(Zmat(2,3)*deg2rad)
  end if
end if
if (nAtoms+nXAtoms > 3) then
  ! Atom #4 ->
  do iAtom=4,nAtoms+nXAtoms
    call ZMatConv(LuWr,iAtom,iErr)
  end do
  if (iErr /= 0) return
end if

! Check for superposed atoms
do i=1,nAtoms+nXAtoms
  if (NAT(i) > 0) then
    do j=i+1,nAtoms+nXAtoms
      if (NAT(j) > 0) then
        r = Zero
        do k=1,3
          r = r+(Coords(k,i)-Coords(k,j))**2
        end do
        if (r < 1.0e-4_wp) then
          iErr = 1
          write(LuWr,*) ' ERROR: Superimposed atoms: ',i,j,'  r=',sqrt(r)
          return
        end if
      end if
    end do
  end if
end do

! Writing

do i=1,nAtoms+nXAtoms
  if (NAT(i) > 0) write(LuOut,999) Symbols(i),(Coords(k,i),k=1,3)
end do

call mma_deallocate(BasReq)
call mma_deallocate(NAT)
call mma_deallocate(Symbols)
call mma_deallocate(iZmat)
call mma_deallocate(Zmat)
call mma_deallocate(Coords)

return

999 format(A5,1X,3(F12.6))

end subroutine ZMatrixConverter_GW

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(ZMatrixConverter_GW)

#endif
