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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine ComputeFuncB2(nOrb2Loc,ipLbl,nComp,Functional,Debug)
! Author: T.B. Pedersen
!
! Purpose: compute Boys localisation functional B2.

use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nOrb2Loc, nComp, ipLbl(nComp)
real(kind=wp) :: Functional
logical(kind=iwp) :: Debug
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iComp, iMO, ip0, j, kij, kji
real(kind=wp) :: Cmp, Tst

Functional = Zero
do iComp=1,nComp
  ip0 = ipLbl(iComp)-1
  do i=1,nOrb2Loc
    Functional = Functional+(Work(ip0+nOrb2Loc*(i-1)+i))**2
  end do
end do

if (Debug) then
  write(u6,*)
  write(u6,*) 'In ComputeFuncB2'
  write(u6,*) '----------------'
  write(u6,*) 'Functional B2 = ',Functional
  write(u6,*) '[Assuming doubly occupied orbitals]'
  do iComp=1,nComp
    ip0 = ipLbl(iComp)-1
    Cmp = Zero
    do iMO=1,nOrb2Loc
      Cmp = Cmp+Work(ip0+nOrb2Loc*(iMO-1)+iMO)
    end do
    Cmp = Two*Cmp
    write(u6,'(A,I5,1X,F15.8)') 'Component, Exp. Val.:',iComp,Cmp
    do j=1,nOrb2Loc-1
      do i=j+1,nOrb2Loc
        kij = ip0+nOrb2Loc*(j-1)+i
        kji = ip0+nOrb2Loc*(i-1)+j
        Tst = Work(kij)-Work(kji)
        if (abs(Tst) > 1.0e-14_wp) then
          write(u6,*) 'ComputeFuncB2: broken symmetry!'
          write(u6,*) '  Component: ',iComp
          write(u6,*) '  i and j  : ',i,j
          write(u6,*) '  Dij      : ',Work(kij)
          write(u6,*) '  Dji      : ',Work(kji)
          write(u6,*) '  Diff.    : ',Tst
          call SysAbendMsg('ComputeFuncB2','Broken symmetry!',' ')
        end if
      end do
    end do
  end do
end if

end subroutine ComputeFuncB2
