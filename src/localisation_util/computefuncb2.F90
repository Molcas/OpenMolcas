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

implicit real*8(a-h,o-z)
integer ipLbl(nComp)
logical Debug
#include "WrkSpc.fh"

Functional = 0.0d0
do iComp=1,nComp
  ip0 = ipLbl(iComp)-1
  do i=1,nOrb2Loc
    Functional = Functional+(Work(ip0+nOrb2Loc*(i-1)+i))**2
  end do
end do

if (Debug) then
  write(6,*)
  write(6,*) 'In ComputeFuncB2'
  write(6,*) '----------------'
  write(6,*) 'Functional B2 = ',Functional
  write(6,*) '[Assuming doubly occupied orbitals]'
  do iComp=1,nComp
    ip0 = ipLbl(iComp)-1
    Cmp = 0.0d0
    do iMO=1,nOrb2Loc
      Cmp = Cmp+Work(ip0+nOrb2Loc*(iMO-1)+iMO)
    end do
    Cmp = 2.0d0*Cmp
    write(6,'(A,I5,1X,F15.8)') 'Component, Exp. Val.:',iComp,Cmp
    do j=1,nOrb2Loc-1
      do i=j+1,nOrb2Loc
        kij = ip0+nOrb2Loc*(j-1)+i
        kji = ip0+nOrb2Loc*(i-1)+j
        Tst = Work(kij)-Work(kji)
        if (abs(Tst) > 1.0d-14) then
          write(6,*) 'ComputeFuncB2: broken symmetry!'
          write(6,*) '  Component: ',iComp
          write(6,*) '  i and j  : ',i,j
          write(6,*) '  Dij      : ',Work(kij)
          write(6,*) '  Dji      : ',Work(kji)
          write(6,*) '  Diff.    : ',Tst
          call SysAbendMsg('ComputeFuncB2','Broken symmetry!',' ')
        end if
      end do
    end do
  end do
end if

end subroutine ComputeFuncB2
