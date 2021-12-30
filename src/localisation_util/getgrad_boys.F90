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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine GetGrad_Boys(nOrb2Loc,ipLbl,nComp,Rmat,GradNorm,Debug)
! Thomas Bondo Pedersen, December 2005.
!
! Purpose: compute R-matrix and gradient norm for Boys functional.

implicit real*8(a-h,o-z)
integer ipLbl(nComp)
real*8 Rmat(nOrb2Loc,nOrb2Loc)
logical Debug
#include "WrkSpc.fh"

call FZero(Rmat,nOrb2Loc**2)
do iComp=1,nComp
  ip0 = ipLbl(iComp)-1
  do j=1,nOrb2Loc
    Rjj = Work(ip0+nOrb2Loc*(j-1)+j)
    do i=1,nOrb2Loc
      Rmat(i,j) = Rmat(i,j)+Work(ip0+nOrb2Loc*(j-1)+i)*Rjj
    end do
  end do
end do

GradNorm = 0.0d0
do i=1,nOrb2Loc-1
  do j=i+1,nOrb2Loc
    GradNorm = GradNorm+(Rmat(i,j)-Rmat(j,i))**2
  end do
end do
GradNorm = 4.0d0*sqrt(GradNorm)

if (Debug) then
  Fun = 0.0d0
  do i=1,nOrb2Loc
    Fun = Fun+Rmat(i,i)
  end do
  write(6,*) 'GetGrad_Boys: functional = Tr(R) = ',Fun
end if

end subroutine GetGrad_Boys
