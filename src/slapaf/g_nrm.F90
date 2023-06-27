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

subroutine G_Nrm(nInter,GNrm,Iter,Grad,mIntEff)

use Slapaf_Info, only: Gx, Degen

implicit real*8(a-h,o-z)
#include "real.fh"
#include "Molcas.fh"
real*8 GNrm(Iter), Grad(nInter,Iter)

! Compute the norm of the cartesian force vector.
!
! |dE/dx|=Sqrt(dE/dx|u|dE/dx)

Fabs = 0.0d0
do i=1,size(Gx,2)
  do j=1,3
    Fabs = Fabs+Degen(j,i)*Gx(j,i,Iter)**2
  end do
end do
Fabs = sqrt(Fabs)
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
write(6,42) Fabs
42 format(/,' Norm of the force vector',F20.15)
#endif
GNrm(iter) = Fabs

! Write out the internal force vector.

mIntEff = 0
do i=1,nInter
  if (abs(Grad(i,Iter)) > 1.0d-6) mIntEff = mIntEff+1
end do
if (mIntEff == 0) mIntEff = 1
#ifdef _DEBUGPRINT_
write(6,*) ' mIntEff=',mIntEff
#endif

return

end subroutine G_Nrm
