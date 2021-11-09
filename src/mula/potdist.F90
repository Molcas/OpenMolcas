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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

subroutine PotDist(F,V,Lambda,PED,NumInt,nOsc)
!  Purpose:
!    To give the fractional contributions of the F-matrix to the
!    potential energy.
!
!  Input:
!    F        : Real*8 two dimensional array -  contains
!               the force constants expressed in internal
!    V        : Real*8 two dimensional array  - contains
!               the eigenvectors of F*G as columns.
!    Lambda   : Real*8 array - contains the eigenvalues
!               of F*G.
!
!  Output:
!    PED      : Real*8 three dimensional array - Potential
!               Energy Distribution for each mode.
!  Uses:
!    Linalg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use Constants, only: Zero
use Definitions, only: wp

!use LinAlg
implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 F(NumInt,NumInt)
real*8 V(NumInt,NumInt)
real*8 Denominator
real*8 Lambda(NumInt)
real*8 PED(nOsc,nOsc,nOsc)

! Initialize.

! Calculate Potential Energy Distribution for each mode.
!PED = Zero
call dcopy_(nOsc*nOsc*nOsc,[Zero],0,PED,1)
do i=1,NumInt
  Denominator = max(1.0e-10_wp,Lambda(i))
  do k=1,NumInt
    do l=1,NumInt
      if (k == l) then
        PED(k,k,i) = (V(k,i)**2*F(k,k))/Denominator
      else
        PED(k,l,i) = (2*V(k,i)*V(l,i)*F(k,l))/Denominator
      end if
    end do
  end do
end do

end subroutine PotDist
