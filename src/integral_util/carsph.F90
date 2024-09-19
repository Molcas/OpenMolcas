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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine CarSph(Win,nab,nijx,Scrt,nScrt,Coeff1,n1,Tr1,Pr1,Coeff2,n2,Tr2,Pr2,Wout,mab)
!***********************************************************************
!                                                                      *
!  Object: to transform the one electron integrals from cartesian      *
!          basis to spherical basis.                                   *
!                                                                      *
! Called from: OneEl                                                   *
!                                                                      *
! Calling    : RecPrt                                                  *
!              DGEMM_   (ESSL)                                         *
!              DGeTMO   (ESSL)                                         *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             February '90                                             *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nab, nijx, nScrt, n1, n2, mab
real(kind=wp), intent(in) :: Win(nab*nijx), Coeff1(nTri_Elem1(n1),nTri_Elem1(n1)), Coeff2(nTri_Elem1(n2),nTri_Elem1(n2))
real(kind=wp), intent(inout) :: Scrt(nScrt)
real(kind=wp), intent(out) :: Wout(mab*nijx)
logical(kind=iwp) :: Tr1, Pr1, Tr2, Pr2
integer(kind=iwp) :: l1, l2, k1, k2

l1 = nTri_Elem1(n1)
k1 = l1
if (Pr1) k1 = 2*n1+1
l2 = nTri_Elem1(n2)
k2 = l2
if (Pr2) k2 = 2*n2+1

if (Tr1 .and. Tr2) then

  ! Starting with a,bIJx transforming to bIJx,A

  call DGEMM_('T','N',l2*nijx,k1,l1,One,Win,l1,Coeff1,l1,Zero,Scrt,l2*nijx)

  ! Transform b,IJxA to IJxAB

  call DGEMM_('T','N',nijx*k1,k2,l2,One,Scrt,l2,Coeff2,l2,Zero,Wout,nijx*k1)

else if (Tr2) then

  ! Transpose from ab,IJ,x to b,IJ,x,a

  call DGeTmO(Win,l1,l1,l2*nijx,Scrt,l2*nijx)

  ! Start transforming b,IJ,x,a to IJ,x,aB

  call DGEMM_('T','N',nijx*l1,k2,l2,One,Scrt,l2,Coeff2,l2,Zero,Wout,nijx*l1)

else

  ! Starting with a,bIJx transforming to AbIJx

  call DGEMM_('T','N',k1,l2*nijx,l1,One,Coeff1,l1,Win,l1,Zero,Scrt,k1)

  ! Transpose to IJxAb

  call DGeTmO(Scrt,k1*l2,k1*l2,nijx,Wout,nijx)

end if

return

end subroutine CarSph
