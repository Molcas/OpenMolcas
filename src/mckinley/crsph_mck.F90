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

subroutine CrSph_mck(Win,nijx,nab,Coeff1,n1,Tr1,Pr1,Wout,mab)
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

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nijx, nab, n1, mab
real(kind=wp), intent(in) :: Win(nab*nijx), Coeff1(nTri_Elem1(n1),nTri_Elem1(n1))
real(kind=wp), intent(_OUT_) :: Wout(mab*nijx)
logical(kind=iwp), intent(in) :: Tr1, Pr1
integer(kind=iwp) :: k1, l1

!iRout = 26
!iPrint = nPrint(iRout)
l1 = nTri_Elem1(n1)
k1 = l1
if (Pr1) k1 = 2*n1+1

if (Tr1) then

  ! Starting with a,bIJx transforming to bIJx,A

  call DGEMM_('T','N',nijx,k1,l1,One,Win,l1,Coeff1,l1,Zero,Wout,nijx)

else

  ! Transpose from ab,IJ,x to b,IJ,x,a

  call DGeTmO(Win,l1,l1,nijx,Wout,nijx)

  ! Start transforming b,IJ,x,a to IJ,x,aB

end if

return

end subroutine CrSph_mck
