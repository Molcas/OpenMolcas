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

implicit none
integer nab, nijx, nScrt, n1, n2, mab
real*8 Win(nab*nijx), Scrt(nScrt), Coeff1((n1+1)*(n1+2)/2,(n1+1)*(n1+2)/2), Coeff2((n2+1)*(n2+2)/2,(n2+1)*(n2+2)/2), Wout(mab*nijx)
logical Tr1, Pr1, Tr2, Pr2
integer l1, k1, l2, k2

l1 = (n1+1)*(n1+2)/2
k1 = l1
if (Pr1) k1 = 2*n1+1
l2 = (n2+1)*(n2+2)/2
k2 = l2
if (Pr2) k2 = 2*n2+1

if (Tr1 .and. Tr2) then

  ! Starting with a,bIJx transforming to bIJx,A

  call DGEMM_('T','N',l2*nijx,k1,l1,1.0d0,Win,l1,Coeff1,l1,0.0d0,Scrt,l2*nijx)

  ! Transform b,IJxA to IJxAB

  call DGEMM_('T','N',nijx*k1,k2,l2,1.0d0,Scrt,l2,Coeff2,l2,0.0d0,Wout,nijx*k1)

else if (Tr2) then

  ! Transpose from ab,IJ,x to b,IJ,x,a

  call DGeTmO(Win,l1,l1,l2*nijx,Scrt,l2*nijx)

  ! Start transforming b,IJ,x,a to IJ,x,aB

  call DGEMM_('T','N',nijx*l1,k2,l2,1.0d0,Scrt,l2,Coeff2,l2,0.0d0,Wout,nijx*l1)
else

  ! Starting with a,bIJx transforming to AbIJx

  call DGEMM_('T','N',k1,l2*nijx,l1,1.0d0,Coeff1,l1,Win,l1,0.0d0,Scrt,k1)

  ! Transpose to IJxAb

  call DGeTmO(Scrt,k1*l2,k1*l2,nijx,Wout,nijx)
end if

return

end subroutine CarSph
