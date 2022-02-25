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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************

real*8 function Vt_lim(rho,drho,ddrho)

implicit real*8(a-h,o-z)
real*8 rho, drho(3), ddrho
#include "real.fh"
parameter(One8=One/Eight)
parameter(One4=One/Four)

rhoinv = One/rho
rhoinv2 = rhoinv**Two
xnorm = drho(1)**2+drho(2)**2+drho(3)**2

Vt_lim = One8*xnorm*rhoinv2-One4*ddrho*rhoinv

end function Vt_lim
