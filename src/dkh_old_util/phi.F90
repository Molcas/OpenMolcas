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
! Copyright (C) 1984,1986, Bernd Artur Hess                            *
!***********************************************************************

real*8 function PHI(M,N)

implicit real*8(A-H,O-Z)
#include "crelop.fh"

if ((mod(N,2) == 1) .or. (mod(M,2) == 1)) goto 10
PHI = 2.d0*GA(M+1)*GA(N+1)/GA(M+N+2)
return
10 PHI = 0.d0

return

end function PHI
