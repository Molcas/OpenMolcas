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
! Copyright (C) Anders Ohrn                                            *
!***********************************************************************

!----------------------------------------------------------------------*
! A function that returns the binomial coefficient. The coefficients   *
! are stored since N and P will not under normal circumstances be      *
! so large.                                                            *
!----------------------------------------------------------------------*
integer function NoverP_Q(N,P)

integer N, P, Bino(22)
data(Bino(i),i=1,21)/1,1,1,1,2,1,1,3,3,1,1,4,6,4,1,1,5,10,10,5,1/
#include "warnings.h"

NoverP_Q = 1
if (N >= 6) then
  write(6,*) 'Must extend NoverP_Q!'
  call Quit(_RC_INTERNAL_ERROR_)
else
  ind = (N+1)*(N+2)/2-(N-P)
  NoverP_Q = Bino(ind)
end if

return

end function NoverP_Q
