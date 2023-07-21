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
! Copyright (C) 2007, Ten-no Research Group                            *
!               2012, Thomas Bondo Pedersen                            *
!***********************************************************************

function SumLeg(K_Lap,W,X,Delta)

implicit real*8(A-H,O-Z)
real*8 W(K_Lap), X(K_Lap)

Dum = 0.0D+00
do I=1,K_Lap
  Dum = Dum+W(I)*FuncLe(X(I),Delta)
end do
SumLeg = Dum

return

end function SumLeg
