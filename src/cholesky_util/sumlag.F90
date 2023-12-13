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

function SumLag(K_Lap,W,X,Delta)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: SumLag
integer(kind=iwp), intent(in) :: K_Lap
real(kind=wp), intent(in) :: W(K_Lap), X(K_Lap), Delta
integer(kind=iwp) :: I
real(kind=wp) :: Dum
real(kind=wp), external :: FuncLa

Dum = Zero
do I=1,K_Lap
  Dum = Dum+W(I)*FuncLa(X(I),Delta)
end do
SumLag = Dum

return

end function SumLag
