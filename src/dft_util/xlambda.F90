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
! Copyright (C) 2010,2012,2017, Francesco Aquilante                    *
!               2015,2017, Alexander Zech                              *
!***********************************************************************

function Xlambda(omega,sigma)

use Constants, only: One
use Definitions, only: wp

implicit none
real(kind=wp) :: Xlambda
real(kind=wp), intent(in) :: omega, sigma

if (sigma*omega > 42.0_wp) then
  Xlambda = One
else
  Xlambda = One-exp(-sigma*omega)
end if

end function Xlambda
