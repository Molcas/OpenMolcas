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

subroutine InitT(K_Lap,T,RIni,RNew)
!-----------------------------------------------------------------------
! Function : Set the initial T()
!-----------------------------------------------------------------------

use ReMez_mod, only: IW
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: K_Lap
real(kind=wp), intent(out) :: T(2*K_Lap)
real(kind=wp), intent(in) :: RIni, RNew
integer(kind=iwp) :: I
real(kind=wp) :: Val
logical(kind=iwp) :: Dbg

Dbg = .false.
Val = (RNew-One)/(RIni-One)

if (Dbg) write(IW,'(A)') 'T(I)'
do I=1,2*K_Lap
  T(I) = One+(T(I)-One)*Val
  if (Dbg) write(IW,'(I3,2X,F20.14)') I,T(I)
end do

return

end subroutine InitT
