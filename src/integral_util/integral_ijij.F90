!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Integral_ijij( &
#                        define _CALLING_
#                        include "int_wrout_interface.fh"
                        )
! calls the proper routines IndSft/PLF
! if IntOrd_jikl == .true. integral order within symblk: jikl
!                    else  integral order within symblk: ijkl

use Definitions, only: wp, iwp

implicit none
#include "int_wrout_interface.fh"
integer(kind=iwp) :: iCmp(4), iMax
integer(kind=iwp), parameter :: inc = 1
integer(kind=iwp), external :: iDAMax_

#include "macros.fh"
unused_var(iSOSym)

iCmp(:) = iSD4(2,:)

if (mSym == 1) then
  ! AOInt(ijkl*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4))
  iMax = iDAMax_(ijkl*iCmp(1)*iCmp(2)*iCmp(3)*iCmp(4),AOInt,inc)
  TInt(1) = AOInt(iMax)
else
  ! SOInt(ijkl*nSOInt)
  iMax = iDAMax_(ijkl*nSOInt,SOInt,inc)
  TInt(1) = SOInt(iMax)
end if

end subroutine Integral_ijij
