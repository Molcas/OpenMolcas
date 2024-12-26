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

subroutine No_Routine( &
#                     define _CALLING_
#                     include "int_wrout_interface.fh"
                     )

use Definitions, only: wp, iwp

implicit none
#include "int_wrout_interface.fh"

#include "macros.fh"
unused_var(ijkl)
unused_var(AOInt(1))
unused_var(SOInt(1))
unused_var(nSOint)
unused_var(iSOSym)
unused_var(TInt)
unused_var(mSym)
unused_var(nSD)
unused_var(iSD4)

end subroutine No_Routine
