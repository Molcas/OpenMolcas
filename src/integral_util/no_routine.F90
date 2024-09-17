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
unused_var(iCmp)
unused_var(iShell)
unused_var(iBas)
unused_var(jBas)
unused_var(kBas)
unused_var(lBas)
unused_var(kOp)
unused_var(Shijij)
unused_var(iAO)
unused_var(iAOst)
unused_var(ijkl)
unused_var(AOInt(1))
unused_var(SOInt(1))
unused_var(nSOint)
unused_var(iSOSym)
unused_var(TInt)
unused_var(mSym)

end subroutine No_Routine
