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

subroutine DumInt( &
#                 define _CALLING_
#                 include "int_interface.fh"
                 )
!***********************************************************************
!                                                                      *
! Object: dummy routine that should never be actually called.          *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Definitions, only: wp, iwp

implicit none
#include "int_interface.fh"

#include "macros.fh"
unused_var(Alpha)
unused_var(nAlpha)
unused_var(Beta)
unused_var(nBeta)
unused_var(Zeta)
unused_var(ZInv)
unused_var(rKappa)
unused_var(P)
unused_var(rFinal)
unused_var(nZeta)
unused_var(nIC)
unused_var(nComp)
unused_var(la)
unused_var(lb)
unused_var(A)
unused_var(RB)
unused_var(nHer)
unused_var(Array)
unused_var(nArr)
unused_var(Ccoor)
unused_var(nOrdOp)
unused_var(lOper)
unused_var(iChO)
unused_var(iStabM)
unused_var(nStabM)
unused_var(PtChrg)
unused_var(nGrid)
unused_var(iAddPot)

call WarningMessage(2,'DumInt should never be called')
call Abend()

end subroutine DumInt
