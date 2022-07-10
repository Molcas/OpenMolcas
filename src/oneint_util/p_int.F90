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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine P_Int( &
#                define _CALLING_
#                include "int_interface.fh"
                )
!***********************************************************************
!                                                                      *
! Object: to compute the multipole moments integrals with the          *
!         Gauss-Hermite quadrature.                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             November '90                                             *
!             Modified to multipole moments November '90               *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "int_interface.fh"
#include "print.fh"
integer(kind=iwp) :: ia, ib, iIC, iPrint, iRout
character(len=80) :: Label

#include "macros.fh"
unused_var(Alpha)
unused_var(Beta)
unused_var(Zeta)
unused_var(ZInv)
unused_var(rKappa)
unused_var(P)
unused_var(A)
unused_var(RB)
unused_var(nHer)
unused_var(Array)
unused_var(Ccoor)
unused_var(nOrdOp)
unused_var(lOper)
unused_var(iChO)
unused_var(iStabM)
unused_var(PtChrg)
unused_var(iAddPot)

iRout = 122
iPrint = nPrint(iRout)
! Observe that this code does not make any sense in case of symmetry!
rFinal(:,:,:,:) = Zero

if (iPrint >= 99) then
  write(u6,*) ' Result in P_Int'
  do ia=1,nTri_Elem1(la)
    do ib=1,nTri_Elem1(lb)
      do iIC=1,nIC
        write(Label,'(A,I2,A,I2,A,I2,A)') ' rFinal(a=',ia,',b=',ib,',iIC=',iIC,')'
        call RecPrt(Label,' ',rFinal(:,ia,ib,iIC),nAlpha,nBeta)
      end do
    end do
  end do
end if

return

end subroutine P_Int
