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

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "oneswi.fh"
#include "print.fh"
#include "int_interface.fh"
! Local variables
character*80 Label
! Statement function for Cartesian index
nElem(i) = (i+1)*(i+2)/2

iRout = 122
iPrint = nPrint(iRout)
! Observe that this code does not make any sense in case of symmetry!
call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,final,1)

if (iPrint >= 99) then
  write(6,*) ' Result in P_Int'
  do ia=1,nElem(la)
    do ib=1,nElem(lb)
      do iIC=1,nIC
        write(Label,'(A,I2,A,I2,A,I2,A)') ' Final(a=',ia,',b=',ib,',iIC=',iIC,')'
        call RecPrt(Label,' ',final(1,ia,ib,iIC),nAlpha,nBeta)
      end do
    end do
  end do
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(Alpha)
  call Unused_real_array(Beta)
  call Unused_real_array(Zeta)
  call Unused_real_array(ZInv)
  call Unused_real_array(rKappa)
  call Unused_real_array(P)
  call Unused_real_array(A)
  call Unused_real_array(RB)
  call Unused_integer(nHer)
  call Unused_real_array(Array)
  call Unused_real_array(Ccoor)
  call Unused_integer(nOrdOp)
  call Unused_integer_array(lOper)
  call Unused_integer_array(iChO)
  call Unused_integer_array(iStabM)
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine P_Int
