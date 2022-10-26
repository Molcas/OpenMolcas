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

subroutine Integral_RICD( &
#                        define _CALLING_
#                        include "int_wrout_interface.fh"
                        )

use Definitions, only: wp, iwp, u6

implicit none
#include "int_wrout_interface.fh"

#include "macros.fh"
unused_var(iShell)
unused_var(MapOrg)
unused_var(Shijij)
unused_var(IJeqKL)
unused_var(SOInt(1))
unused_var(nSOint)
unused_var(iSOSym)
unused_var(nSkal)

if (mSym == 1) then
  ! note that iTOffs is being abused for something else
  call PLF_RICD(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iAO,iAOst,iBas,jBas,kBas,lBas,kOp,TInt,iTOffs(0,0,1),iTOffs(0,0,2), &
                iTOffs(0,0,0),iTOffs(0,0,3),iTOffs(0,0,4))
else
  write(u6,*) 'Integral_RICD: fatal error!'
  call Abend()
end if

return

end subroutine Integral_RICD
