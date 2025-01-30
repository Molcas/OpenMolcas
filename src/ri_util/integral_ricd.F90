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

use Int_Options, only: iTOffs
use Definitions, only: wp, iwp, u6

implicit none
#include "int_wrout_interface.fh"
integer(kind=iwp) :: iAO(4), iAOst(4), iBas, iCmp(4), iShell(4), jBas, kBas, kOp(4), lBas

#include "macros.fh"
unused_var(iShell)
unused_var(SOInt(1))
unused_var(nSOint)
unused_var(iSOSym)

iCmp(:) = iSD4(2,:)
iShell(:) = iSD4(11,:)
iAO(:) = iSD4(7,:)
iAOst(:) = iSD4(8,:)
iBas = iSD4(19,1)
jBas = iSD4(19,2)
kBas = iSD4(19,3)
lBas = iSD4(19,4)

if (mSym == 1) then
  ! note that iTOffs is being abused for something else
  kOp(:) = 0
  call PLF_RICD(AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),iAO,iAOst,iBas,jBas,kBas,lBas,kOp,TInt,iTOffs(2),iTOffs(3),iTOffs(1), &
                iTOffs(4),iTOffs(5))
else
  write(u6,*) 'Integral_RICD: fatal error!'
  call Abend()
end if

end subroutine Integral_RICD
