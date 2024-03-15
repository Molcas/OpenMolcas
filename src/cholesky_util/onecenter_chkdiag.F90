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
! Copyright (C) Thomas Bondo Pedersen                                  *
!               Francesco Aquilante                                    *
!***********************************************************************

subroutine OneCenter_ChkDiag(Diag,l_D,Stat,DoPrint)

use Cholesky, only: iRS2F, nBasT, nnBstRT
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: l_D
real(kind=wp), intent(inout) :: Diag(l_D)
real(kind=wp), intent(out) :: Stat(7)
logical(kind=iwp), intent(in) :: DoPrint
#include "Molcas.fh"
integer(kind=iwp) :: ia, ib, krs
real(kind=wp) :: Err(4)
character(len=LenIn8) :: BName(maxbfn)
character(len=LenIn) :: ctmp1, ctmp2
real(kind=wp), external :: dDot_

call Get_cArray('Unique Basis Names',BName,LENIN8*nBasT)

do krs=1,nnBstRT(1)
  ia = iRS2F(1,krs)
  ctmp1 = BName(ia)(1:LENIN)
  ib = iRS2F(2,krs)
  ctmp2 = BName(ib)(1:LENIN)
  if (ctmp1 /= ctmp2) Diag(krs) = Zero
end do

if (DoPrint) call Cho_Head('Analysis of Difference (1-Center only)','=',80,u6)
call Statistics(Diag,l_D,Stat,1,2,3,4,5,6,7)
if (DoPrint) call Cho_PrtSt(Diag,l_D,Stat)

! Set Err array.
! --------------

Err(1) = Stat(3)
Err(2) = Stat(4)
Err(3) = Stat(1)
Err(4) = sqrt(dDot_(nnBstRT(1),Diag,1,Diag,1)/real(nnBstRT(1),kind=wp))

if (DoPrint) then
  write(u6,'(/,1X,A,ES15.6)') 'Minimum error   : ',Err(1)
  write(u6,'(1X,A,ES15.6)') 'Maximum error   : ',Err(2)
  write(u6,'(1X,A,ES15.6)') 'Average error   : ',Err(3)
  write(u6,'(1X,A,ES15.6)') 'RMS error       : ',Err(4)
end if

return

end subroutine OneCenter_ChkDiag
