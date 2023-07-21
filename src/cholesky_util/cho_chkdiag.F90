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

use ChoArr, only: iRS2F

implicit real*8(a-h,o-z)
real*8 Diag(l_D), Stat(7)
logical DoPrint
#include "Molcas.fh"
#include "cholesky.fh"
#include "choorb.fh"
character*(LENIN8) Name(maxbfn)
character*(LENIN) ctmp1, ctmp2
real*8 Err(4)

call Get_cArray('Unique Basis Names',Name,LENIN8*nBasT)

do krs=1,nnBstRT(1)
  ia = iRS2F(1,krs)
  ctmp1 = Name(ia)(1:LENIN)
  ib = iRS2F(2,krs)
  ctmp2 = Name(ib)(1:LENIN)
  if (ctmp1 /= ctmp2) Diag(krs) = 0.0d0
end do

if (DoPrint) call Cho_Head('Analysis of Difference (1-Center only)','=',80,6)
call Statistics(Diag,l_D,Stat,1,2,3,4,5,6,7)
if (DoPrint) call Cho_PrtSt(Diag,l_D,Stat)

! Set Err array.
! --------------

Err(1) = Stat(3)
Err(2) = Stat(4)
Err(3) = Stat(1)
Err(4) = sqrt(dDot_(nnBstRT(1),Diag,1,Diag,1)/dble(nnBstRT(1)))

if (DoPrint) then
  write(6,'(/,1X,A,1P,D15.6)') 'Minimum error   : ',Err(1)
  write(6,'(1X,A,1P,D15.6)') 'Maximum error   : ',Err(2)
  write(6,'(1X,A,1P,D15.6)') 'Average error   : ',Err(3)
  write(6,'(1X,A,1P,D15.6)') 'RMS error       : ',Err(4)
end if

return

end subroutine OneCenter_ChkDiag
