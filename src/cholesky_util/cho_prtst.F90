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

subroutine Cho_PrtSt(Vec,lVec,Stat)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: lVec
real(kind=wp), intent(in) :: Vec(lVec), Stat(7)
real(kind=wp), external :: dDot_

write(u6,'(/,1X,A,I15)') 'No. of elements: ',lVec
write(u6,'(1X,A,1P,D15.6)') 'Frobenius norm : ',sqrt(dDot_(lVec,Vec,1,Vec,1))
write(u6,'(1X,A,1P,D15.6)') 'Minimum value  : ',Stat(3)
write(u6,'(1X,A,1P,D15.6)') 'Maximum value  : ',Stat(4)
write(u6,'(1X,A,1P,D15.6)') 'Mean value     : ',Stat(1)
write(u6,'(1X,A,1P,D15.6)') 'Mean abs. value: ',Stat(2)
write(u6,'(1X,A,1P,D15.6)') 'Max. abs. value: ',Stat(5)
write(u6,'(1X,A,1P,D15.6)') 'Biased variance: ',Stat(6)
write(u6,'(1X,A,1P,D15.6,A)') 'Standard dev.  : ',Stat(7),' (unbiased variance)'

end subroutine Cho_PrtSt
