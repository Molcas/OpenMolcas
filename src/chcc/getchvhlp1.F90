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

subroutine GetChVHlp1(cGrp,deGrp,yes,NL2,L2Status)
! this routine does:
! check, if Cholesky vector block cGrp,deGrp is actually situated
! in the memory as one of L2x
!
! descrition of parameters:
! cGrp, deGrp - required Group of c, delta in L2(m,c',de') (I)
! yes    -  0 - not loaded  (O)
!           x - loaded, value of L2x (x=1-4)
! NL2         - Number of L2 arrays really reserved in memory (I)
! L2Status    - L2 status matrix (I)
!               L2Status(i,1) - cGrp
!               L2Status(i,2) - deGrp
!               L2Status(i,3) - Pos

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: cGrp, deGrp, NL2, L2Status(4,3)
integer(kind=iwp), intent(out) :: yes
integer(kind=iwp) :: i

yes = 0
do i=1,NL2
  if ((cGrp == L2Status(i,1)) .and. (deGrp == L2Status(i,2))) yes = i
end do

return

end subroutine GetChVHlp1
