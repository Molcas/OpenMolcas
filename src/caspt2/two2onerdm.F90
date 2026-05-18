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
! Copyright (C) 2014, Naoki Nakatani                                   *
!***********************************************************************

subroutine TWO2ONERDM(NA,NE,G2,G1)
! Compute 1-RDM from 2-RDM
! Written by N. Nakatani, Oct. 2014

use definitions, only: iwp, wp
use constants, only: Zero

implicit none
integer(kind=iwp), intent(in) :: NA, NE
real(kind=wp), intent(in) :: G2(NA,NA,NA,NA)
real(kind=wp), intent(out) :: G1(NA,NA)
integer(kind=iwp) I, J, K
real(kind=wp) G1TMP

do I=1,NA
  do J=1,NA
    G1TMP = Zero
    do K=1,NA
      G1TMP = G1TMP+G2(K,K,J,I)
    end do
    G1(J,I) = G1TMP/dble(NE-1)
  end do
end do

end subroutine TWO2ONERDM
