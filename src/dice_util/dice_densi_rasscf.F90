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
! Copyright (C) 2017, Quan Phung                                       *
!***********************************************************************

subroutine DICE_DENSI_RASSCF(jRoot,D,DS,PS,PA,PT)

use rasscf_global, only: NAC, NACPAR, NACPR2
use general_data, only: NACTEL
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: jRoot
real(kind=wp), intent(out) :: D(NACPAR), DS(NACPAR), PS(NACPR2), PA(NACPR2), PT(NAC,NAC,NAC,NAC)
integer(kind=iwp) :: I, IJ_pack, IJKL_pack, J, K, L, LLIM
real(kind=wp) :: D1sum

D(:) = Zero
DS(:) = Zero
PS(:) = Zero
PA(:) = Zero

if (NACTEL <= 1) then
  write(u6,*) 'Dice does not allow 1 electron.'
  return
end if

call dice_load2pdm(NAC,PT,jRoot)
IJ_pack = 1
do J=1,NAC
  do I=1,J
    D1sum = Zero
    do K=1,NAC
      D1sum = D1sum+PT(K,K,I,J)
    end do
    D(IJ_pack) = D1sum/(NACTEL-1)
    IJ_pack = IJ_pack+1
  end do
end do

IJKL_pack = 0
do I=1,NAC
  do J=1,I
    do K=1,I
      LLIM = K
      if (K == I) LLIM = J
      do L=1,LLIM
        IJKL_pack = IJKL_pack+1
        if (K == L) then
          PS(IJKL_pack) = Half*PT(L,K,J,I)
        else
          PS(IJKL_pack) = Half*(PT(L,K,J,I)+PT(K,L,J,I))
          PA(IJKL_pack) = Half*(PT(L,K,J,I)-PT(K,L,J,I))
        end if
      end do
    end do
  end do
end do

return

end subroutine DICE_DENSI_RASSCF
