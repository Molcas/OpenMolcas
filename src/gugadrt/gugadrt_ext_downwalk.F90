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

subroutine gugadrt_ext_downwalk()

use gugadrt_global, only: iseg_downwei, ng_sm, nlsm_ext, norb_ext, nu_ae
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: im, imi, imij, imj, iwmij(8)

nu_ae(1) = 1
do im=1,ng_sm
  nu_ae(1+im) = 1+im
  nu_ae(9+im) = 9+im
  nu_ae(17+im) = 17+im
end do

iwmij = 0
iseg_downwei(nu_ae(1)) = 1
do imi=1,ng_sm
  iseg_downwei(nu_ae(1+imi)) = nlsm_ext(imi)
  do imj=imi,ng_sm
    imij = Mul(imi,imj)
    if (imij /= 1) then
      iwmij(imij) = iwmij(imij)+nlsm_ext(imi)*nlsm_ext(imj)
      cycle
    end if
    iwmij(1) = iwmij(1)+nlsm_ext(imi)*(nlsm_ext(imi)-1)/2
  end do
end do
do im=1,ng_sm
  iseg_downwei(nu_ae(9+im)) = iwmij(im)
  iseg_downwei(nu_ae(17+im)) = iwmij(im)
end do
iseg_downwei(nu_ae(18)) = iseg_downwei(nu_ae(18))+norb_ext

return

end subroutine gugadrt_ext_downwalk
