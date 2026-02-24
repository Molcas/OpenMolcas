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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************

subroutine GetDDgMat(DDg,GDMat,Gtuvx)

use Index_Functions, only: iTri, nTri_Elem
use rasscf_global, only: lRoots, NAC
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: DDG(lRoots,lRoots,lRoots,lRoots)
real(kind=wp), intent(in) :: GDMat(nTri_Elem(lRoots),NAC,NAC), Gtuvx(NAC,NAC,NAC,NAC)
integer(kind=iwp) :: iI, iIJ, iJ, iK, iKL, iL, iv, ix

DDG(:,:,:,:) = Zero
do iL=1,lRoots
  do iK=1,lRoots
    iKL = iTri(iL,iK)
    do iJ=1,lRoots
      do iI=1,lRoots
        iIJ = iTri(iI,iJ)
        do ix=1,NAC
          do iv=1,NAC
            DDG(iI,iJ,iK,iL) = DDG(iI,iJ,iK,iL)+sum(GDMat(iKL,iv,ix)*GDMat(iIJ,:,:)*Gtuvx(:,:,iv,ix))
          end do
        end do
      end do
    end do
  end do
end do

return

end subroutine GetDDgMat
