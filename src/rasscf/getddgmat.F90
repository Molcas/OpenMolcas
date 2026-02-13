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

use rasscf_global, only: lRoots, NAC
use Constants, only: Zero

implicit none
real*8, dimension(lRoots,lRoots,lRoots,lRoots) :: DDG
real*8, dimension(NAC,NAC,NAC,NAC) :: Gtuvx
real*8, dimension(LRoots*(LRoots+1)/2,NAC,NAC) :: GDMat
integer iI, iJ, iK, iL, it, iu, iv, ix, iII, iJJ, iKK, iLL
#include "warnings.h"

do iI=1,lRoots
  do iJ=1,lRoots
    if (iJ > iI) then
      iJJ = iI
      iII = iJ
    else
      iII = iI
      iJJ = iJ
    end if
    do iK=1,lRoots
      do iL=1,lRoots
        if (iL > iK) then
          iLL = iK
          iKK = iL
        else
          iLL = iL
          iKK = iK
        end if
        DDG(iI,iJ,iK,iL) = Zero
        do it=1,NAC
          do iu=1,NAC
            do iv=1,NAC
              do ix=1,NAC
                DDG(iI,iJ,iK,iL) = DDG(iI,iJ,iK,iL)+GDMat(iII*(iII-1)/2+iJJ,it,iu)*GDMat(iKK*(iKK-1)/2+iLL,iv,ix)*Gtuvx(it,iu,iv,ix)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do

return

end subroutine GetDDgMat
