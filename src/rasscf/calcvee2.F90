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

subroutine CalcVee2(Vee,GD,Gtuvx)

use Index_Functions, only: nTri_Elem
use rasscf_global, only: lRoots, NAC
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Vee(lRoots), GD(nTri_Elem(lRoots),NAC,NAC), Gtuvx(NAC,NAC,NAC,NAC)
integer(kind=iwp) :: I, III, v, x

Vee(:) = Zero
do I=1,lRoots
  III = nTri_Elem(I)
  do x=1,nac
    do v=1,nac
      Vee(I) = Vee(I)+GD(III,v,x)*sum(GD(III,:,:)*Gtuvx(:,:,v,x))
    end do
  end do
end do
Vee(:) = Half*Vee(:)

end subroutine CalcVee2
