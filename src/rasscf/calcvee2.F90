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

use rasscf_global, only: lRoots, NAC
use Constants, only: Zero, Half

implicit none
real*8, dimension(LRoots*(LRoots+1)/2,NAC,NAC) :: GD
real*8, dimension(lRoots) :: Vee
real*8, dimension(NAC,NAC,NAC,NAC) :: Gtuvx
integer I, t, u, v, x, III
#include "warnings.h"

do I=1,lRoots
  Vee(I) = Zero
  III = I*(I+1)/2
  do t=1,nac
    do u=1,nac
      do v=1,nac
        do x=1,nac
          Vee(I) = Vee(I)+GD(III,t,u)*GD(III,v,x)*Gtuvx(t,u,v,x)
        end do
      end do
    end do
  end do
  Vee(I) = Vee(I)*Half
end do

end subroutine CalcVee2
