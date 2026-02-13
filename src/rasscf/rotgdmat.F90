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

subroutine RotGDMat(R,GD)

use rasscf_global, only: lRoots, NAC
use Constants, only: Zero

implicit none
real*8, dimension(LRoots*(LRoots+1)/2,NAC,NAC) :: GD, GD2
real*8, dimension(lRoots,lRoots) :: R
integer I, J, K, L, p, q, iIJ, iKL, ip, iq
#include "warnings.h"

do p=1,nac
  do q=1,nac
    do I=1,lRoots
      do J=1,I
        iIJ = (I-1)*I/2+J
        GD2(iIJ,p,q) = Zero
        do K=1,lRoots
          do L=1,lRoots
            if (K > L) then
              iKL = (K-1)*K/2+L
              ip = p
              iq = q
            else
              iKL = (L-1)*L/2+K
              ip = q
              iq = p
            end if
            GD2(iIJ,p,q) = GD2(iIJ,p,q)+GD(iKL,ip,iq)*R(I,K)*R(J,L)
          end do
        end do
      end do
    end do
  end do
end do

do p=1,nac
  do q=1,nac
    do I=1,lRoots
      do J=1,I
        iIJ = (I-1)*I/2+J
        GD(iIJ,p,q) = GD2(iIJ,p,q)
      end do
    end do
  end do
end do

end subroutine RotGDMat
