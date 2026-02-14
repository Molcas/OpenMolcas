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
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: R(lRoots,lRoots), GD(lRoots*(lRoots+1)/2,NAC,NAC)
integer(kind=iwp) :: I, iIJ, iKL, ip, iq, J, K, L, p, q
real(kind=wp), allocatable :: GD2(:,:,:)

call mma_allocate(GD2,lRoots*(lRoots+1)/2,NAC,NAC,Label='GD2')

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

call mma_deallocate(GD2)

end subroutine RotGDMat
