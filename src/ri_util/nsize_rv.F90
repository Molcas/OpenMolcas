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

function nSize_Rv(kS,lS,nShBf,nShell,nIrrep,iOff,nVec)
!***********************************************************************
!                                                                      *
!     Compute the size of Rv(nu,mu,K) and the offsets to the           *
!     different symmetry blocks.                                       *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nSize_Rv
integer(kind=iwp), intent(in) :: kS, lS, nShell, nIrrep, nShBf(0:nIrrep-1,nShell), nVec(0:nIrrep-1)
integer(kind=iwp), intent(out) :: iOff(0:nIrrep-1)
integer(kind=iwp) :: kIrrep, klIrrep, lIrrep, nJ, nK, nKL, nL

nSize_Rv = 0
iOff(:) = 0

if (nIrrep == 1) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (kS /= lS) then
    nK = nShBf(0,kS)
    nL = nShBf(0,lS)
    nKL = nK*nL
  else
    nK = nShBf(0,kS)
    nKL = nTri_Elem(nK)
  end if

  nJ = nVec(0)
  nSize_Rv = nJ*nkl
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do klIrrep=0,nIrrep-1
    iOff(klIrrep) = nSize_Rv

    nKL = 0
    if (kS /= lS) then
      do kIrrep=0,nIrrep-1
        nK = nShBf(kIrrep,kS)
        lIrrep = Mul(klIrrep+1,kIrrep+1)-1
        nL = nShBf(lIrrep,lS)
        nKL = nKL+nK*nL
      end do
    else
      do kIrrep=0,nIrrep-1
        nK = nShBf(kIrrep,kS)
        lIrrep = Mul(klIrrep+1,kIrrep+1)-1
        nL = nShBf(lIrrep,lS)

        if (kIrrep > lIrrep) then
          nKL = nKL+nK*nL
        else if (kIrrep == lIrrep) then
          nKL = nKL+nTri_Elem(nK)
        end if

      end do
    end if

    nJ = nVec(klIrrep)
    nSize_Rv = nSize_Rv+nJ*nKL

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end function nSize_Rv
