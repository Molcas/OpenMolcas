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

function nSize_3C(kS,lS,nShBf,nShell,nIrrep,iOff,nBas_Aux)
!***********************************************************************
!                                                                      *
!     Compute the size of ({nu,mu}|K) and the offsets to the           *
!     different symmetry blocks.                                       *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: Mul
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nSize_3C
integer(kind=iwp), intent(in) :: kS, lS, nShell, nIrrep, nShBf(0:nIrrep-1,nShell), nBas_Aux(0:nIrrep-1)
integer(kind=iwp), intent(out) :: iOff(3,0:nIrrep-1)
integer(kind=iwp) :: kIrrep, klIrrep, lIrrep, nJ, nK, nKL, nL

nSize_3C = 0
iOff(:,:) = 0

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
  iOff(1,0) = nKL

  nJ = nBas_Aux(0)-1
  iOff(2,0) = nJ
  nSize_3C = nJ*nkl
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do klIrrep=0,nIrrep-1
    iOff(3,klIrrep) = nSize_3C

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
          nKL = nKL+nTri_ELem(nK)
        end if

      end do
    end if
    iOff(1,klIrrep) = nKL

    nJ = nBas_Aux(klIrrep)
    if (klIrrep == 0) nJ = nJ-1
    iOff(2,klIrrep) = nJ
    nSize_3C = nSize_3C+nJ*nKL

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end function nSize_3C
