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
! Copyright (C) 1992, Roland Lindh                                     *
!               1995, Martin Schuetz                                   *
!***********************************************************************

subroutine AlloK2_Funi(nr_of_Densities)
!***********************************************************************
!                                                                      *
!  Object: Allocate space for K2 entities.                             *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, Sweden. November '92                 *
!             Martin Schuetz, Dept. of Theoretical Chemistry,          *
!             University of Lund, Sweden. Jun '95                      *
!***********************************************************************

use iSD_data, only: iSD
use k2_arrays, only: MaxDe, nDeDe_DFT
use Symmetry_Info, only: nIrrep
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nr_of_Densities
integer(kind=iwp) :: iAO, iBas, iCmp, iDeSiz, iS, iShell, iSmLbl, jAO, jBas, jCmp, jS, jShell, nSkal, nSO
integer(kind=iwp), external :: MemSO1

call Nr_Shells(nSkal)

! determine memory size nDeDe, MaxDe, and MaxDRC
nDeDe_DFT = 0
MaxDe = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Double loop over shells. These loops decide the integral type

do iS=1,nSkal
  !iAng = iSD(1,iS)
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  !iPrim = iSD(5,iS)
  iAO = iSD(7,iS)
  !mdci = iSD(10,iS)
  iShell = iSD(11,iS)

  do jS=1,iS
    !jAng = iSD(1,jS)
    jCmp = iSD(2,jS)
    jBas = iSD(3,jS)
    !jPrim = iSD(5,jS)
    jAO = iSD(7,jS)
    !mdcj = iSD(10,jS)
    jShell = iSD(11,jS)

    !iDeSiz = 1+iPrim*jPrim+(iBas*jBas+1)*iCmp*jCmp
    iDeSiz = iBas*jBas*iCmp*jCmp
    MaxDe = max(MaxDe,iDeSiz)
    iSmLbl = 1
    nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
    if (nSO > 0) nDeDe_DFT = nDeDe_DFT+nr_of_Densities*iDeSiz*nIrrep

  end do
end do

return

end subroutine AlloK2_Funi
