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

module MPProp_globals

use Definitions, only: wp, iwp

implicit none
private

#include "Molcas.fh"

integer(kind=iwp), parameter :: mxMltPl = 10
integer(kind=iwp) :: iAtMltPlAd(0:mxMltPl), iAtBoMltPlAd(0:mxMltPl), iMltPlAd(0:mxMltPl), iAtPolAd, iAtBoPolAd, iQnuc, &
                     iAtMltPlTotAd(0:mxMltPl), iAtBoMltPlTotAd(0:mxMltPl), iAtBoMltPlAdCopy(0:mxMltpl)
real(kind=wp) :: CordMltPl(3,0:mxMltPl), EneV
character(len=180) :: Title
character(len=8) :: Method
integer(kind=iwp), allocatable :: iAtomType(:), iAtomPar(:), nAtomPBas(:), iAtPrTab(:,:)
real(kind=wp), allocatable :: Cor(:,:,:), Frac(:,:)
character(len=LenIn), allocatable :: Labe(:)
character(len=LenIn*2+2), allocatable :: Cen_Lab(:)
logical(kind=iwp), allocatable :: BondMat(:,:)

public :: BondMat, Cen_Lab, Cor, CordMltPl, EneV, Frac, iAtBoMltPlAd, iAtBoMltPlAdCopy, iAtBoMltPlTotAd, iAtBoPolAd, iAtMltPlAd, &
          iAtMltPlTotAd, iAtomPar, iAtomType, iAtPolAd, iAtPrTab, iMltPlAd, iQnuc, Labe, Method, mxMltPl, nAtomPBas, Title

end module MPProp_globals
