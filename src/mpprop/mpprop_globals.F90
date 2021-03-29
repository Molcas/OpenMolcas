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

integer(kind=iwp), parameter :: mxAtomMP = 500, mxCen = mxAtomMP*(mxAtomMP+1)/2, mxMltPl = 10
integer(kind=iwp) :: iAtPrTab(maxbfn,mxAtomMP), iAtomType(mxAtomMP), iAtomPar(mxAtomMP), nAtomPBas(mxAtomMP), iBondPar(mxCen), &
                     NUB(mxAtomMP), NBI(mxAtomMP,mxAtomMP), iAtMltPlAd(0:mxMltPl), iAtBoMltPlAd(0:mxMltPl), iMltPlAd(0:mxMltPl), &
                     iAtPolAd, iAtBoPolAd, iQnuc, iAtMltPlTotAd(0:mxMltPl), iAtBoMltPlTotAd(0:mxMltPl), iAtBoMltPlAdCopy(0:mxMltpl)
real(kind=wp) :: Cor(3,mxAtomMP,mxAtomMP), Frac(mxAtomMP,mxAtomMP), CordMltpl(3,0:99), EneV
character(len=LenIn) :: Labe(mxAtomMP)
character(len=LenIn*2+2) :: Cen_Lab(mxCen)
character(len=180) :: Title
character(len=8) :: Method
logical(kind=iwp) :: BondMat(mxAtomMP,mxAtomMP)

public :: BondMat, Cen_Lab, Cor, CordMltPl, EneV, Frac, iAtBoMltPlAd, iAtBoMltPlAdCopy, iAtBoMltPlTotAd, iAtBoPolAd, iAtMltPlAd, &
          iAtMltPlTotAd, iAtomPar, iAtomType, iAtPolAd, iAtPrTab, iBondPar, iMltPlAd, iQnuc, Labe, Method, mxAtomMP, mxCen, &
          mxMltPl, nAtomPBas, NBI, NUB, Title

end module MPProp_globals
