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

! A type for multipole array data, where each element
! is an array of different size (number of components)
type MltPlArr
  real(kind=wp), allocatable :: M(:,:)
end type MltPlArr

integer(kind=iwp), parameter :: mxMltPl = 10
real(kind=wp) :: EneV
character(len=180) :: Title
character(len=8) :: Method
integer(kind=iwp), allocatable :: iAtomType(:), iAtomPar(:), nAtomPBas(:), iAtPrTab(:,:)
real(kind=wp), allocatable :: Cor(:,:,:), CordMltPl(:,:), Frac(:,:), AtPol(:,:), AtBoPol(:,:), Qnuc(:)
character(len=LenIn), allocatable :: Labe(:)
character(len=LenIn*2+2), allocatable :: Cen_Lab(:)
logical(kind=iwp), allocatable :: BondMat(:,:)
type(MltPlArr), allocatable :: AtBoMltPl(:), AtBoMltPlCopy(:), AtBoMltPlTot(:), AtMltPl(:), AtMltPlTot(:), MltPl(:)

public :: AtPol, AtBoPol, BondMat, Cen_Lab, Cor, CordMltPl, EneV, Frac, iAtomPar, iAtomType, iAtPrTab, Labe, Method, mxMltPl, &
          nAtomPBas, Qnuc, Title

public :: AtBoMltPl, AtBoMltPlCopy, AtBoMltPlTot, AtMltPl, AtMltPlTot, MltPl

end module MPProp_globals
