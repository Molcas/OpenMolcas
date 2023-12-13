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

use Data_Structures, only: Alloc1DArray_Type, Alloc2DArray_Type
use Definitions, only: wp, iwp

implicit none
private

#include "LenIn.fh"

real(kind=wp) :: EneV
character(len=180) :: Title
character(len=8) :: Method
integer(kind=iwp), allocatable :: iAtomType(:), iAtomPar(:), nAtomPBas(:), iAtPrTab(:,:)
real(kind=wp), allocatable :: Cor(:,:,:), CordMltPl(:,:), Frac(:,:), AtPol(:,:), AtBoPol(:,:), Qnuc(:)
character(len=LenIn), allocatable :: Labe(:)
character(len=LenIn*2+2), allocatable :: Cen_Lab(:)
logical(kind=iwp), allocatable :: BondMat(:,:)
type(Alloc1DArray_Type), allocatable :: AtBoMltPlTot(:), AtMltPlTot(:)
type(Alloc2DArray_Type), allocatable :: AtBoMltPl(:), AtBoMltPlCopy(:), AtMltPl(:), MltPl(:)

public :: AtBoMltPl, AtBoMltPlCopy, AtBoMltPlTot, AtBoPol, AtMltPl, AtMltPlTot, AtPol, BondMat, Cen_Lab, Cor, CordMltPl, EneV, &
          Frac, iAtomPar, iAtomType, iAtPrTab, Labe, Method, MltPl, nAtomPBas, Qnuc, Title

end module MPProp_globals
