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

! This is just an encapsulation of the common blocks in
! src/Include/rasdim.fh
! src/Include/general.fh
! into a data module

module general_data

use Definitions, only: wp, iwp

implicit none

#include "rasdim.fh"
Private :: MaxBfn,MaxBfn_Aux,MxAO,mxAtom,mxroot,mxNemoAtom,Mxdbsc,lCache,mxact,mxina,mxbas,mxOrb,mxSym,mxGAS, &
           LENIN,LENIN1,LENIN2,LENIN3,LENIN4,LENIN5,LENIN6,LENIN8
Private :: mxRef,mxIter,mxCiIt,mxSxIt,mxTit

#include "general.fh"

integer(kind=iwp), allocatable :: CleanMask(:)
real(kind=wp), allocatable :: CRPROJ(:), CRVEC(:)

end module general_data
