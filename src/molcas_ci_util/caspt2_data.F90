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

module caspt2_data

use Definitions, only: wp, iwp

implicit none
private

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"

real(kind=wp), allocatable, target:: CMO_Internal(:)
real(kind=wp), pointer:: CMO(:)
real(kind=wp), allocatable:: FIMO(:)
real(kind=wp), allocatable:: FAMO(:)
real(kind=wp), allocatable:: FIFA(:)
real(kind=wp), allocatable:: HONE(:)
real(kind=wp), allocatable:: DREF(:)
real(kind=wp), allocatable:: PREF(:)
real(kind=wp), allocatable:: DMIX(:,:)
real(kind=wp), allocatable:: DWGT(:,:)
real(kind=wp), allocatable:: CMOPT2(:)
real(kind=wp), allocatable:: TAT(:)
integer(kind=iwp) :: NTAT=0
real(kind=wp), allocatable:: TORB(:)
integer(kind=iwp) :: NTORB=0

public :: CMO, CMO_Internal, jState, mState, nActEl, nG3, FIMO, FAMO, FIFA, HONE, DREF, PREF, DMIX, DWGT, CMOPT2
public :: TAT, NTAT, TORB, NTORB

end module caspt2_data
