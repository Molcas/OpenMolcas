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

#include "caspt2.fh"
#include "pt2_guga.fh"

! UNIT numbers:
! IDCIEX, IDTCEX, LUCIEX, LUDMAT, LUDRA, LUDRATOT, LUH0T, LUHLF1, LUHLF2, LUHLF3, LUINTA, LUINTM, LUONEM, LURHS, LUSBT, LUSOLV

integer(kind=iwp) :: IDCIEX, IDTCEX, LUCIEX, LUDMAT, LUDRA, LUDRATOT, LUH0T(4), LUHLF1, LUHLF2, LUHLF3, LUINTA, LUINTM, LUONEM, &
                     LURHS(8), LUSBT, LUSOLV, NCMO = 0, NDREF = 0, NPREF = 0, NTAT = 0, NTORB = 0
integer(kind=iwp), allocatable :: IDSCT(:), LISTS(:)
real(kind=wp), allocatable:: CMOPT2(:), DMIX(:,:), DREF(:), DWGT(:,:), FAMO(:), FIFA(:), FIMO(:), HONE(:), PREF(:), TAT(:), TORB(:)
real(kind=wp), allocatable, target:: CMO_Internal(:)
real(kind=wp), pointer:: CMO(:)

public :: CMO, CMO_Internal, CMOPT2, DMIX, DREF, DWGT, FAMO, FIFA, FIMO, HONE, IDCIEX, IDSCT, IDTCEX, jState, LISTS, LUCIEX, &
          LUDMAT, LUDRA, LUDRATOT, LUH0T, LUHLF1, LUHLF2, LUHLF3, LUINTA, LUINTM, LUONEM, LURHS, LUSBT, LUSOLV, mState, nActEl, &
          NCMO, NDREF, nG3, NPREF, NTAT, NTORB, PREF, TAT, TORB

end module caspt2_data
