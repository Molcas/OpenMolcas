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

module RASSI_AUX

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: ipglob, mTRA, NASHT_Save
logical(kind=iwp) :: AO_Mode = .false.
integer(kind=iwp), allocatable :: jDisk_TDM(:,:), JOB_INDEX(:), TocM(:)
real(kind=wp), allocatable :: CMO1(:), CMO2(:), DMAB(:)
integer(kind=iwp), parameter :: MULT = 37, NHASH = 997

public :: AO_Mode, CMO1, CMO2, DMAB, iDisk_TDM, ipglob, jDisk_TDM, JOB_INDEX, mTRA, MULT, NASHT_Save, NHASH, TocM

contains

function iDisk_TDM(I,J,K)

  use Index_Functions, only: nTri_Elem

  integer(kind=iwp) :: iDisk_TDM
  integer(kind=iwp) :: I, J, K
  integer(kind=iwp) :: I_Max, J_Min

  I_Max = max(I,J)
  J_Min = min(I,J)
  iDisk_TDM = jDisk_TDM(K,nTri_Elem(I_Max-1)+J_Min)

end function iDisk_TDM

end module RASSI_AUX
