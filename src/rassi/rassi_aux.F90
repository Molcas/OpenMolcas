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

implicit none

logical :: AO_Mode = .false.
integer, allocatable :: TocM(:), jDisk_TDM(:,:), JOB_INDEX(:)
real*8, allocatable :: CMO1(:), CMO2(:), DMAB(:)
integer NASHT_Save, mTRA
integer :: JOB1_old = -1, JOB2_old = -1
integer :: ipglob

contains

integer function iDisk_TDM(I,J,K)

  integer I, J, I_Max, J_Min, K

  I_Max = max(I,J)
  J_Min = min(I,J)
  iDisk_TDM = jDisk_TDM(K,I_Max*(I_Max-1)/2+J_Min)

end function iDisk_TDM

end module RASSI_AUX
