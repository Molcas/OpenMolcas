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

module ipPage

integer, parameter :: Max_CI_Vectors = 40
integer, parameter :: On_Disk = 0, In_Memory = 1, Null_Vector = 2
integer, parameter :: dWrite = 0, write = 1, read = 2

! ip_Mem : memory pointer
! n  : Length of CI-vector
! ida: disk address

integer :: n(0:Max_CI_Vectors)
integer :: ida(0:Max_CI_Vectors)
integer :: Status(0:Max_CI_Vectors)

integer :: n_CI_Vectors = 0
integer :: iDisk_Addr_End = 0
integer :: Lu_ip = -99
logical :: DiskBased = .false.

type Vector
  real*8, allocatable :: Vec(:)
end type Vector

type(Vector) :: W(0:Max_CI_Vectors)

end module ipPage
