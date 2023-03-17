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

module write_pdft_job
  implicit none
  private
  integer :: iwjob

  ! read from / write to HDF5 file
  logical :: hasHDF5ref=.false.

  ! reference wave function is of MPS type (aka "DMRG wave function")
  logical :: hasMPSref=.false.

  public :: iwjob, hasHDF5ref, hasMPSref

end module write_pdft_job