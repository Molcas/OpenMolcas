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
      SubRoutine Free_iSD()
      use iSD_data, only: iSD, nSkal_iSD
      use stdalloc, only: mma_deallocate
      Implicit None
!
      If (Allocated(iSD)) Call mma_deallocate(iSD)
      nSkal_iSD=0
!
      Return
      End SubRoutine Free_iSD
