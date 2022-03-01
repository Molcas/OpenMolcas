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

subroutine ccsort_rea(lun,length,vector)
! this routine reads length-R8 numbers from open unformatted file
! with number lun form the given position as one record
!
! lun    - Logical unit number of file, where mediate is stored (Input)
! length - # of R8 numbers to be read (Input)
! vector - space, where numbers are stored after reading (Output)

integer lun, length, i
real*8 vector(1:length)

read(lun) (vector(i),i=1,length)

return

end subroutine ccsort_rea
