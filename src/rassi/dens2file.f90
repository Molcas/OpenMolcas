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
  subroutine dens2file(array1,array2,array3,adim,lu,adr,iOpt)
  implicit none

  integer, intent(in) :: adim, lu, iOpt
  integer, intent(inout) :: adr
  real*8 , intent(inout) :: array1(adim),array2(adim),array3(adim)

    call ddafile(lu,iOpt,array1,adim,adr)
    call ddafile(lu,iOpt,array2,adim,adr)
    call ddafile(lu,iOpt,array3,adim,adr)
!
  end subroutine dens2file
