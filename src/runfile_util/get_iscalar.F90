!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2003, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine gets scalar integer data from the runfile.              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Get_iScalar
!
!> @brief
!>   Get scalar data from runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine gets scalar integer data from the runfile.
!>
!> @param[in]  Label Name of field
!> @param[out] Data  Data to get from runfile
!>
!> @see ::Put_iScalar
!***********************************************************************

subroutine Get_iScalar(Label,iData)

use RunFile_data, only: IS_cache, lw, nTocIS, num_IS_init
use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(out) :: iData
integer(kind=iwp) :: i
character(len=lw) :: CmpLab

CmpLab = Label
call UpCase(CmpLab)

do i=1,num_IS_init
  if (IS_cache(i)%lab == CmpLab) then
    iData = IS_cache(i)%val
    return
  end if
end do

call Get_iScalar_(Label,iData)
num_IS_init = num_IS_init+1

if (num_IS_init > nTocIS) then
# ifdef _DEBUGPRINT__
  do i=1,nTocIS
    write(u6,*) IS_cache(i)%lab,IS_cache(i)%val,CmpLab
  end do
# endif
  call Abend()
end if

IS_cache(num_IS_init)%lab = CmpLab
IS_cache(num_IS_init)%val = iData

return

end subroutine Get_iScalar
