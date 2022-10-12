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
! This routine gets scalar double data from the runfile.               *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************
!  Get_dScalar
!
!> @brief
!>   Get scalar data form runfile
!> @author Per-Olof Widmark
!>
!> @details
!> This routine gets scalar double data from the runfile.
!>
!> @param[in]  Label Name of field
!> @param[out] Data  Data to get from runfile
!>
!> @see ::Put_dScalar
!***********************************************************************

subroutine Get_dScalar(Label,rData)

use RunFile_data, only: DS_cache, nTocDS, lw, num_DS_init
use Definitions, only: wp, iwp

implicit none
character(len=*), intent(in) :: Label
real(kind=wp), intent(out) :: rData
character(len=lw) :: CmpLab
integer(kind=iwp) :: i

CmpLab = Label
call UpCase(CmpLab)

do i=1,num_DS_init
  if (DS_cache(i)%lab == CmpLab) then
    rData = DS_cache(i)%val
    return
  end if
end do

call Get_dScalar_(Label,rData)
num_DS_init = num_DS_init+1

if (num_DS_init > nTocDS) then
# ifdef _DEBUGPRINT_
  do i=1,nTocDS
    write(u6,*) DS_cache(i)%lab,DS_cache(i)%val,CmpLab
  end do
# endif
  call Abend()
end if

DS_cache(num_DS_init)%lab = CmpLab
DS_cache(num_DS_init)%val = rData

return

end subroutine Get_dScalar
