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
! Copyright (C) Per Ake Malmqvist                                      *
!               Jesper Wisborg Krogh                                   *
!***********************************************************************
!  Chk4NAN
!
!> @brief
!>   Check whether an array contains NANs
!> @author Per &Aring;ke Malmqvist
!> @modified_by Jesper Wisborg Krogh
!>
!> @details
!> The routine checks the elements of the input array for NANs. If any are found,
!> up to 100 of these will be listed. Upon return \p iErr contains the number of
!> NANs found.
!>
!> @param[in]  nDim  Total dimension of array
!> @param[in]  Array Array to be checked
!> @param[out] Ierr  Return code
!***********************************************************************

subroutine Chk4NAN(nDim,Array,Ierr)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nDim
real(kind=wp), intent(in) :: Array(nDim)
integer(kind=iwp), intent(out) :: Ierr
integer(kind=iwp) :: I, iCount
real(kind=wp) :: CheckSum
character(len=16) :: str16

ICOUNT = 0
CHECKSUM = sum(ARRAY)
write(STR16,'(G16.7)') CHECKSUM
call NORMAL(STR16)
if (STR16(1:1) == 'N') then
  write(u6,*) '!!! WARNING !!!'
  write(u6,*) 'NANs encountered'
  write(u6,*)
  write(u6,*) ' The numbers in the array will now be checked.'
  write(u6,*) ' There are ',NDIM,' elements.'
  do I=1,nDim
    write(STR16,'(G16.7)') ARRAY(I)
    call NORMAL(STR16)
    if (STR16(1:1) == 'N') then
      ICOUNT = ICOUNT+1
      if (iCount <= 100) write(u6,*) ' Element nr.',I,' is ',ARRAY(I)
    end if
  end do
  if (ICOUNT > 100) write(u6,*) ' ...too many. I give up here.'
  write(u6,*) 'There were a total of ',iCount,' NANs'
end if

iErr = iCount

return

end subroutine Chk4NAN
