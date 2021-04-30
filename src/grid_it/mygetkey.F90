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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************
!  MyGetKey
!
!> @brief
!>   General purpose routine to read arbitrary data from user input.
!> @author V. Veryazov
!>
!> @details
!> The routine read a line from unit \p InUnit (ignoring molcas comments
!> and blank lines), and return a value or an array.
!>
!> Parameter \p What specifies the type of data:
!> * ``I`` -- read an integer and return it as \p IValue
!> * ``R`` -- read a real and return it as \p RValue
!> * ``S`` -- read a string and return it as \p SValue
!> * ``U`` -- recognize integer/real/string and return corresponding value
!> * ``A`` -- read integer array and return \p IArr
!> * ``D`` -- read real array and return \p RArr
!>
!> @param[in]     InUnit   Unit number
!> @param[in,out] What     Type of input
!> @param[out]    IValue   Integer Value
!> @param[out]    RValue   Real Value
!> @param[out]    SValue   String value
!> @param[in]     N        Size of array
!> @param[out]    IArr     Integer Array
!> @param[out]    RArr     Real Array

function MyGetKey(InUnit,What,IValue,RValue,SValue,N,IArr,RArr)
!***********************************************************************
! Adapted from SAGIT to work with OpenMolcas (October 2020)            *
!***********************************************************************

#include "intent.fh"

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: MyGetKey
integer(kind=iwp), intent(in) :: InUnit, N
character, intent(inout) :: What
integer(kind=iwp), intent(out) :: IValue
real(kind=wp), intent(out) :: RValue
character(len=*), intent(out) :: SValue
integer(kind=iwp), intent(_OUT_) :: IArr(*)
real(kind=wp), intent(_OUT_) :: RArr(*)
integer(kind=iwp) :: i, iptr, istatus
logical(kind=iwp) :: retry
character(len=120) :: KWord

MyGetKey = 0
iptr = 1
i = 1
retry = .true.
do while (retry)
  read(InUnit,'(A)',iostat=istatus) KWord
  if (istatus /= 0) then
    MyGetKey = 1
    exit
  end if
  if ((KWord(1:1) == '*') .or. (KWord == ' ')) cycle
  retry = .false.
  call UpCase(KWord)
  select case (What)
    case ('I')
      read(KWord,*,iostat=istatus) IValue
      if (istatus > 0) MyGetKey = 1
    case ('R')
      read(KWord,*,iostat=istatus) RValue
      if (istatus > 0) MyGetKey = 1
    case ('A')
      read(KWord,*,iostat=istatus) (IArr(i),i=iptr,N)
      if (istatus > 0) then
        MyGetKey = 1
      else if (istatus < 0) then
        iptr = i
        retry = .true.
      end if
    case ('D')
      read(KWord,*,iostat=istatus) (RArr(i),i=iptr,N)
      if (istatus > 0) then
        MyGetKey = 1
      else if (istatus < 0) then
        iptr = i
        retry = .true.
      end if
    case ('S')
      call NoBlanks(SValue,KWord)
    case ('U')
      read(KWord,*,iostat=istatus) IValue
      if (istatus > 0) then
        read(KWord,*,iostat=istatus) RValue
        if (istatus > 0) then
          call NoBlanks(SValue,KWord)
          What = 'S'
        else
          What = 'R'
        end if
      else
        What = 'I'
      end if
    case default
  end select
end do

return

end function MyGetKey
