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
!  StrToChk
!
!> @brief
!>   Compute a checksum for a character string. Used for example to compute checksum for basis set
!> @author Per-Olof Widmark
!>
!> @details
!> This routine takes a character string and compute an
!> integer number. This number is not necessarily unique
!> but the chance that two different string will have the
!> same checksum is small. Whitespace is ignored so that
!> '``HelloDolly``' and '``Hello Dolly``'
!> will have the same checksum.
!> By default (\p iOpt = ``0``) the checksum is case insensitive,
!> but can be made case sensitive by setting \p iOpt = ``1``.
!>
!> @param[in]  String String to be checksummed
!> @param[out] Chk    Checksum
!> @param[in]  iOpt   Bitswitch, ``1`` for case sensitive checksum
!***********************************************************************
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
! Written: August 2003                                                 *
!                                                                      *
!***********************************************************************

subroutine StrToChk(String,Chk,iOpt)

!----------------------------------------------------------------------*
! Arguments                                                            *
!----------------------------------------------------------------------*
character*(*) String
integer Chk
integer iOpt
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
logical sensitive
integer i, k, m, n

!----------------------------------------------------------------------*
! Process option switch                                                *
!----------------------------------------------------------------------*
if (btest(iOpt,0)) then
  sensitive = .true.
else
  sensitive = .false.
end if
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
n = 0
k = 1
do i=1,len(String)
  k = mod(k+12,17)+1
  m = ichar(String(i:i))
  if (m == 32) cycle
  if (m == 9) cycle
  if (.not. sensitive) then
    if ((m >= ichar('a')) .and. (m <= ichar('z'))) then
      m = m-ichar('a')+ichar('A')
    end if
  end if
  n = n+k*m
end do
Chk = n

!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return

end subroutine StrToChk
