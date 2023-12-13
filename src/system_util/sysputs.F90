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
! Copyright (C) 2001, Valera Veryazov                                  *
!***********************************************************************
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!       general purpose routine for printing nice messages             *
!       with simple reformatting:                                      *
!           \n and long string converted to new line                   *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     V.Veryazov University of Lund, 2001                              *
!                                                                      *
!***********************************************************************

subroutine SysPuts(str,str1,str2)

use Definitions, only: iwp

implicit none
character(len=*), intent(in) :: str, str1, str2
integer(kind=iwp) :: i, icr, icr1, ii, iii, ij, iLongEnough, ipos, iTooLong, j, mlen, mleni
character(len=len(str)+len(str1)+len(str2)) :: Junk

! because of bug in g77 we can't just concatenate strings and
! had to have limited length of the string
iTooLong = 60
iLongEnough = 50
Junk = str//trim(str1)//trim(str2)
mlen = len_trim(Junk)
mleni = mlen
ipos = 1
! check is '\n' == <CR>?
icr = len('\n')-1
icr1 = 0
do
  j = 100000
  ii = index(Junk(ipos:mleni),'\n')
  if (ii > 0) j = min(j,ii)
  iii = index(Junk(ipos:mleni),';')
  if (iii > 0) j = min(j,iii)
  if ((j == ii) .and. (ii > 0)) icr1 = icr
  if ((j == iii) .and. (iii > 0)) icr1 = 0
  if (j == 100000) j = 0
  i = j
  if ((i > iTooLong) .or. ((i == 0) .and. (mlen > iTooLong))) then
    ij = index(Junk(ipos+iLongEnough:mleni),' ')
    if (ij == 0) then
      call SysDumpStr(Junk(ipos:mleni))
      if (i == 0) return
    else
      ij = ij+iLongEnough-1
      call SysDumpStr(Junk(ipos:ipos+ij))
    end if
    ipos = ipos+ij+1
    mlen = mlen-ij-1
    cycle
  end if
  if (i == 0) then
    call SysDumpStr(Junk(ipos:mleni))
    return
  end if
  call SysDumpStr(Junk(ipos:ipos+i-2))
  ipos = ipos+i+icr1
  mlen = mlen-i-icr1
  if (mlen <= 0) exit
end do

return

end subroutine SysPuts
