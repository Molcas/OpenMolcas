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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine CD_Tester_Final(irc,NumCho,n,Thr,Err,Verbose)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: NumCho, n
real(kind=wp), intent(in) :: Thr, Err(6)
logical(kind=iwp), intent(in) :: Verbose
character(len=*), parameter :: SecNam = 'CD_Tester_Final'

irc = 0

if (Verbose) then
  write(u6,*)
  write(u6,*) 'Final results from ',SecNam,':'
  write(u6,*) 'Matrix dimension: ',n
  write(u6,*) 'Number of vecs. : ',NumCho
  write(u6,*) 'Threshold       : ',Thr
  write(u6,*) 'Min. Diag. err. : ',Err(1)
  write(u6,*) 'Max. Diag. err. : ',Err(2)
  write(u6,*) 'RMS  Diag. err. : ',Err(3)
  write(u6,*) 'Min. Matr. err. : ',Err(4)
  write(u6,*) 'Max. Matr. err. : ',Err(5)
  write(u6,*) 'RMS  Matr. err. : ',Err(6)
end if

if ((NumCho < 0) .or. (NumCho > n)) then
  irc = -1
  if (Verbose) write(u6,*) '>>> NumCho out of bounds!'
  return
end if

if (abs(Err(1)) > Thr) then
  irc = irc+1
  if (Verbose) write(u6,*) '>>> LARGE MINIMUM DIAGONAL ERROR: ',Err(1)
end if
if (abs(Err(2)) > Thr) then
  irc = irc+1
  if (Verbose) write(u6,*) '>>> LARGE MAXIMUM DIAGONAL ERROR: ',Err(2)
end if
if (abs(Err(3)) > Thr) then
  irc = irc+1
  if (Verbose) write(u6,*) '>>> LARGE RMS     DIAGONAL ERROR: ',Err(3)
end if
if (abs(Err(4)) > Thr) then
  irc = irc+1
  if (Verbose) write(u6,*) '>>> LARGE MINIMUM MATRIX   ERROR: ',Err(4)
end if
if (abs(Err(5)) > Thr) then
  irc = irc+1
  if (Verbose) write(u6,*) '>>> LARGE MAXIMUM MATRIX   ERROR: ',Err(5)
end if
if (abs(Err(6)) > Thr) then
  irc = irc+1
  if (Verbose) write(u6,*) '>>> LARGE RMS     MATRIX   ERROR: ',Err(6)
end if

if (Verbose) call xFlush(u6)

end subroutine CD_Tester_Final
