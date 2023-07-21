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

implicit none
integer irc, NumCho, n
real*8 Thr, Err(6)
logical Verbose
character*15 SecNam
parameter(SecNam='CD_Tester_Final')

irc = 0

if (Verbose) then
  write(6,*)
  write(6,*) 'Final results from ',SecNam,':'
  write(6,*) 'Matrix dimension: ',n
  write(6,*) 'Number of vecs. : ',NumCho
  write(6,*) 'Threshold       : ',Thr
  write(6,*) 'Min. Diag. err. : ',Err(1)
  write(6,*) 'Max. Diag. err. : ',Err(2)
  write(6,*) 'RMS  Diag. err. : ',Err(3)
  write(6,*) 'Min. Matr. err. : ',Err(4)
  write(6,*) 'Max. Matr. err. : ',Err(5)
  write(6,*) 'RMS  Matr. err. : ',Err(6)
end if

if ((NumCho < 0) .or. (NumCho > n)) then
  irc = -1
  if (Verbose) write(6,*) '>>> NumCho out of bounds!'
  return
end if

if (abs(Err(1)) > Thr) then
  irc = irc+1
  if (Verbose) write(6,*) '>>> LARGE MINIMUM DIAGONAL ERROR: ',Err(1)
end if
if (abs(Err(2)) > Thr) then
  irc = irc+1
  if (Verbose) write(6,*) '>>> LARGE MAXIMUM DIAGONAL ERROR: ',Err(2)
end if
if (abs(Err(3)) > Thr) then
  irc = irc+1
  if (Verbose) write(6,*) '>>> LARGE RMS     DIAGONAL ERROR: ',Err(3)
end if
if (abs(Err(4)) > Thr) then
  irc = irc+1
  if (Verbose) write(6,*) '>>> LARGE MINIMUM MATRIX   ERROR: ',Err(4)
end if
if (abs(Err(5)) > Thr) then
  irc = irc+1
  if (Verbose) write(6,*) '>>> LARGE MAXIMUM MATRIX   ERROR: ',Err(5)
end if
if (abs(Err(6)) > Thr) then
  irc = irc+1
  if (Verbose) write(6,*) '>>> LARGE RMS     MATRIX   ERROR: ',Err(6)
end if

if (Verbose) call xFlush(6)

end subroutine CD_Tester_Final
