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

subroutine diisra(wrk,wrksize,diispoint,num,mapd1,mapi1,poss10,mapd2,mapi2,poss20,mapd3,mapi3,poss30,mapd4,mapi4,poss40)
! this routine reads num vectors of amplitudes
! form prepaired diispoint(1) - diispoint(num) files
!
! diispoint - stack of lun numbers (I)
! num       - number of vertors to be read (1-4) (I)
! mapd1     - direct matrix of vector 1 (I)
! mapi1     - inverse matrix of vector 1 (I)
! poss10    - initial position of vector 1 (I)
! mapd2     - direct matrix of vector 2 (I)
! mapi2     - inverse matrix of vector 2 (I)
! poss20    - initial position of vector 2 (I)
! mapd3     - direct matrix of vector 3 (I)
! mapi3     - inverse matrix of vector 3 (I)
! poss30    - initial position of vector 3 (I)
! mapd4     - direct matrix of vector 4 (I)
! mapi4     - inverse matrix of vector 4 (I)
! poss40    - initial position of vector 4 (I)
! if there is less than 4 vectors required use any map's and poss's

#include "wrk.fh"
integer diispoint(1:4)
integer num
integer mapd1(0:512,1:6)
integer mapd2(0:512,1:6)
integer mapd3(0:512,1:6)
integer mapd4(0:512,1:6)
integer mapi1(1:8,1:8,1:8)
integer mapi2(1:8,1:8,1:8)
integer mapi3(1:8,1:8,1:8)
integer mapi4(1:8,1:8,1:8)
integer poss10, poss20, poss30, poss40
! help variables
integer lun, rc

if (num == 1) then
  lun = diispoint(1)
  call getmediate(wrk,wrksize,lun,poss10,mapd1,mapi1,rc)
else if (num == 2) then
  lun = diispoint(1)
  call getmediate(wrk,wrksize,lun,poss10,mapd1,mapi1,rc)
  lun = diispoint(2)
  call getmediate(wrk,wrksize,lun,poss20,mapd2,mapi2,rc)
else if (num == 3) then
  lun = diispoint(1)
  call getmediate(wrk,wrksize,lun,poss10,mapd1,mapi1,rc)
  lun = diispoint(2)
  call getmediate(wrk,wrksize,lun,poss20,mapd2,mapi2,rc)
  lun = diispoint(3)
  call getmediate(wrk,wrksize,lun,poss30,mapd3,mapi3,rc)
else if (num == 4) then
  lun = diispoint(1)
  call getmediate(wrk,wrksize,lun,poss10,mapd1,mapi1,rc)
  lun = diispoint(2)
  call getmediate(wrk,wrksize,lun,poss20,mapd2,mapi2,rc)
  lun = diispoint(3)
  call getmediate(wrk,wrksize,lun,poss30,mapd3,mapi3,rc)
  lun = diispoint(4)
  call getmediate(wrk,wrksize,lun,poss40,mapd4,mapi4,rc)
end if

return

end subroutine diisra
