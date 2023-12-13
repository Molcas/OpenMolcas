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
! Copyright (C) 2009, Giovanni Ghigo                                   *
!***********************************************************************

!  InterSystem Crossing rate evaluation: Reduction of States
!  Author: Giovanni Ghigo
!          Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
!          07 Jan-09 - XX Jan-09

subroutine LogEVec(iPrint,nOsc,max_nOrd,minQ,nMaxQ,nMat,lVec,nYes)
! Generate Logical Vector of useful States

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, nOsc, max_nOrd, minQ, nMaxQ(nOsc), nMat(0:max_nOrd,nOsc)
integer(kind=iwp), intent(out) :: lVec(0:max_nOrd), nYes
integer(kind=iwp) :: iOrd, iOsc, nSumQ

if (iPrint >= 3) then
  write(u6,*) ' Original number of States=',max_nOrd+1
end if

do iOrd=0,max_nOrd
  lVec(iOrd) = 1
  nSumQ = 0
  do iOsc=1,nOsc
    if (nMat(iOrd,iOsc) > nMaxQ(iOsc)) lVec(iOrd) = 0
    nSumQ = nSumQ+nMat(iOrd,iOsc)
  end do
  if (nSumQ < minQ) lVec(iOrd) = 0
end do
nYes = 0
do iOrd=0,max_nOrd
  nYes = nYes+lVec(iOrd)
end do

if (iPrint >= 3) then
  write(u6,*) ' Selected number of States=',nYes
end if

return

end subroutine LogEVec
