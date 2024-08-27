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
! Copyright (C) 1990, IBM                                              *
!               1991, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
function TstFnc(iCoSet,iIrrep,iBsFnc,nStab)
!***********************************************************************
!                                                                      *
! Object: to establish if a function is a basis function of a          *
!         irreducible representation.                                  *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             September 1991                                           *
!***********************************************************************

use Symmetry_Info, only: nIrrep, iOper, iChTbl

implicit none
logical TstFnc
integer iIrrep, iBsFnc, nStab
integer iCoSet(0:7,0:7), iAcc(0:7)
integer nCoSet, i, j, k, n
integer, external :: iPrmt

TstFnc = .true.
nCoSet = nIrrep/nStab
iAcc(0:nCoSet-1) = 0

#ifdef _DEBUGPRINT_
write(6,*) 'TstFnc'
write(6,*)
write(6,*) 'Coset:'
do i=0,nCoSet-1
  write(6,'(8I4)') (iCoSet(i,j),j=0,nStab-1)
end do

write(6,*)
write(6,*) 'iOper:'
write(6,'(8I4)') (iOper(i),i=0,nIrrep-1)
write(6,*)
write(6,*) 'iBsFnc=',iBsFnc
write(6,*)
write(6,*) 'iChTbl:'
write(6,'(8I4)') (iChTbl(iIrrep,i),i=0,nIrrep-1)
#endif

! Loop over operators

do i=0,nIrrep-1

  ! Find index of the generated center

  n = -1
  do j=0,nCoSet-1
    if (n >= 0) cycle
    do k=0,nStab-1
      if (iOper(i) == iCoSet(j,k)) n = j
    end do
  end do

  if ((n < 0) .or. (n > nCoSet-1)) then
    call WarningMessage(2,'TstFnc: n < 0 .or. n > nCoSet-1')
    write(6,*) ' Coset index',n,' is wrong!'
    call Abend()
  end if

  iAcc(n) = iAcc(n)+iChTbl(iIrrep,i)*iPrmt(i,iBsFnc)

end do
do i=0,nCoSet-1
  if (iAcc(i) == 0) TstFnc = .false.
end do

return

end function TstFnc
