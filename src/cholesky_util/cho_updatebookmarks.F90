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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_UpdateBookmarks(iRS,nSym,nRS,nVec,delta,nBkmVec,BkmThr)

implicit none
integer iRS
integer nRS
integer nSym
integer nVec(nSym)
real*8 delta(nSym)
integer nBkmVec(nSym,nRS)
real*8 BkmThr(nSym,nRS)
integer iSym

do iSym=1,nSym
  nBkmVec(iSym,iRS) = nVec(iSym)
end do

do iSym=1,nSym
  BkmThr(iSym,iRS) = delta(iSym)
end do

end subroutine Cho_UpdateBookmarks
