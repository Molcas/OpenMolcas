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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 11, 2022, created this file.               *
!*****************************************************************

! Subroutine relating to generalized 1-e density matrix (GD) called
! in CMSNewton
!     transform GD^KL_tu from leading with state indices
!     to leading with orbital indices in mode 1, and vice
!     versa in mode 2.
subroutine TransposeMat(Matout,Matin,nElem,nRow_in,nCol_in)

implicit none

integer nElem, nRow_in, nCol_in, iRow, iCol, iOff1, iOff2
real*8 Matin(nElem), Matout(nElem)

if (nRow_in*nCol_in /= nElem) then
  write(6,*) 'Error in TransposeMat()'
  write(6,*) 'nRow_in*nCol_in != nElem'
end if

do iCol=1,nCol_in
  iOff1 = (iCol-1)*nRow_in
  do iRow=1,nRow_in
    iOff2 = (iRow-1)*nCol_in
    Matout(iOff2+iCol) = Matin(iOff1+iRow)
  end do
end do

return

end subroutine TransposeMat
