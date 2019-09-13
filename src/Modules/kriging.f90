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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************
      module kriging
        real*8, allocatable, protected :: x(:,:), y(:), dy(:)
        contains
!
      Subroutine Setup_Kriging(nPoints,nInter,x_,dy_,y_)
!
        Real*8 x_(nInter,nPoints),dy_(nInter,nPoints),y_(nPoints)
!
!x is the n-dimensional internal coordinates
          x(:,:) = x_(:,:)
          ! write(6,*) 'x',x
!y is the energy
          y(:) = y_(:)
          ! write(6,*) 'y',y
!dy it's a vector of Grad-y (eq. (5)  ref. gradients of
! the energy with respect to the internal coordinates
          do i=1,nInter
            do j=1,nPoints
              dy((i-1)*nPoints+j) = dy_(i,j)
            enddo
          enddo
!
        return
      end subroutine Setup_kriging
!
      end module kriging
