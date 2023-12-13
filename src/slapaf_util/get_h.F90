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
! Copyright (C) 2010, Roland Lindh                                     *
!***********************************************************************

subroutine Get_H(F,nX)
!***********************************************************************
!                                                                      *
!     Object: to get the force constant matrix in cartesians           *
!                                                                      *
!                                                                      *
!     Author: Roland Lindh                                             *
!             Department of Quantum Chemistry                          *
!             Uppsala University, Sweden                               *
!             October 2010                                             *
!***********************************************************************

use Slapaf_Info, only: Coor, mTROld, nDimBC, Smmtrc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nX
real(kind=wp), intent(out) :: F(nX**2)
integer(kind=iwp) :: i, ii, iix, ij, ik, ix, j, jj, jjx, jx, k, nAtom, nBMx, nDoF, nHss, nInter
real(kind=wp) :: tmp_ij
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: BMtrx(:,:), H(:,:), Tmp2(:)

!#define _DEBUGPRINT_
nDoF = nDimBC
nInter = nDimBC-mTROld
nAtom = size(Coor,2)

call mma_allocate(Tmp2,nX**2,Label='Tmp2')
call mma_allocate(H,nInter,nInter,Label='H')
! If there is an updated Hessian in the runfile, use it;
! otherwise use the last computed one.
! (Hss_upd must be removed every time Hss_Q is added)
call Qpg_dArray('Hss_upd',Found,nHss)
if (Found .and. (nHss == nInter**2)) then
  call Get_dArray('Hss_upd',H,nInter**2)
else
  call Get_dArray('Hss_Q',H,nInter**2)
end if
call mma_allocate(BMtrx,nX,nInter,Label='BMtrx')
! If there is an old BMx stored, use it;
! otherwise use the current BMx
call Qpg_dArray('BMxOld',Found,nBMx)
if (Found .and. (nBMx == nX*nInter)) then
  call Get_dArray('BMxOld',BMtrx,nX*nInter)
else
  call Get_dArray('BMtrx',BMtrx,nX*nInter)
end if

!                                                                      *
!***********************************************************************
!                                                                      *
! The Hessian matrix (zero elements over translations and
! rotations are explicitly excluded).

#ifdef _DEBUGPRINT_
call RecPrt('BMtrx',' ',BMtrx,nX,nInter)
#endif
ii = 0
do i=1,nAtom
  do ix=1,3

    if (Smmtrc(ix,i)) then
      iix = (i-1)*3+ix
      ii = ii+1
      do j=1,nInter
        ij = (j-1)*nDoF+ii

        tmp_ij = Zero
        do k=1,nInter
          tmp_ij = tmp_ij+BMtrx(iix,k)*H(k,j)
        end do
        Tmp2(ij) = tmp_ij

      end do
    end if

  end do
end do

do i=1,nDoF

  jj = 0
  do j=1,nAtom
    do jx=1,3

      if (Smmtrc(jx,j)) then
        jjx = (j-1)*3+jx
        jj = jj+1
        ij = (jj-1)*nDoF+i
        tmp_ij = Zero
        do k=1,nInter
          ik = (k-1)*nDoF+i
          tmp_ij = tmp_ij+Tmp2(ik)*BMtrx(jjx,k)
        end do
        F(ij) = tmp_ij
      end if

    end do
  end do

end do

#ifdef _DEBUGPRINT_
call RecPrt('Hessian (cartesian)',' ',F,nDoF,nDoF)
#endif
call Put_dArray('FC-Matrix',F,nDoF**2)

call mma_deallocate(BMtrx)
call mma_deallocate(H)
call mma_deallocate(Tmp2)

!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Get_H
