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

subroutine Shell_MxDens(Dens,DMax,nSkal)
!***********************************************************************
!                                                                      *
!  Subroutine Shell_MxDens:   returns max density values for each      *
!                             shell pair...                            *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Symmetry_Info, only: nIrrep
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSkal
real(kind=wp), intent(in) :: dens(*)
real(kind=wp), intent(out) :: dmax(nskal,nskal)
integer(kind=iwp) :: i, ia, ie, ij, ijOff, irp, iSh, j, ja, je, jSh, m, n
integer(kind=iwp), external :: nbfshl

ijoff = 0
dmax(:,:) = Zero
do irp=0,nirrep-1
  ie = 0
  do iSh=1,nSkal
    n = nbfshl(ish,irp)
    ia = ie+1
    ie = ie+n
    je = 0
    do jSh=1,iSh
      m = nbfshl(jsh,irp)
      ja = je+1
      je = je+m
      do i=ia,ie
        ij = nTri_Elem(i-1)+ja+ijoff
        do j=ja,min(i,je)
          dmax(jsh,ish) = max(dmax(jsh,ish),abs(dens(ij)))
          ij = ij+1
        end do
      end do
      dmax(ish,jsh) = dmax(jsh,ish)
    end do
  end do
  ijoff = ijoff+nTri_Elem(ie)
end do

end subroutine Shell_MxDens
