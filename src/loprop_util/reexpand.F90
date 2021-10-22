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

subroutine ReExpand(rMP,nij,nElem,A,B,ij,lMax)

use Index_Functions, only: nTri3_Elem
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nij, nElem, ij, lMax
real(kind=wp), intent(inout) :: rMP(nij,nElem)
real(kind=wp), intent(in) :: A(3), B(3)
integer(kind=iwp) :: iElem, ix, iy, iz, jElem, jx, jy, jz, k, l
real(kind=wp) :: ABx, ABx_, ABy, ABy_, ABz, ABz_, temp
#include "itmax.fh"
#include "binom.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
!call RecPrt('A',' ',A,1,3)
!call RecPrt('B',' ',B,1,3)
!call RecPrt('rMP',' ',rMP,nij,nElem)
do l=lMax,0,-1
  iElem = nTri3_Elem(l)
  do ix=l,0,-1
    ABx = A(1)-B(1)
    do iy=l-ix,0,-1
      ABy = A(2)-B(2)
      iz = l-ix-iy
      ABz = A(3)-B(3)
      iElem = iElem+1
      !write(u6,*)
      !write(u6,*)
      !write(u6,*) 'ix,iy,iz=',ix,iy,iz
      !write(u6,*)

      temp = Zero
      do jx=0,ix
        do jy=0,iy
          do jz=0,iz
            !write(u6,*) 'jx,jy,jz=',jx,jy,jz

            if (ix-jx == 0) then
              ABx_ = One
            else
              ABx_ = ABx**(ix-jx)
            end if
            if (iy-jy == 0) then
              ABy_ = One
            else
              ABy_ = ABy**(iy-jy)
            end if
            if (iz-jz == 0) then
              ABz_ = One
            else
              ABz_ = ABz**(iz-jz)
            end if
            k = jx+jy+jz
            jElem = nTri3_Elem(k)+(jy+jz)*(jy+jz+1)/2+jz+1
            !write(u6,*) 'jElem=',jElem
            !write(u6,*) binom(ix,jx),binom(iy,jy),binom(iz,jz),rMP(ij,jElem),ABx_,ABy_,ABz_
            temp = temp+binom(ix,jx)*binom(iy,jy)*binom(iz,jz)*rMP(ij,jElem)*ABx_*ABy_*ABz_

          end do
        end do
      end do
      rMP(ij,iElem) = temp

    end do
  end do
end do
!call RecPrt('rMP',' ',rMP,nij,nElem)

return

end subroutine ReExpand
