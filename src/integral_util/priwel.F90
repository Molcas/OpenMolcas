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

subroutine priwel(k,alfa,beta,r0,a,gri,nz,isum,grin)

use welcom, only: iPot3, kMax
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: k, nz, iSum
real(kind=wp), intent(in) :: alfa(nz), Beta, r0, a(nz)
real(kind=wp), intent(inout) :: gri(nz,isum)
real(kind=wp), intent(out) :: grin(nz,0:k,k/2+1,k/4+1)
integer(kind=iwp) :: i, iDiv, indst, iPot3i, iv(kmax), ix, ix2, ixS, ixyS, ixyz, iy, iy2, iyS, iz, j, jj, l

call binte(k,alfa,beta,r0,a,grin,nz)
!call RecPrt(' In PriWel: Grin',' ',Grin,nz,(k+1)*(k/2+1)*(k/4+1))

! distribute the integrals into gri

indst = 1
gri(:,1) = grin(:,0,1,1)
do i=1,k
  ipot3i = ipot3(i)
  gri(:,indst+1:indst+iPot3i) = Zero
  do j=1,ipot3i
    jj = j
    do l=i,1,-1
      idiv = ipot3(l-1)
      ixyz = (jj-1)/idiv+1
      iv(l) = ixyz
      jj = jj-(ixyz-1)*idiv
    end do

    ! the potency vector for this integral is now ready
    ! ... now analyze it

    ix = 0
    iy = 0
    iz = 0
    do l=1,i
      if (iv(l) == 1) ix = ix+1
      if (iv(l) == 2) iy = iy+1
      if (iv(l) == 3) iz = iz+1
    end do
    ix2 = (ix/2)*2
    if (ix2 /= ix) cycle
    iy2 = (iy/2)*2
    if (iy2 /= iy) cycle
    ixs = max(ix,iy)
    iys = min(ix,iy)
    ixys = (ixs+iys)/2+1
    iys = iys/2+1
    gri(:,j+indst) = grin(:,i,ixys,iys)
  end do
  indst = indst+ipot3i
end do
!call RecPrt(' In PriWel:gri',' ',gri,nz,isum)

return

end subroutine priwel
