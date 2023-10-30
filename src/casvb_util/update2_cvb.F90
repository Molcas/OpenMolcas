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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine update2_cvb(orbs,cvb,orbsp,cvbp,sorbs,dxorg,ic,norb,nvb,nprorb,npr,orbopt,strucopt,sym,iorts,nort)

use casvb_global, only: ipr, wdx
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ic, norb, nvb, nprorb, npr, nort, iorts(2,nort)
real(kind=wp), intent(out) :: orbs(norb,norb), cvb(nvb), sorbs(norb,norb)
real(kind=wp), intent(in) :: orbsp(norb,norb), cvbp(nvb), dxorg(npr)
logical(kind=iwp), intent(in) :: orbopt, strucopt, sym
integer(kind=iwp) :: ierr, ij, iorb, iort, j, jorb, k, korb, l, lorb
real(kind=wp) :: dum(1), fac, sdidj
real(kind=wp), allocatable :: sorbsinv(:,:)

call free2all_cvb(dxorg,wdx,1)
if ((ipr(3) >= 3) .and. (ic == 1)) then
  write(u6,'(/,a)') ' Update vector :'
  call vecprint_cvb(wdx,npr)
end if

orbs(:,:) = orbsp(:,:)
cvb(:) = cvbp(:)

if (orbopt) then
  call mxattb_cvb(orbsp,orbsp,norb,norb,norb,sorbs)

  ij = 0
  do iorb=1,norb
    do jorb=1,norb
      if (jorb /= iorb) then
        ij = ij+1
        orbs(:,iorb) = orbs(:,iorb)+wdx(ij)*orbsp(:,jorb)
      end if
    end do
  end do

  ! 2nd-order correction for orthogonality constraints:
  call mma_allocate(sorbsinv,norb,norb,label='sorbsinv')
  sorbsinv(:,:) = sorbs(:,:)
  call mxinv_cvb(sorbsinv,norb)
  do iort=1,nort
    iorb = iorts(1,iort)
    jorb = iorts(2,iort)
    sdidj = zero
    do k=1,norb-1
      korb = k
      if (korb >= iorb) korb = korb+1
      do l=1,norb-1
        lorb = l
        if (lorb >= jorb) lorb = lorb+1
        sdidj = sdidj+sorbs(korb,lorb)*wdx(k+(iorb-1)*(norb-1))*wdx(l+(jorb-1)*(norb-1))
      end do
    end do
    fac = -Half*sdidj
    do j=1,norb
      orbs(:,iorb) = orbs(:,iorb)+fac*orbsp(:,j)*sorbsinv(j,jorb)
      orbs(:,jorb) = orbs(:,jorb)+fac*orbsp(:,j)*sorbsinv(j,iorb)
    end do
  end do
  call mma_deallocate(sorbsinv)
end if

if (strucopt) then
  cvb(:) = cvb(:)+wdx(nprorb+1:nprorb+nvb)
  call scalstruc_cvb(orbs,cvb)
  call cvbnrm_cvb(cvb)
end if

ierr = 0
call nize_cvb(orbs,norb,dum,norb,0,ierr)
if (sym) call symtriz_cvb(orbs,cvb)

return

end subroutine update2_cvb
