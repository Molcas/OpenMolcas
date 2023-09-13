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

subroutine update2_cvb(orbs,cvb,orbsp,cvbp,sorbs,dxorg,ic,norb,nvb,nprorb,npr,orbopt,strucopt,sym,dx,iorts,nort,sorbsinv)

implicit real*8(a-h,o-z)
logical orbopt, strucopt, sym
#include "print_cvb.fh"
dimension orbs(norb,norb), cvb(nvb)
dimension orbsp(norb,norb), cvbp(nvb)
dimension sorbs(norb,norb)
dimension dxorg(npr), dx(npr)
dimension iorts(2,nort)
dimension sorbsinv(norb,norb)
dimension dum(1)
save zero
data zero/0d0/

call free2all_cvb(dxorg,dx,1)
if ((ip(3) >= 3) .and. (ic == 1)) then
  write(6,'(/,a)') ' Update vector :'
  call vecprint_cvb(dx,npr)
end if

call fmove_cvb(orbsp,orbs,norb*norb)
call fmove_cvb(cvbp,cvb,nvb)

if (orbopt) then
  call mxattb_cvb(orbsp,orbsp,norb,norb,norb,sorbs)

  ij = 0
  do iorb=1,norb
    do jorb=1,norb
      if (jorb /= iorb) then
        ij = ij+1
        do i=1,norb
          orbs(i,iorb) = orbs(i,iorb)+dx(ij)*orbsp(i,jorb)
        end do
      end if
    end do
  end do

  ! 2nd-order correction for orthogonality constraints:
  call fmove_cvb(sorbs,sorbsinv,norb*norb)
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
        sdidj = sdidj+sorbs(korb,lorb)*dx(k+(iorb-1)*(norb-1))*dx(l+(jorb-1)*(norb-1))
      end do
    end do
    fac = -.5d0*sdidj
    do i=1,norb
      do j=1,norb
        orbs(i,iorb) = orbs(i,iorb)+fac*orbsp(i,j)*sorbsinv(j,jorb)
        orbs(i,jorb) = orbs(i,jorb)+fac*orbsp(i,j)*sorbsinv(j,iorb)
      end do
    end do
  end do
end if

if (strucopt) then
  call addvec(cvb,cvb,dx(nprorb+1),nvb)
  call scalstruc_cvb(orbs,cvb)
  call cvbnrm_cvb(cvb)
end if

call nize_cvb(orbs,norb,dum,norb,0,0)
if (sym) call symtriz_cvb(orbs,cvb)

return

end subroutine update2_cvb
