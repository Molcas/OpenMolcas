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

subroutine serber_cvb(bikcof,nel,nalf,nbet,ndet,ifns,minspn,maxspn,nkspn,locca,lnocca,xspin,ialfs,ibets,ianti)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nel, nalf, nbet, ndet, ifns, minspn(0:nel), maxspn(0:nel), nkspn(0:nel), locca(nel), lnocca(nel), &
                     xspin((nel+1)*(nalf+1)), ialfs(nalf), ibets(nbet), ianti(ifns)
real(kind=wp) :: bikcof(ndet,ifns)
integer(kind=iwp) :: ia, iachek, iantival, ib, indx, iorb, k, l, nc, rc
real(kind=wp) :: dum(1)

! Serber spin functions

! For each Rumer spin function we determine (minus) the "antisymmetry"
! number, as explained in [Spin Eigenfunctions, R. Pauncz, Sec. 5.5].

! Spin function weight arrays:
do iorb=0,nel
  minspn(iorb) = max(iorb-nalf,0)
  maxspn(iorb) = min(iorb/2,nbet)
end do
call weight_cvb(xspin,minspn,maxspn,nbet,nel)
if (ifns /= xspin((nel+1)*(nbet+1))) then
  write(u6,*) ' Discrepancy in IFNS:',ifns,xspin((nel+1)*(nbet+1))
  call abend_cvb()
end if
call imove_cvb(maxspn,nkspn,nel+1)
call occupy_cvb(nkspn,nel,locca,lnocca)
! Loop:
indx = 1
! Determine pairings
do
  do ib=1,nbet
    ibets(ib) = locca(ib)
    do ia=nalf,1,-1
      ialfs(ib) = lnocca(ia)
      if (ialfs(ib) < ibets(ib)) then
        do iachek=1,ib-1
          if (ialfs(iachek) == ialfs(ib)) exit
        end do
        if (iachek >= ib) exit
      end if
    end do
  end do

  ianti(indx) = 0
  do ib=1,nbet
    if ((mod(ialfs(ib),2) == 1) .and. (ialfs(ib) == ibets(ib)-1)) ianti(indx) = ianti(indx)-1
  end do
  call loind_cvb(nel,nbet,nkspn,minspn,maxspn,locca,lnocca,indx,xspin,rc)
  if (rc == 0) exit
end do

! Sort according to decreasing values of IANTI:
nc = 0
do iantival=nbet,0,-1
  do k=1,ifns
    if (ianti(k) == -iantival) then
      nc = nc+1
      ianti(k) = nc
    end if
  end do
end do

do k=1,ifns
  if (ianti(k) /= k) then
    do l=1,ifns
      if (ianti(l) == k) exit
    end do
    if (l > ifns) then
      write(u6,*) ' Error - swap function not found!',k,ianti(k)
      call abend_cvb()
    end if
    call dswap_(ndet,bikcof(1,k),1,bikcof(1,l),1)
    ianti(l) = ianti(k)
    ianti(k) = k
  end if
end do

! Orthonormalize:
call schmidtn_cvb(bikcof,ifns,dum,ndet,0)

return

end subroutine serber_cvb
