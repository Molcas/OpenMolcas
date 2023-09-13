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

implicit real*8(a-h,o-w,y-z),integer(x)
dimension bikcof(ndet,ifns)
dimension minspn(0:nel), maxspn(0:nel), nkspn(0:nel)
dimension locca(nel), lnocca(nel)
dimension xspin((nel+1)*(nalf+1))
dimension ialfs(nalf), ibets(nbet), ianti(ifns)
dimension dum(1)

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
  write(6,*) ' Discrepancy in IFNS:',ifns,xspin((nel+1)*(nbet+1))
  call abend_cvb()
end if
call imove_cvb(maxspn,nkspn,nel+1)
call occupy_cvb(nkspn,nel,locca,lnocca)
! Loop:
index = 1
! Determine pairings
2100 continue
do ib=1,nbet
  ibets(ib) = locca(ib)
  do ia=nalf,1,-1
    ialfs(ib) = lnocca(ia)
    if (ialfs(ib) < ibets(ib)) then
      do iachek=1,ib-1
        if (ialfs(iachek) == ialfs(ib)) goto 2300
      end do
      goto 2500
    end if
2300 continue
  end do
2500 continue
end do

ianti(index) = 0
do ib=1,nbet
  if ((mod(ialfs(ib),2) == 1) .and. (ialfs(ib) == ibets(ib)-1)) ianti(index) = ianti(index)-1
end do
call loind_cvb(nel,nbet,nkspn,minspn,maxspn,locca,lnocca,index,xspin,*2100)

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
      if (ianti(l) == k) goto 5200
    end do
    write(6,*) ' Error - swap function not found!',k,ianti(k)
    call abend_cvb()
5200 call dswap_(ndet,bikcof(1,k),1,bikcof(1,l),1)
    ianti(l) = ianti(k)
    ianti(k) = k
  end if
end do

! Orthonormalize:
call schmidtn_cvb(bikcof,ifns,dum,ndet,0)

return

end subroutine serber_cvb
