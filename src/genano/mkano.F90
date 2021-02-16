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
! Copyright (C) Per-Olof Widmark                                       *
!***********************************************************************
!***********************************************************************
!                                                                      *
!======================================================================*
!                                                                      *
! Author: Per-Olof Widmark                                             *
!         IBM Sweden                                                   *
!                                                                      *
!***********************************************************************

subroutine MkAno()

use Genano_globals, only: MxLqn, nPrim, Ssym, tDSym, thr, rowise, Title, symlab
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: Lu, i, iLqn, ind, iPtr, j, k, l, n, nBig, nOccNo, nTri, nContr(0:MxLqn)
real(kind=wp) :: ChkSum, p, t
real(kind=wp), allocatable :: s(:), u(:), c(:), d(:), q(:), OccNo(:)
integer(kind=iwp), external :: IsFreeUnit

Lu = IsFreeUnit(17)
call molcas_open(Lu,'ANO')
iPtr = 1
write(u6,*)
call CollapseOutput(1,'   Contraction coefficients')
write(u6,'(3x,a)') '   ------------------------'

call mma_allocate(OccNo,sum(nPrim),label='OccNo')

! Loop over angular quantum number

nOccNo = 0
ChkSum = Zero
do iLqn=0,MxLqn
  n = nPrim(iLqn)
  nContr(iLqn) = 0
  if (n <= 0) cycle
  write(u6,*)
  write(u6,*) '*** Contraction coefficients for the ',symlab(iLqn*(iLqn+1)+1)(3:3),' shell ***'
  write(u6,*)
  nTri = n*(n+1)/2
  call mma_allocate(s,nTri,label='s')
  call mma_allocate(u,n*n,label='u')
  u(:) = Zero
  do i=1,n
    ind = i+(i-1)*n
    u(ind) = One
  end do
  do i=1,nTri
    s(i) = Ssym(iPtr-1+i)
  end do
  !call TriPrt('s matrix','(6f12.6)',s,n)
  !write(u6,*) 'u matrix'
  !call sqprt(u,n)
  call Jacob(s,u,n,n)
  !call TriPrt('s matrix','(6f12.6)',s,n)
  !write(u6,*) 'u matrix'
  !call sqprt(u,n)
  !write(u6,*)
  do i=1,n
    s(i) = sqrt(s(i*(i+1)/2))
  end do
  call mma_allocate(c,n*n,label='c')
  call mma_allocate(q,n*n,label='c')
  do i=1,n
    do j=1,n
      t = Zero
      p = Zero
      do k=1,n
        t = t+u(i+(k-1)*n)/s(k)*u(j+(k-1)*n)
        p = p+u(i+(k-1)*n)*s(k)*u(j+(k-1)*n)
      end do
      c(j+(i-1)*n) = t
      q(j+(i-1)*n) = p
    end do
  end do
  call mma_deallocate(s)
  call mma_deallocate(u)
  call mma_allocate(d,nTri,label='d')
  do i=1,n
    do j=1,i
      t = Zero
      do k=1,n
        do l=1,n
          ind = min(k,l)+max(k,l)*(max(k,l)-1)/2
          p = q(i+(k-1)*n)*tDsym(ind+iPtr-1)*q(j+(l-1)*n)
          !write(u6,'(a,f12.6)') '   p:',p
          t = t+p
        end do
      end do
      d(j+i*(i-1)/2) = t
    end do
  end do
  call mma_deallocate(q)
  !call TriPrt('d matrix','(6f12.6)',d,n)
  !write(u6,*) 'c matrix'
  !call sqprt(c,n)
  call Jacob(d,c,n,n)
  !call TriPrt('d matrix','(6f12.6)',d,n)
  !write(u6,*) 'c matrix'
  !call sqprt(c,n)
  !write(u6,*)
  nBig = 0
  do i=1,n
    d(i) = d(i*(i+1)/2)
    if (d(i) > thr) nBig = nBig+1
  end do
  !nBig = n
  nContr(iLqn) = nBig
  call Sort_genano(d,c,n,n)
  do i=1,nContr(iLqn)
    OccNo(i+nOccNo) = d(i)
  end do
  nOccNo = nOccNo+nContr(iLqn)
  !write(u6,*) 'c matrix'
  !call sqprt(c,n)
  call NOphase(c,n)
  !write(u6,*) 'c matrix'
  !call sqprt(c,n)
  !write(u6,*)
  write(u6,'(a,10(1x,f9.4))') 'occ    ',(d(i),i=1,nBig)
  write(u6,'(a,10(1x,f9.4))') 'lg(occ)',(log10(d(i)+1.0e-12_wp),i=1,nBig)
  write(u6,*)

  do j=1,n
    write(u6,'(7x,10(1x,f9.4))') (c(j+(i-1)*n),i=1,nBig)
  end do
  do i=1,nBig
    do j=1,n
      do k=1,n
        ChkSum = ChkSum+d(i)*c(j+(i-1)*n)*c(k+(i-1)*n)
      end do
    end do
  end do
  !write(Lu,*) 'Coefficients for the ',symlab(iLqn*(iLqn+1)+1),' shell.'
  !write(Lu,*)
  write(Lu,'(2i5)') n,nBig
  if (rowise) then
    if (n*nBig > 0) then
      do i=1,nBig
        write(Lu,'(10(1x,f15.10))') (c(j+(i-1)*n),j=1,n)
      end do
    end if
  else
    if (n*nBig > 0) then
      do j=1,n
        write(Lu,'(10(1x,f15.10))') (c(j+(i-1)*n),i=1,nBig)
      end do
    end if
  end if
  call mma_deallocate(c)
  call mma_deallocate(d)
  !write(Lu,*)
  !write(Lu,*)
  iPtr = iPtr+(2*iLqn+1)*nTri
end do
call CollapseOutput(0,'   Contraction coefficients')
write(u6,*)
!write(u6,'(a,f12.6)') 'Check sum is',ChkSum
call Add_Info('GENANO_CHKSUM',[ChkSum],1,5)
close(Lu)
Lu = IsFreeUnit(18)
call molcas_open(Lu,'FIG')
call FigOpn(Lu)
call FigPrt(Lu,Title,MxLqn,nContr,OccNo)
call FigCls(Lu)
close(Lu)
call mma_deallocate(OccNo)

return

end subroutine MkAno
