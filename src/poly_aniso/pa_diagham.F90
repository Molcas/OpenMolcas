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

subroutine pa_diagham(exch,npair,i_pair,nneq,neq,nexch,nmax,lmax,eso,HLIN1,HLIN3,HLIN9,HDIP,HKEX,HDMO,HITO,Dipol,DM_exchange, &
                      AnisoLines1,AnisoLines3,AnisoLines9,KE,JITO_exchange,WLIN1,WLIN3,WLIN9,WLIN,WDIP,WKEX,WDMO,WITO,W,Z)
! this function builds and diagonalizes the interaction Hamiltonians

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: exch, npair, i_pair(npair,2), nneq, neq(nneq), nexch(nneq), nmax, lmax
real(kind=wp), intent(in) :: eso(nneq,nmax)
complex(kind=wp), intent(in) :: HLIN1(npair,nmax,nmax,nmax,nmax), HLIN3(npair,nmax,nmax,nmax,nmax), &
                                HLIN9(npair,nmax,nmax,nmax,nmax), HDIP(npair,nmax,nmax,nmax,nmax), &
                                HKEX(npair,nmax,nmax,nmax,nmax), HDMO(npair,nmax,nmax,nmax,nmax), HITO(npair,nmax,nmax,nmax,nmax)
logical(kind=iwp), intent(in) :: Dipol, DM_exchange, AnisoLines1, AnisoLines3, AnisoLines9, KE, JITO_exchange
real(kind=wp), intent(out) :: wlin1(exch), wlin3(exch), wlin9(exch), wlin(exch), wdip(exch), wkex(exch), wdmo(exch), wito(exch), &
                              w(exch)
complex(kind=wp), intent(out) :: Z(exch,exch)
integer(kind=iwp) :: i, i1, i2, info, is1, is2, j, js1, js2, l, lb, lb1, lb2, lp, lwork, nb, nb1, nb2
integer(kind=iwp), allocatable :: ibas(:,:), icoord(:), intc(:), nind(:,:)
real(kind=wp), allocatable :: rwork(:)
complex(kind=wp), allocatable :: HTOT(:,:), work(:)
integer(kind=iwp), external :: norder

! allocate memory and initialize variables:
if (exch >= 0) then
  call mma_allocate(HTOT,exch,exch,'HTOT')
  call mma_allocate(WORK,(2*exch-1),'WORK')
  if (lmax >= 0) call mma_allocate(ibas,exch,lmax,'ibas')
  call mma_allocate(rwork,(3*exch-2),'rwork')
end if

if (lmax >= 0) then
  call mma_allocate(nind,lmax,2,'nind')
  call mma_allocate(intc,lmax,'intc')
  call mma_allocate(icoord,lmax,'icoord')
end if
call zcopy_(exch*exch,[cZero],0,Z,1)
call dcopy_(exch,[Zero],0,w,1)
call dcopy_(exch,[Zero],0,wlin,1)
call dcopy_(exch,[Zero],0,wlin1,1)
call dcopy_(exch,[Zero],0,wlin3,1)
call dcopy_(exch,[Zero],0,wlin9,1)
call dcopy_(exch,[Zero],0,wdip,1)
call dcopy_(exch,[Zero],0,wkex,1)
call dcopy_(exch,[Zero],0,wdmo,1)
call dcopy_(exch,[Zero],0,wito,1)

! generate the tables:
l = 0
call icopy(2*lmax,[0],0,nind,1)
do i=1,nneq
  do j=1,neq(i)
    l = l+1
    nind(l,1) = i
    nind(l,2) = j
  end do
end do

call icopy(lmax,[0],0,intc,1)
intc(1) = 1
if (lmax > 1) then
  do i=2,lmax
    i1 = nind(i-1,1)
    intc(i) = intc(i-1)*nexch(i1)
  end do
end if

call icopy(exch*lmax,[0],0,ibas,1)
do nb=1,exch
  nb1 = nb-1
  do i=1,lmax
    ibas(nb,lmax-i+1) = nb1/intc(lmax-i+1)
    nb1 = nb1-ibas(nb,lmax-i+1)*intc(lmax-i+1)
  end do
end do
! build the interaction Hamiltonians
!----------------------------------------------------------------------!
if (AnisoLines1) then
  call zcopy_(exch*exch,[cZero],0,HTOT,1)
  do nb1=1,exch
    do lp=1,npair
      call icopy(lmax,[0],0,icoord,1)
      do i=1,lmax
        icoord(i) = ibas(nb1,i)
      end do

      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1)
      i2 = nind(lb2,1)
      is1 = icoord(lb1)+1
      is2 = icoord(lb2)+1
      do js1=1,nexch(i1)
        icoord(lb1) = js1-1
        do js2=1,nexch(i2)
          icoord(lb2) = js2-1
          nb2 = norder(icoord,intc,lmax)
          HTOT(nb1,nb2) = HTOT(nb1,nb2)+HLIN1(lp,is1,js1,is2,js2)
        end do ! js2
      end do ! js1
    end do ! lp

    l = 0
    lb = 1
    do i=1,nneq
      do j=1,neq(i)
        l = l+1
        if (l == lb) then
          HTOT(nb1,nb1) = HTOT(nb1,nb1)+eso(i,ibas(nb1,lb)+1)*cOne
          if ((lb+1) <= (lmax)) lb = lb+1
        end if
      end do !j
    end do !i
  end do !nb1
  ! diagonalize
  info = 0
  lwork = 0
  lwork = 2*exch-1
  call dcopy_((3*exch-2),[Zero],0,rwork,1)
  call zcopy_((2*exch-1),[cZero],0,WORK,1)
  call zheev('n','u',exch,htot,exch,wlin1,work,lwork,rwork,info)
end if

!----------------------------------------------------------------------!
if (AnisoLines3) then
  call zcopy_(exch*exch,[cZero],0,HTOT,1)
  do nb1=1,exch
    do lp=1,npair
      call icopy(lmax,[0],0,icoord,1)
      do i=1,lmax
        icoord(i) = ibas(nb1,i)
      end do

      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1)
      i2 = nind(lb2,1)
      is1 = icoord(lb1)+1
      is2 = icoord(lb2)+1
      do js1=1,nexch(i1)
        icoord(lb1) = js1-1
        do js2=1,nexch(i2)
          icoord(lb2) = js2-1
          nb2 = norder(icoord,intc,lmax)
          HTOT(nb1,nb2) = HTOT(nb1,nb2)+HLIN3(lp,is1,js1,is2,js2)
        end do ! js2
      end do ! js1
    end do ! lp

    l = 0
    lb = 1
    do i=1,nneq
      do j=1,neq(i)
        l = l+1
        if (l == lb) then
          HTOT(nb1,nb1) = HTOT(nb1,nb1)+eso(i,ibas(nb1,lb)+1)*cOne
          if ((lb+1) <= (lmax)) lb = lb+1
        end if
      end do !j
    end do !i
  end do !nb1
  ! diagonalize
  info = 0
  lwork = 0
  lwork = 2*exch-1
  call dcopy_((3*exch-2),[Zero],0,rwork,1)
  call zcopy_((2*exch-1),[cZero],0,WORK,1)
  call zheev('n','u',exch,htot,exch,wlin3,work,lwork,rwork,info)
end if

!----------------------------------------------------------------------!
if (AnisoLines9) then
  call zcopy_(exch*exch,[cZero],0,HTOT,1)
  do nb1=1,exch
    do lp=1,npair
      call icopy(lmax,[0],0,icoord,1)
      do i=1,lmax
        icoord(i) = ibas(nb1,i)
      end do

      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1)
      i2 = nind(lb2,1)
      is1 = icoord(lb1)+1
      is2 = icoord(lb2)+1
      do js1=1,nexch(i1)
        icoord(lb1) = js1-1
        do js2=1,nexch(i2)
          icoord(lb2) = js2-1
          nb2 = norder(icoord,intc,lmax)
          HTOT(nb1,nb2) = HTOT(nb1,nb2)+HLIN9(lp,is1,js1,is2,js2)
        end do ! js2
      end do ! js1
    end do ! lp

    l = 0
    lb = 1
    do i=1,nneq
      do j=1,neq(i)
        l = l+1
        if (l == lb) then
          HTOT(nb1,nb1) = HTOT(nb1,nb1)+eso(i,ibas(nb1,lb)+1)*cOne
          if ((lb+1) <= (lmax)) lb = lb+1
        end if
      end do !j
    end do !i
  end do !nb1
  ! diagonalize
  info = 0
  lwork = 0
  lwork = 2*exch-1
  call dcopy_((3*exch-2),[Zero],0,rwork,1)
  call zcopy_((2*exch-1),[cZero],0,WORK,1)
  call zheev('n','u',exch,htot,exch,wlin9,work,lwork,rwork,info)
end if

!----------------------------------------------------------------------!
if (AnisoLines1 .or. AnisoLines3 .or. AnisoLines9) then
  call zcopy_(exch*exch,[cZero],0,HTOT,1)
  do nb1=1,exch
    do lp=1,npair
      call icopy(lmax,[0],0,icoord,1)
      do i=1,lmax
        icoord(i) = ibas(nb1,i)
      end do

      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1)
      i2 = nind(lb2,1)
      is1 = icoord(lb1)+1
      is2 = icoord(lb2)+1
      do js1=1,nexch(i1)
        icoord(lb1) = js1-1
        do js2=1,nexch(i2)
          icoord(lb2) = js2-1
          nb2 = norder(icoord,intc,lmax)
          HTOT(nb1,nb2) = HTOT(nb1,nb2)+HLIN1(lp,is1,js1,is2,js2)+HLIN3(lp,is1,js1,is2,js2)+HLIN9(lp,is1,js1,is2,js2)
        end do ! js2
      end do ! js1
    end do ! lp

    l = 0
    lb = 1
    do i=1,nneq
      do j=1,neq(i)
        l = l+1
        if (l == lb) then
          HTOT(nb1,nb1) = HTOT(nb1,nb1)+eso(i,ibas(nb1,lb)+1)*cOne
          if ((lb+1) <= (lmax)) lb = lb+1
        end if
      end do !j
    end do !i
  end do !nb1
  ! diagonalize
  info = 0
  lwork = 0
  lwork = 2*exch-1
  call dcopy_((3*exch-2),[Zero],0,rwork,1)
  call zcopy_((2*exch-1),[cZero],0,WORK,1)
  call zheev('n','u',exch,htot,exch,wlin,work,lwork,rwork,info)
end if

!----------------------------------------------------------------------!
if (Dipol) then
  call zcopy_(exch*exch,[cZero],0,HTOT,1)
  do nb1=1,exch
    do lp=1,npair
      call icopy(lmax,[0],0,icoord,1)
      do i=1,lmax
        icoord(i) = ibas(nb1,i)
      end do

      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1)
      i2 = nind(lb2,1)
      is1 = icoord(lb1)+1
      is2 = icoord(lb2)+1
      do js1=1,nexch(i1)
        icoord(lb1) = js1-1
        do js2=1,nexch(i2)
          icoord(lb2) = js2-1
          nb2 = norder(icoord,intc,lmax)
          HTOT(nb1,nb2) = HTOT(nb1,nb2)+HDIP(lp,is1,js1,is2,js2)
        end do ! js2
      end do ! js1
    end do ! lp
    l = 0
    lb = 1
    do i=1,nneq
      do j=1,neq(i)
        l = l+1
        if (l == lb) then
          ! kind=8, complex double precision
          HTOT(nb1,nb1) = HTOT(nb1,nb1)+eso(i,ibas(nb1,lb)+1)*cOne
          if ((lb+1) <= (lmax)) lb = lb+1
        end if
      end do !j
    end do !i
  end do !nb1
  ! diagonalize
  info = 0
  lwork = 0
  lwork = 2*exch-1
  call dcopy_((3*exch-2),[Zero],0,rwork,1)
  call zcopy_((2*exch-1),[cZero],0,WORK,1)
  call zheev('n','u',exch,htot,exch,wdip,work,lwork,rwork,info)
end if

!----------------------------------------------------------------------!
if (DM_exchange) then
  call zcopy_(exch*exch,[cZero],0,HTOT,1)
  do nb1=1,exch
    do lp=1,npair
      call icopy(lmax,[0],0,icoord,1)
      do i=1,lmax
        icoord(i) = ibas(nb1,i)
      end do

      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1)
      i2 = nind(lb2,1)
      is1 = icoord(lb1)+1
      is2 = icoord(lb2)+1
      do js1=1,nexch(i1)
        icoord(lb1) = js1-1
        do js2=1,nexch(i2)
          icoord(lb2) = js2-1
          nb2 = norder(icoord,intc,lmax)
          HTOT(nb1,nb2) = HTOT(nb1,nb2)+HDMO(lp,is1,js1,is2,js2)
        end do ! js2
      end do ! js1
    end do ! lp

    l = 0
    lb = 1
    do i=1,nneq
      do j=1,neq(i)
        l = l+1
        if (l == lb) then
          HTOT(nb1,nb1) = HTOT(nb1,nb1)+eso(i,ibas(nb1,lb)+1)*cOne
          if ((lb+1) <= (lmax)) lb = lb+1
        end if
      end do !j
    end do !i
  end do !nb1
  ! diagonalize
  info = 0
  lwork = 0
  lwork = 2*exch-1
  call dcopy_((3*exch-2),[Zero],0,rwork,1)
  call zcopy_((2*exch-1),[cZero],0,WORK,1)
  call zheev('n','u',exch,htot,exch,wdmo,work,lwork,rwork,info)
end if

!----------------------------------------------------------------------!
if (KE) then
  call zcopy_(exch*exch,[cZero],0,HTOT,1)
  do nb1=1,exch
    do lp=1,npair
      call icopy(lmax,[0],0,icoord,1)
      do i=1,lmax
        icoord(i) = ibas(nb1,i)
      end do

      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1)
      i2 = nind(lb2,1)
      is1 = icoord(lb1)+1
      is2 = icoord(lb2)+1
      do js1=1,nexch(i1)
        icoord(lb1) = js1-1
        do js2=1,nexch(i2)
          icoord(lb2) = js2-1
          nb2 = norder(icoord,intc,lmax)
          HTOT(nb1,nb2) = HTOT(nb1,nb2)+HKEX(lp,is1,js1,is2,js2)
        end do ! js2
      end do ! js1
    end do ! lp
  end do !nb1
  ! diagonalize
  info = 0
  lwork = 0
  lwork = 2*exch-1
  call dcopy_((3*exch-2),[Zero],0,rwork,1)
  call zcopy_((2*exch-1),[cZero],0,WORK,1)
  call zheev('n','u',exch,htot,exch,wkex,work,lwork,rwork,info)
end if

!----------------------------------------------------------------------!
if (JITO_exchange) then
  call zcopy_(exch*exch,[cZero],0,HTOT,1)
  do nb1=1,exch
    do lp=1,npair
      call icopy(lmax,[0],0,icoord,1)
      do i=1,lmax
        icoord(i) = ibas(nb1,i)
      end do

      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1)
      i2 = nind(lb2,1)
      is1 = icoord(lb1)+1
      is2 = icoord(lb2)+1
      do js1=1,nexch(i1)
        icoord(lb1) = js1-1
        do js2=1,nexch(i2)
          icoord(lb2) = js2-1
          nb2 = norder(icoord,intc,lmax)
          HTOT(nb1,nb2) = HTOT(nb1,nb2)+HITO(lp,is1,js1,is2,js2)
        end do ! js2
      end do ! js1
    end do ! lp
  end do !nb1
  ! diagonalize
  info = 0
  lwork = 0
  lwork = 2*exch-1
  call dcopy_((3*exch-2),[Zero],0,rwork,1)
  call zcopy_((2*exch-1),[cZero],0,WORK,1)
  call zheev('n','u',exch,htot,exch,wito,work,lwork,rwork,info)
end if

!----------------------------------------------------------------------!
!cccccccc  total Hamiltonian cccccccc
call zcopy_(exch*exch,[cZero],0,HTOT,1)
do nb1=1,exch
  do lp=1,npair
    call icopy(lmax,[0],0,icoord,1)
    do i=1,lmax
      icoord(i) = ibas(nb1,i)
    end do

    lb1 = i_pair(lp,1)
    lb2 = i_pair(lp,2)
    i1 = nind(lb1,1)
    i2 = nind(lb2,1)
    is1 = icoord(lb1)+1
    is2 = icoord(lb2)+1
    do js1=1,nexch(i1)
      icoord(lb1) = js1-1
      do js2=1,nexch(i2)
        icoord(lb2) = js2-1
        nb2 = norder(icoord,intc,lmax)

        HTOT(nb1,nb2) = HTOT(nb1,nb2)+HLIN1(lp,is1,js1,is2,js2)+HLIN3(lp,is1,js1,is2,js2)+HLIN9(lp,is1,js1,is2,js2)+ &
                        HDIP(lp,is1,js1,is2,js2)+HKEX(lp,is1,js1,is2,js2)+HITO(lp,is1,js1,is2,js2)
      end do ! js2
    end do ! js1
  end do ! lp

  if (.not. KE) then
    l = 0
    lb = 1
    do i=1,nneq
      do j=1,neq(i)
        l = l+1
        if (l == lb) then
          !kind=8, complex double precision
          HTOT(nb1,nb1) = HTOT(nb1,nb1)+eso(i,ibas(nb1,lb)+1)*cOne
          if ((lb+1) <= (lmax)) lb = lb+1
        end if
      end do !j
    end do !i
  end if ! .not. KE
end do !nb1
! diagonalize
info = 0
lwork = 0
lwork = 2*exch-1
call dcopy_((3*exch-2),[Zero],0,rwork,1)
call zcopy_((2*exch-1),[cZero],0,WORK,1)
call zheev('v','u',exch,htot,exch,w,work,lwork,rwork,info)
if (info == 0) then
  call zcopy_(exch*exch,htot,1,Z,1)
else
  write(u6,'(A,i10)') 'DIAG:  non-zero Return: INFO=',info
end if

!----------------------------------------------------------------------!
! deallocate memory:
if (exch >= 0) then
  call mma_deallocate(HTOT)
  call mma_deallocate(WORK)
  if (lmax >= 0) call mma_deallocate(ibas)
  call mma_deallocate(rwork)
end if

if (lmax >= 0) then
  call mma_deallocate(nind)
  call mma_deallocate(intc)
  call mma_deallocate(icoord)
end if

return

end subroutine pa_diagham
