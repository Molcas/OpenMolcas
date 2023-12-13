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

subroutine XDR_Local_Ham(nbas,isize,jsize,imethod,paratyp,dkhorder,xorder,inS,inK,inV,inpVp,inUL,inUS,nbl,ibl,Lmap,DoFullLT,clight)
! Local (Atom/Block) relativistic transformation of Hamiltonian

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nbas, isize, jsize, imethod, paratyp, dkhorder, xorder, nbl, ibl(nbl), Lmap(nbas)
real(kind=wp), intent(in) :: inS(isize), inV(isize), inpVp(isize), clight
real(kind=wp), intent(inout) :: inK(isize)
real(kind=wp), intent(out) :: inUL(jsize), inUS(jsize)
logical(kind=iwp), intent(in) :: doFullLT
integer(kind=iwp) :: i, j, k, iblock, mbl, ks, kL
real(kind=wp), allocatable :: sK(:,:), sS(:,:), sV(:,:), spVp(:,:), sH(:,:), sB(:,:), sKL(:,:), sSL(:,:), sVL(:,:), spVpL(:,:), &
                              sULL(:,:), sUSL(:,:), tmp(:,:)

! Convert triangle matrices to square matrices

call mma_allocate(sK,nbas,nbas,label='skin')
call mma_allocate(sS,nbas,nbas,label='sSS')
call mma_allocate(sV,nbas,nbas,label='sV')
call mma_allocate(spVp,nbas,nbas,label='spVp')
call mma_allocate(sH,nbas,nbas,label='sHam')
call mma_allocate(sB,nbas,nbas,label='sSav')
k = 0
do i=1,nbas
  do j=1,i
    k = k+1
    sK(j,i) = inK(k)
    sS(j,i) = inS(k)
    sV(j,i) = inV(k)
    spVp(j,i) = inpVp(k)
    if (i /= j) then
      sK(j,i) = inK(k)
      sS(j,i) = inS(k)
      sV(j,i) = inV(k)
      spVp(j,i) = inpVp(k)
    end if
  end do
end do
inUL(:) = Zero
inUS(:) = Zero
if (.not. DoFullLT) then
  sH(:,:) = sK(:,:)+sV(:,:)
end if

! Cycle for each local blocks

ks = 0
do iblock=1,nbl
  mbl = ibl(iblock)
  call mma_allocate(sKL,mbl,mbl,label='skinL')
  call mma_allocate(sSL,mbl,mbl,label='sSSL')
  call mma_allocate(sVL,mbl,mbl,label='sVL')
  call mma_allocate(spVpL,mbl,mbl,label='spVpL')
  call mma_allocate(sULL,mbl,mbl,label='ULlco')
  call mma_allocate(sUSL,mbl,mbl,label='USlco')

  ! Copy block matrices

  do i=1,mbl
    do j=1,mbl
      sKL(j,i) = sK(Lmap(j+ks),Lmap(i+ks))
      sSL(j,i) = sS(Lmap(j+ks),Lmap(i+ks))
      sVL(j,i) = sV(Lmap(j+ks),Lmap(i+ks))
      spVpL(j,i) = spVp(Lmap(j+ks),Lmap(i+ks))
    end do
  end do

  ! Calculate relativistic one-electron Hamiltonian for each blocks

  if (imethod == 2) then

    ! Call X2C driver

    call x2c_ts1e(mbl,sSL,sKL,sVL,spVpL,sULL,sUSL,clight)
  else if (imethod == 3) then

    ! Call BSS driver

    call bss_ts1e(mbl,sSL,sKL,sVL,spVpL,sULL,sUSL,clight)
  else if (imethod == 1) then

    ! Call arbitrary order DKH driver

    call dkh_ts1e(mbl,sSL,sKL,sVL,spVpL,sULL,sUSL,clight,dkhorder,xorder,paratyp)
  end if

  ! Copy back to full matrix

  do i=1,mbl
    do j=1,mbl
      kL = Lmap(j+ks)-1+(Lmap(i+ks)-1)*nbas
      inUL(kL+1) = sULL(j,i)
      inUS(kL+1) = sUSL(j,i)
      if (.not. DoFullLT) then
        sH(Lmap(j+ks),Lmap(i+ks)) = sVL(j,i)
      else
        sB(Lmap(j+ks),Lmap(i+ks)) = sVL(j,i)
      end if
    end do
  end do

  ! End cycle for blocks

  call mma_deallocate(sKL)
  call mma_deallocate(sSL)
  call mma_deallocate(sVL)
  call mma_deallocate(spVpL)
  call mma_deallocate(sULL)
  call mma_deallocate(sUSL)
  ks = ks+mbl
end do

! Apply transformation construct from each blocks

if (DoFullLT) then
  call mma_allocate(tmp,nbas,nbas,label='Tempm')
  call dmxma(nbas,'C','N',inUS,sK,sS,Two*clight)
  call dmxma(nbas,'N','N',sS,inUS,sH,-Two*clight)
  call dmxma(nbas,'N','N',sS,inUL,tmp,One)
  call daxpy_(nbas*nbas,One,tmp,1,sH,1)
  call dmxma(nbas,'C','N',inUL,sK,sS,One)
  call dmxma(nbas,'N','N',sS,inUS,tmp,Two*clight)
  call daxpy_(nbas*nbas,One,tmp,1,sH,1)

  call dmxma(nbas,'C','N',inUL,sV,sS,One)
  call dmxma(nbas,'N','N',sS,inUL,tmp,One)
  call daxpy_(nbas*nbas,One,tmp,1,sH,1)
  call dmxma(nbas,'C','N',inUS,spVp,sS,One)
  call dmxma(nbas,'N','N',sS,inUS,tmp,One)
  call daxpy_(nbas*nbas,One,tmp,1,sH,1)
  call mma_deallocate(tmp)
  ks = 0
  do iblock=1,nbl
    mbl = ibl(iblock)
    do i=1,mbl
      do j=1,mbl
        kL = Lmap(j+ks)-1+(Lmap(i+ks)-1)*nbas
        sH(Lmap(j+ks),Lmap(i+ks)) = sB(Lmap(j+ks),Lmap(i+ks))
      end do
    end do
    ks = ks+mbl
  end do
end if
sV(:,:) = sH(:,:)

! Copy relativistic one-electron Hamiltonian back to inK

k = 1
do i=1,nbas
  do j=1,i
    inK(k) = sV(j,i)
    k = k+1
  end do
end do

! Free temp memories

call mma_deallocate(sK)
call mma_deallocate(sS)
call mma_deallocate(sV)
call mma_deallocate(spVp)
call mma_deallocate(sH)
call mma_deallocate(sB)

return

end subroutine XDR_Local_Ham
