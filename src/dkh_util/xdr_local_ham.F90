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

subroutine XDR_Local_Ham(nbas,isize,jsize,imethod,paratyp,dkhorder,xorder,inS,inK,inV,inpVp,inUL,inUS,inDX,nbl,ibl,Lmap,DoFullLT, &
                         clight)
! Local (Atom/Block) relativistic transformation of Hamiltonian

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nbas, isize, jsize, imethod, paratyp, dkhorder, xorder, inDX(nbas), nbl, ibl(nbl), Lmap(nbas)
real(kind=wp), intent(in) :: inS(isize), inV(isize), inpVp(isize), clight
real(kind=wp), intent(inout) :: inK(isize)
real(kind=wp), intent(out) :: inUL(jsize), inUS(jsize)
logical(kind=iwp), intent(in) :: doFullLT
#include "WrkSpc.fh"
integer(kind=iwp) :: nn, i, j, k, iblock, mbl, ks, kL, jS, jK, jV, jpVp, jSL, jKL, jVL, jpVpL, jULL, jUSL, jH, jtmp, jB

! Convert triangle matrices to square matrices

nn = nbas*nbas+4
call getmem('skin ','ALLOC','REAL',jK,nn)
call getmem('sSS  ','ALLOC','REAL',jS,nn)
call getmem('sV   ','ALLOC','REAL',jV,nn)
call getmem('spVp ','ALLOC','REAL',jpVp,nn)
call getmem('sHam ','ALLOC','REAL',jH,nn)
call getmem('sSav ','ALLOC','REAL',jB,nn)
k = 0
do i=1,nbas
  do j=1,i
    k = k+1
    Work(jK+j-1+(i-1)*nbas) = inK(k)
    Work(jS+j-1+(i-1)*nbas) = inS(k)
    Work(jV+j-1+(i-1)*nbas) = inV(k)
    Work(jpVp+j-1+(i-1)*nbas) = inpVp(k)
    if (i /= j) then
      Work(jK+i-1+(j-1)*nbas) = inK(k)
      Work(jS+i-1+(j-1)*nbas) = inS(k)
      Work(jV+i-1+(j-1)*nbas) = inV(k)
      Work(jpVp+i-1+(j-1)*nbas) = inpVp(k)
    end if
  end do
end do
do i=1,jsize
  inUL(i) = Zero
  inUS(i) = Zero
end do
if (.not. DoFullLT) then
  do i=0,nbas*nbas-1
    Work(jH+i) = Work(jK+i)+Work(jV+i)
  end do
end if

! Cycle for each local blocks

ks = 0
do iblock=1,nbl
  mbl = ibl(iblock)
  call getmem('skinL','ALLOC','REAL',jKL,mbl*mbl+4)
  call getmem('sSSL ','ALLOC','REAL',jSL,mbl*mbl+4)
  call getmem('sVL  ','ALLOC','REAL',jVL,mbl*mbl+4)
  call getmem('spVpL','ALLOC','REAL',jpVpL,mbl*mbl+4)
  call getmem('ULlco','ALLOC','REAL',jULL,mbl*mbl+4)
  call getmem('USlco','ALLOC','REAL',jUSL,mbl*mbl+4)

  ! Copy block matrices

  do i=1,mbl
    do j=1,mbl
      k = (j-1)+(i-1)*mbl
      kL = Lmap(j+ks)-1+(Lmap(i+ks)-1)*nbas
      Work(jKL+k) = Work(jK+kL)
      Work(jSL+k) = Work(jS+kL)
      Work(jVL+k) = Work(jV+kL)
      Work(jpVpL+k) = Work(jpVp+kL)
    end do
  end do

  ! Calculate relativistic one-electron Hamiltonian for each blocks

  if (imethod == 2) then

  ! Call X2C driver

    call x2c_ts1e(mbl,Work(jSL),Work(jKL),Work(jVL),Work(jpVpL),Work(jULL),Work(jUSL),clight)
  else if (imethod == 3) then

  ! Call BSS driver

    call bss_ts1e(mbl,Work(jSL),Work(jKL),Work(jVL),Work(jpVpL),Work(jULL),Work(jUSL),clight)
  else if (imethod == 1) then

  ! Call arbitrary order DKH driver

    call dkh_ts1e(mbl,Work(jSL),Work(jKL),Work(jVL),Work(jpVpL),Work(jULL),Work(jUSL),clight,dkhorder,xorder,paratyp)
  end if

  ! Copy back to full matrix

  do i=1,mbl
    do j=1,mbl
      k = (j-1)+(i-1)*mbl
      kL = Lmap(j+ks)-1+(Lmap(i+ks)-1)*nbas
      inUL(kL+1) = Work(jULL+k)
      inUS(kL+1) = Work(jUSL+k)
      if (.not. DoFullLT) then
        Work(jH+kL) = Work(jVL+k)
      else
        Work(jB+kL) = Work(jVL+k)
      end if
    end do
  end do

  ! End cycle for blocks

  call getmem('skinL','FREE','REAL',jKL,mbl*mbl+4)
  call getmem('sSSL ','FREE','REAL',jSL,mbl*mbl+4)
  call getmem('sVL  ','FREE','REAL',jVL,mbl*mbl+4)
  call getmem('spVpL','FREE','REAL',jpVpL,mbl*mbl+4)
  call getmem('ULlco','FREE','REAL',jULL,mbl*mbl+4)
  call getmem('USlco','FREE','REAL',jUSL,mbl*mbl+4)
  ks = ks+mbl
end do

! Apply transformation construct from each blocks

if (DoFullLT) then
  call getmem('Tempm ','ALLOC','REAL',jtmp,nn)
  call dmxma(nbas,'C','N',inUS,Work(jK),Work(jS),Two*clight)
  call dmxma(nbas,'N','N',Work(jS),inUS,Work(jH),-Two*clight)
  call dmxma(nbas,'N','N',Work(jS),inUL,Work(jtmp),One)
  call daxpy_(nbas*nbas,One,Work(jtmp),1,Work(jH),1)
  call dmxma(nbas,'C','N',inUL,Work(jK),Work(jS),One)
  call dmxma(nbas,'N','N',Work(jS),inUS,Work(jtmp),Two*clight)
  call daxpy_(nbas*nbas,One,Work(jtmp),1,Work(jH),1)

  call dmxma(nbas,'C','N',inUL,Work(jV),Work(jS),One)
  call dmxma(nbas,'N','N',Work(jS),inUL,Work(jtmp),One)
  call daxpy_(nbas*nbas,One,Work(jtmp),1,Work(jH),1)
  call dmxma(nbas,'C','N',inUS,Work(jpVp),Work(jS),One)
  call dmxma(nbas,'N','N',Work(jS),inUS,Work(jtmp),One)
  call daxpy_(nbas*nbas,One,Work(jtmp),1,Work(jH),1)
  call getmem('Tempm ','FREE','REAL',jtmp,nn)
  ks = 0
  do iblock=1,nbl
    mbl = ibl(iblock)
    do i=1,mbl
      do j=1,mbl
        kL = Lmap(j+ks)-1+(Lmap(i+ks)-1)*nbas
        Work(jH+kL) = Work(jB+kL)
      end do
    end do
    ks = ks+mbl
  end do
end if
do i=0,nbas*nbas-1
  Work(jV+i) = Work(jH+i)
end do

! Copy relativistic one-electron Hamiltonian back to inK

k = 0
do i=1,nbas
  do j=1,i
    k = k+1
    inK(k) = Work(jV+j-1+(i-1)*nbas)
  end do
end do

! Free temp memories

call getmem('skin ','FREE','REAL',jK,nn)
call getmem('sSS  ','FREE','REAL',jS,nn)
call getmem('sV   ','FREE','REAL',jV,nn)
call getmem('spVp ','FREE','REAL',jpVp,nn)
call getmem('sHam ','FREE','REAL',jH,nn)
call getmem('sSav ','FREE','REAL',jB,nn)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(inDX)

end subroutine XDR_Local_Ham
