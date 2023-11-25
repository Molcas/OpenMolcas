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
! Copyright (C) 2019, Thomas J. Duignan                                *
!               2021, Rulin Feng                                       *
!***********************************************************************

subroutine repmat(idbg,bInt,sInt,donorm)

use Basis_Info, only: dbsc, nBas, nCnttp
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: idbg
real(kind=wp), intent(in) :: bInt(*)
real(kind=wp), intent(_OUT_) :: sInt(*)
logical(kind=iwp), intent(in) :: donorm
#include "Molcas.fh"
#include "rinfo.fh"
integer(kind=iwp) :: ip, istart, jp, kp, nAngrMax, nAtomTot, nc, np, nrBasTot
real(kind=wp) :: finish, kpp
logical(kind=iwp) :: New_Center, New_l, New_m, Old_Center, Old_l
integer(kind=iwp), allocatable :: icaddr(:), ihelp(:,:), mcaddr(:), numb(:), numc(:)
real(kind=wp), allocatable :: fin(:,:), mag(:,:), pa(:), scr(:,:), u2c(:,:), u2ct(:,:)
integer(kind=iwp) :: i, ia, iBas, ibasL, ic, icnt, iCnttp, icon, iCont, ilarge, indb, indbL, ipbasL, iPrim, iPrint, ismall, iSym, &
                     j, jbas, jbasL, jcon, jpbasL, jPrim, k, ka, kbias, kc, kcL, la, mp, nSym, numck, numcl
real(kind=wp) :: sum_

! contracted basis, atomic basis functions
!
! symmetry info
!
! lant(i): number of atoms in i:th symmetry bf
! expand the coefficient matrix into symmetry basis set
! auxiliary
! icaddr(i): adresses in coeff for a symmetry adapted function

!do i=1,12640
!  write(67,'(es25.14)') bint(i)
!end do
nSym = nIrrep
iPrint = 0
if ((iprint >= 10) .or. (idbg > 0)) then
  write(idbg,*) ' in repmat',nsym
  write(idbg,*) nSym,(nBas(i),i=0,nsym-1)
  write(idbg,*) nSym,(nrBas(i),i=1,nsym)
end if

nAtomTot = 0
do iCnttp=1,nCnttp
  nAtomTot = nAtomTot+dbsc(iCnttp)%nCntr
end do
nAngrMax = 0
do ia=1,nAtomTot
  nAngrMax = max(nAngrMax,nAngr(ia)+1)
end do
call mma_allocate(ihelp,nAtomTot,nAngrMax,label='ihelp')

! set up pointer

k = 0
ia = 0  ! center index
ka = 0  ! shell index
do iCnttp=1,nCnttp
  do icnt=1,dbsc(iCnttp)%nCntr
    ia = ia+1
    do la=1,nAngr(ia)+1
      ka = ka+1
      ihelp(ia,la) = k
      k = k+nPrimr(ka)*nBasisr(ka)
    end do
  end do
end do
if ((iPrint >= 10) .or. (idbg > 0)) then
  write(idbg,*) ' Help vector'
  ia = 0
  do iCnttp=1,nCnttp
    do icnt=1,dbsc(iCnttp)%nCntr
      ia = ia+1
      write(idbg,'(10i5)') (ihelp(ia,j),j=1,nAngr(ia)+1)
    end do
  end do
end if

nrBasTot = 0
do iSym=1,nSym
  nrBasTot = nrBasTot+nrBas(iSym)
end do
call mma_allocate(icaddr,nrBasTot,label='icaddr')
call mma_allocate(mcaddr,nrBasTot,label='mcaddr')
call mma_allocate(numb,nrBasTot,label='numb')
call mma_allocate(numc,nrBasTot,label='numc')
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over irreps

k = 0
do iSym=1,nSym
  numck = 1
  numcl = 0
  kbias = 0

  ! Loop over basis functions in irrep

  do iCont=1,nrBas(iSym)
    k = k+1
    if (iCont > 1) then
      New_Center = icent(k) /= icent(k-1)
      New_l = lnang(k) /= lnang(k-1)
      New_m = lmag(k) /= lmag(k-1)
      if (New_m) kbias = kbias-numc(k-1)
      if ((New_Center) .or. New_l) kbias = 0
    end if
    kbias = kbias+1
    ka = 0
    ia = 0
    do iCnttp=1,nCnttp
      do iCnt=1,dbsc(iCnttp)%nCntr
        ia = ia+1
        do la=1,nAngr(ia)+1
          ka = ka+1
          Old_Center = icent(k) == ia
          Old_l = lnang(k) == (la-1)
          if (idbg > 0) write(idbg,*) ' at numck',k,ia,icent(k),la-1,lnang(k),ia,New_Center,New_l
          if (Old_Center .and. Old_l) then
            numc(k) = nBasisr(ka)
            numb(k) = nPrimr(ka)
            icaddr(k) = ihelp(ia,la)+(kbias-1)*nPrimr(ka)
            if ((k > 1) .and. (kbias == 1)) numck = numcl+numck
            mcaddr(k) = numck
            numcl = nPrimr(ka)
          end if
        end do  ! la
      end do    ! iCnt
    end do      ! iCnttp
  end do        ! iCont
end do          ! iSym
call mma_deallocate(ihelp)
k = 0
if ((iPrint >= 10) .or. (idbg > 0)) then
  ic = 0
  ip = 0
  do iSym=1,nSym
    write(idbg,*) ' symmetry',iSym
    write(idbg,*) ' numb'
    write(idbg,'(20i4)') (numb(i+ic),i=1,nrBas(iSym))
    write(idbg,*) ' numc'
    write(idbg,'(20i4)') (numc(i+ic),i=1,nrBas(iSym))
    write(idbg,*) ' Pointer to contraction vector'
    write(idbg,'(20i4)') (icaddr(i+ic),i=1,nrBas(iSym))
    write(idbg,*) ' mcaddr'
    write(idbg,'(20i4)') (mcaddr(i+ic),i=1,nrBas(iSym))
    ic = ic+nrBas(iSym)
    ip = ip+nBas(iSym-1)
  end do
end if
! debugdebug
!write(u6,*) (rCof(i),i=1,4)
! debugdebug

! transform

if (donorm) then
  kc = 0
  kcL = 0
  ibasL = 0
  indbL = 0
  kp = 0
  do iSym=1,nSym
    ! loop over contracted
    do iBas=1,nrBas(iSym)
      ibasL = ibasL+1
      jbasL = kcL
      do jbas=1,ibas
        jbasL = jbasL+1
        kc = kc+1
        sum_ = Zero
        ipbasL = mcaddr(ibasL)-1
        ! loop over primitives in this contracted function
        do iPrim=1,numb(ibasL)
          ipbasL = ipbasL+1
          jpbasL = mcaddr(jbasL)-1
          do jPrim=1,numb(jbasL)
            jpbasL = jpbasL+1
            ilarge = max(ipbasL,jpbasL)
            ismall = min(ipbasL,jpbasL)
            kp = kp+1
            indb = indbL+(ilarge*(ilarge-1))/2+ismall
            sum_ = sum_+bint(indb)*rCof(icaddr(ibasL)+iPrim)*rCof(icaddr(jbasL)+jPrim)
            if (idbg > 0) then
              write(idbg,*) indb,icaddr(ibasL)+iPrim,icaddr(jbasL)+jPrim
              write(idbg,*) bint(indb),rCof(icaddr(ibasL)+iPrim),rCof(icaddr(jbasL)+jPrim)
            end if
          end do  ! jprim
        end do    ! iprim
        sint(kc) = sum_
        !write(66,'(es25.14)') sum_
      end do
    end do
    kcL = kcL+nrBas(iSym)
    indbL = indbL+(nBas(iSym-1)*(nBas(iSym-1)+1))/2
    !if (idbg > 0) write(idbg,*) ipbasL,jpbasL,kp
  end do

else

  np = nBas(0)
  nc = nrBas(1)
  ! Scratch to square the mag ints and contract
  call mma_allocate(mag,np,np)
  call mma_allocate(u2c,np,nc)
  call mma_allocate(u2ct,nc,np)
  call mma_allocate(scr,np,nc)
  call mma_allocate(fin,nc,nc)
  call mma_allocate(pa,np)

  mag(:,:) = Zero
  u2c(:,:) = Zero
  u2ct(:,:) = Zero
  scr(:,:) = Zero
  fin(:,:) = Zero
  pa(:) = Zero

  ! Square the uncontracted ints
  mp = 0
  do ip=1,np
    do jp=1,ip
      mp = mp+1
      mag(jp,ip) = bint(mp)
      ! Suboptimal
      mag(ip,jp) = bint(mp)
    end do
  end do

  ! Construct contraction matrix
  do icon=1,nc
    istart = mcaddr(icon)-1
    do iprim=1,numb(icon)
      !write(u6,'(A,4I5,F14.10)') 'u2c loop',icon,iprim,istart,icaddr(icon),rCof(icaddr(icon)+iprim)
      u2c(istart+iprim,icon) = rCof(icaddr(icon)+iprim)
    end do
  end do

  ! U2C.T * UNCON * U2C => CON
  u2ct(:,:) = transpose(u2c)

  ! Since dgemm_ is causing trouble and this is trivial
  do iprim=1,np
    do icon=1,nc
      kpp = Zero
      pa(:) = u2ct(icon,:)*mag(:,iprim)
      do ip=1,np
        kpp = kpp+pa(ip)
      end do
      scr(iprim,icon) = kpp
    end do
  end do

  kp = 0
  do icon=1,nc
    do jcon=1,nc
      pa(:) = u2c(:,icon)*scr(:,jcon)
      kpp = Zero
      do ip=1,np
        kpp = kpp+pa(ip)
      end do
      fin(icon,jcon) = kpp
    end do
  end do

  kp = 0
  do icon=1,nc
    do jcon=1,icon
      kp = kp+1
      sint(kp) = fin(icon,jcon)
    end do
  end do

  call mma_deallocate(mag)
  call mma_deallocate(u2c)
  call mma_deallocate(u2ct)
  call mma_deallocate(scr)
  call mma_deallocate(fin)
  call mma_deallocate(pa)

  call cpu_time(finish)

end if !donorm

call mma_deallocate(icaddr)
call mma_deallocate(mcaddr)
call mma_deallocate(numb)
call mma_deallocate(numc)

return

end subroutine repmat
