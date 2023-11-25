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

subroutine permvb2_cvb(v1,iperm,vb,iapr,ixapr,v2,ialg)

use casvb_global, only: nalf, nbet, nda, ndb, ndet, ndetvb, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
! V1 is dimensioned either NDET or NDETVB according to CI/VB
! V2 is dimensioned NDET/NDA or NDETVB according to CI/VB
real(kind=wp), intent(inout) :: v1(*)
integer(kind=iwp), intent(in) :: iperm(norb), iapr(ndetvb), ixapr(nda+1)
logical(kind=iwp), intent(in) :: vb
real(kind=wp), intent(_OUT_) :: v2(*)
integer(kind=iwp), intent(inout) :: ialg
integer(kind=iwp) :: i, ia, ialf, iat, iato, iatold, ib, ibet, iboff, ibt, ibto, ibtold, inboff, indx, ineg, ioffs, ioffs1, &
                     ioffs2, iorb, iprm, ixa, ixato, rc
logical(kind=iwp) :: done
integer(kind=iwp), allocatable :: inda(:), indb(:), inewocc(:), inocc2(:), locc(:), lunocc(:), maxgrph(:), mingrph(:), negs(:), &
                                  nk(:), xalf(:,:), xbet(:,:)
real(kind=wp), allocatable :: phsa(:), phsb(:)
integer(kind=iwp), external :: indget_cvb
real(kind=wp), external :: party_cvb

call mma_allocate(negs,norb,label='negs')

! Some tests of permutation
! Valid?
negs(:) = 0
do i=1,norb
  iprm = abs(iperm(i))
  if ((iprm < 1) .or. (iprm > norb)) then
    write(u6,*) ' Illegal orbital permutation!'
    call abend_cvb()
  end if
  negs(iprm) = negs(iprm)+1
end do
do iorb=1,norb
  if (negs(iorb) /= 1) then
    write(u6,*) ' Illegal orbital permutation!'
    call abend_cvb()
  end if
end do
! Return if identity
do iorb=1,norb
  if (iperm(iorb) /= iorb) exit
end do
if (iorb <= norb) then

  call mma_allocate(mingrph,[0,norb],label='mingrph')
  call mma_allocate(maxgrph,[0,norb],label='maxgrph')
  call mma_allocate(nk,[0,norb],label='nk')
  call mma_allocate(locc,norb+1,label='locc')
  call mma_allocate(lunocc,norb+1,label='lunocc')
  call mma_allocate(inewocc,norb,label='inewocc')
  call mma_allocate(inocc2,norb,label='inocc2')
  call mma_allocate(inda,nda,label='inda')
  call mma_allocate(phsa,nda,label='phsa')
  call mma_allocate(indb,ndb,label='indb')
  call mma_allocate(phsb,ndb,label='phsb')

  ! Use IALG=2 if only phase changes
  do iorb=1,norb
    if (abs(iperm(iorb)) /= iorb) exit
  end do
  if (iorb > norb) ialg = 2
  negs(:) = 0
  do i=1,norb
    if (iperm(i) < 0) negs(abs(iperm(i))) = 1
  end do
  ! Alpha loop:
  inocc2(:) = 0
  do iorb=0,norb
    mingrph(iorb) = max(iorb-norb+nalf,0)
    maxgrph(iorb) = min(iorb,nalf)
  end do
  call mma_allocate(xalf,[0,norb],[0,nalf],label='xalf')
  call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
  nk(:) = maxgrph(:)
  call occupy_cvb(nk,norb,locc,lunocc)
  indx = 1
  do
    inewocc(:) = 0
    do ialf=1,nalf
      inewocc(abs(iperm(locc(ialf)))) = ialf
    end do
    ineg = 0
    ia = 0
    do iorb=1,norb
      if (inewocc(iorb) /= 0) then
        ia = ia+1
        inocc2(ia) = inewocc(iorb)
        inewocc(iorb) = 1
        if (negs(iorb) == 1) ineg = ineg+1
      end if
    end do
    if (mod(ineg,2) == 0) then
      phsa(indx) = party_cvb(inocc2,nalf)
    else
      phsa(indx) = -party_cvb(inocc2,nalf)
    end if
    inda(indx) = indget_cvb(inewocc,nalf,norb,xalf)

    call loind_cvb(norb,nalf,nk,mingrph,maxgrph,locc,lunocc,indx,xalf,rc)
    if (rc == 0) exit
  end do
  call mma_deallocate(xalf)
  ! Beta loop:
  inocc2(:) = 0
  do iorb=0,norb
    mingrph(iorb) = max(iorb-norb+nbet,0)
    maxgrph(iorb) = min(iorb,nbet)
  end do
  call mma_allocate(xbet,[0,norb],[0,nbet],label='xbet')
  call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)
  nk(:) = maxgrph(:)
  call occupy_cvb(nk,norb,locc,lunocc)
  indx = 1
  do
    inewocc(:) = 0
    do ibet=1,nbet
      inewocc(abs(iperm(locc(ibet)))) = ibet
    end do
    ineg = 0
    ib = 0
    do iorb=1,norb
      if (inewocc(iorb) /= 0) then
        ib = ib+1
        inocc2(ib) = inewocc(iorb)
        inewocc(iorb) = 1
        if (negs(iorb) == 1) ineg = ineg+1
      end if
    end do
    if (mod(ineg,2) == 0) then
      phsb(indx) = party_cvb(inocc2,nbet)
    else
      phsb(indx) = -party_cvb(inocc2,nbet)
    end if
    indb(indx) = indget_cvb(inewocc,nbet,norb,xbet)

    call loind_cvb(norb,nbet,nk,mingrph,maxgrph,locc,lunocc,indx,xbet,rc)
    if (rc == 0) exit
  end do
  call mma_deallocate(xbet)

  if (vb) then
    v2(1:ndetvb) = Zero
    do ia=1,nda
      iato = inda(ia)
      do ixa=ixapr(ia),ixapr(ia+1)-1
        ib = iapr(ixa)
        ibto = indb(ib)
        done = .false.
        do ixato=ixapr(iato),ixapr(iato+1)-1
          if (iapr(ixato) == ibto) then
            done = .true.
            exit
          end if
        end do
        if (.not. done) then
          ! Shouldn't get here ...
          write(u6,'(a,100i3)') ' Error, VB determinants not closed under permutation :',iperm
          call abend_cvb()
        end if
        v2(ixa) = phsa(ia)*phsb(ib)*v1(ixato)
      end do
    end do
    v1(1:ndetvb) = v2(1:ndetvb)
  else if (ialg == 1) then
    ! Brute force strategy if enough memory (x1.5 faster):
    do ib=1,ndb
      iboff = (ib-1)*nda
      inboff = (indb(ib)-1)*nda
      do ia=1,nda
        v2(ia+iboff) = phsa(ia)*phsb(ib)*v1(inda(ia)+inboff)
      end do
    end do
    v1(1:ndet) = v2(1:ndet)
  else if (ialg == 2) then
    ! More-or-less in-place update of V1:
    do ia=1,nda
      if (ia == inda(ia)) then
        if (phsa(ia) == -One) then
          ioffs = ia-nda
          do ib=1,ndb
            v1(ib*nda+ioffs) = -v1(ib*nda+ioffs)
          end do
        end if
      else if (inda(ia) /= 0) then
        ! Cyclic permutation involving IA:
        ioffs = ia-nda
        do ib=1,ndb
          v2(ib) = v1(ib*nda+ioffs)
        end do
        iat = ia
        do
          if (phsa(iat) == One) then
            ioffs1 = iat-nda
            ioffs2 = inda(iat)-nda
            do ib=1,ndb
              v1(ib*nda+ioffs1) = v1(ib*nda+ioffs2)
            end do
          else
            ioffs1 = iat-nda
            ioffs2 = inda(iat)-nda
            do ib=1,ndb
              v1(ib*nda+ioffs1) = -v1(ib*nda+ioffs2)
            end do
          end if
          iatold = iat
          iat = inda(iat)
          inda(iatold) = 0
          if (inda(iat) == ia) exit
        end do
        if (phsa(iat) == One) then
          ioffs = iat-nda
          do ib=1,ndb
            v1(ib*nda+ioffs) = v2(ib)
          end do
        else
          ioffs = iat-nda
          do ib=1,ndb
            v1(ib*nda+ioffs) = -v2(ib)
          end do
        end if
        inda(iat) = 0
      end if
    end do
    do ib=1,ndb
      if (ib == indb(ib)) then
        if (phsb(ib) == -One) then
          ioffs = (ib-1)*nda
          v1(ioffs+1:ioffs+nda) = -v1(ioffs+1:ioffs+nda)
        end if
      else if (indb(ib) /= 0) then
        ! Cyclic permutation involving IB:
        ioffs = (ib-1)*nda
        v2(1:nda) = v1(ioffs+1:ioffs+nda)
        ibt = ib
        do
          if (phsb(ibt) == One) then
            ioffs1 = (ibt-1)*nda
            ioffs2 = (indb(ibt)-1)*nda
            v1(ioffs1+1:ioffs1+nda) = v1(ioffs2+1:ioffs2+nda)
          else
            ioffs1 = (ibt-1)*nda
            ioffs2 = (indb(ibt)-1)*nda
            v1(ioffs1+1:ioffs1+nda) = -v1(ioffs2+1:ioffs2+nda)
          end if
          ibtold = ibt
          ibt = indb(ibt)
          indb(ibtold) = 0
          if (indb(ibt) == ib) exit
        end do
        if (phsb(ibt) == One) then
          ioffs = (ibt-1)*nda
          v1(ioffs+1:ioffs+nda) = v2(1:nda)
        else
          ioffs = (ibt-1)*nda
          v1(ioffs+1:ioffs+nda) = -v2(1:nda)
        end if
        indb(ibt) = 0
      end if
    end do
  end if

  call mma_deallocate(mingrph)
  call mma_deallocate(maxgrph)
  call mma_deallocate(nk)
  call mma_deallocate(locc)
  call mma_deallocate(lunocc)
  call mma_deallocate(inewocc)
  call mma_deallocate(inocc2)
  call mma_deallocate(inda)
  call mma_deallocate(phsa)
  call mma_deallocate(indb)
  call mma_deallocate(phsb)

end if

call mma_deallocate(negs)

return

end subroutine permvb2_cvb
