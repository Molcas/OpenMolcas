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

!***********************************************************************
!*                                                                     *
!*  DEV2B  := calculate two-electron Hessian                           *
!*                                                                     *
!***********************************************************************
subroutine dev2b_2_cvb(v1,v2,cfrom,hessorb,hesst,oaa2,aa1,gx,grad2)
! Calculates V1 EijEkl CFROM and V2 EijEkl CFROM

use casvb_global, only: absym, i1alf, i1bet, iafrm, iapr, iato, ibfrm, ibpr, ibto, ixapr, ixbpr, n1a, n1b, nda, ndb, norb, nprorb, &
                        phato, phbto, projcas, sc
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: v1(nda,ndb), v2(nda,ndb), cfrom(nda,ndb), oaa2, aa1, gx(norb,norb), grad2(nprorb)
real(kind=wp), intent(inout) :: hessorb(nprorb,nprorb)
real(kind=wp), intent(out) :: hesst(norb*norb,norb*norb)
integer(kind=iwp) :: ia, iax, iaxtmp, ib, ibx, ibxtmp, iorb, ip1, ip2, itmp, ixa, ixb, jax, jbx, ji, jorb, kax, kbx, korb, lax, &
                     lbx, li, lk, lorb
real(kind=wp) :: phase, res1, res2, t1, t2, tcof, term

do ip1=1,norb*norb
  iorb = (ip1-1)/norb+1
  jorb = ip1-(iorb-1)*norb
  do ip2=1,ip1
    korb = (ip2-1)/norb+1
    lorb = ip2-(korb-1)*norb
    res1 = Zero
    res2 = Zero
    if (projcas .and. (.not. sc)) then
      ! 1) Alpha excitation
      do ia=1,n1a
        iaxtmp = i1alf(ia,iorb)
        jax = iato(jorb,iaxtmp)
        if (jax /= 0) then
          iax = iato(iorb,iaxtmp)
          tcof = phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
          ! I -> J
          do ixa=ixapr(iax),ixapr(iax+1)-1
            ibx = iapr(ixa)
            ! 2. alpha k -> l
            itmp = iafrm(korb,jax)
            kax = iato(lorb,itmp)
            if (kax /= 0) then
              phase = phato(korb,itmp)*phato(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1+v1(kax,ibx)*term
              res2 = res2+v2(kax,ibx)*term
            end if
            ! 2. alpha l -> k
            itmp = iafrm(lorb,iax)
            lax = iato(korb,itmp)
            if (lax /= 0) then
              phase = phato(lorb,itmp)*phato(korb,itmp)
              term = tcof*phase*cfrom(lax,ibx)
              res1 = res1-v1(jax,ibx)*term
              res2 = res2-v2(jax,ibx)*term
            end if
            ! 2. beta  k -> l
            itmp = ibfrm(korb,ibx)
            kbx = ibto(lorb,itmp)
            if (kbx /= 0) then
              phase = phbto(korb,itmp)*phbto(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1+v1(jax,kbx)*term
              res2 = res2+v2(jax,kbx)*term
            end if
            ! 2. beta  l -> k
            itmp = ibfrm(lorb,ibx)
            lbx = ibto(korb,itmp)
            if (lbx /= 0) then
              phase = phbto(lorb,itmp)*phbto(korb,itmp)
              term = tcof*phase*cfrom(iax,lbx)
              res1 = res1-v1(jax,ibx)*term
              res2 = res2-v2(jax,ibx)*term
            end if
          end do
          ! J -> I
          do ixa=ixapr(jax),ixapr(jax+1)-1
            ibx = iapr(ixa)
            ! 2. alpha k -> l
            itmp = iafrm(korb,jax)
            kax = iato(lorb,itmp)
            if (kax /= 0) then
              phase = phato(korb,itmp)*phato(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1-v1(kax,ibx)*term
              res2 = res2-v2(kax,ibx)*term
            end if
            ! 2. alpha l -> k
            itmp = iafrm(lorb,iax)
            lax = iato(korb,itmp)
            if (lax /= 0) then
              phase = phato(lorb,itmp)*phato(korb,itmp)
              term = tcof*phase*cfrom(lax,ibx)
              res1 = res1+v1(jax,ibx)*term
              res2 = res2+v2(jax,ibx)*term
            end if
            ! 2. beta  k -> l
            itmp = ibfrm(korb,ibx)
            kbx = ibto(lorb,itmp)
            if (kbx /= 0) then
              phase = phbto(korb,itmp)*phbto(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1-v1(jax,kbx)*term
              res2 = res2-v2(jax,kbx)*term
            end if
            ! 2. beta  l -> k
            itmp = ibfrm(lorb,ibx)
            lbx = ibto(korb,itmp)
            if (lbx /= 0) then
              phase = phbto(lorb,itmp)*phbto(korb,itmp)
              term = tcof*phase*cfrom(iax,lbx)
              res1 = res1+v1(jax,ibx)*term
              res2 = res2+v2(jax,ibx)*term
            end if
          end do
        end if
      end do

      if (absym(3)) then
        res1 = Two*res1
        res2 = Two*res2
      else
        ! 2) Beta excitation
        do ib=1,n1b
          ibxtmp = i1bet(ib,iorb)
          jbx = ibto(jorb,ibxtmp)
          if (jbx /= 0) then
            ibx = ibto(iorb,ibxtmp)
            tcof = phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
            ! I -> J
            do ixb=ixbpr(ibx),ixbpr(ibx+1)-1
              iax = ibpr(ixb)
              ! 2. beta  k -> l
              itmp = ibfrm(korb,jbx)
              kbx = ibto(lorb,itmp)
              if (kbx /= 0) then
                phase = phbto(korb,itmp)*phbto(lorb,itmp)
                term = tcof*phase*cfrom(iax,ibx)
                res1 = res1+v1(iax,kbx)*term
                res2 = res2+v2(iax,kbx)*term
              end if
              ! 2. beta  l -> k
              itmp = ibfrm(lorb,ibx)
              lbx = ibto(korb,itmp)
              if (lbx /= 0) then
                phase = phbto(lorb,itmp)*phbto(korb,itmp)
                term = tcof*phase*cfrom(iax,lbx)
                res1 = res1-v1(iax,jbx)*term
                res2 = res2-v2(iax,jbx)*term
              end if
              ! 2. alpha k -> l
              itmp = iafrm(korb,iax)
              kax = iato(lorb,itmp)
              if (kax /= 0) then
                phase = phato(korb,itmp)*phato(lorb,itmp)
                term = tcof*phase*cfrom(iax,ibx)
                res1 = res1+v1(kax,jbx)*term
                res2 = res2+v2(kax,jbx)*term
              end if
              ! 2. alpha l -> k
              itmp = iafrm(lorb,iax)
              lax = iato(korb,itmp)
              if (lax /= 0) then
                phase = phato(lorb,itmp)*phato(korb,itmp)
                term = tcof*phase*cfrom(lax,ibx)
                res1 = res1-v1(iax,jbx)*term
                res2 = res2-v2(iax,jbx)*term
              end if
            end do
            ! J -> I
            do ixb=ixbpr(jbx),ixbpr(jbx+1)-1
              iax = ibpr(ixb)
              ! 2. beta  k -> l
              itmp = ibfrm(korb,jbx)
              kbx = ibto(lorb,itmp)
              if (kbx /= 0) then
                phase = phbto(korb,itmp)*phbto(lorb,itmp)
                term = tcof*phase*cfrom(iax,ibx)
                res1 = res1-v1(iax,kbx)*term
                res2 = res2-v2(iax,kbx)*term
              end if
              ! 2. beta  l -> k
              itmp = ibfrm(lorb,ibx)
              lbx = ibto(korb,itmp)
              if (lbx /= 0) then
                phase = phbto(lorb,itmp)*phbto(korb,itmp)
                term = tcof*phase*cfrom(iax,lbx)
                res1 = res1+v1(iax,jbx)*term
                res2 = res2+v2(iax,jbx)*term
              end if
              ! 2. alpha k -> l
              itmp = iafrm(korb,iax)
              kax = iato(lorb,itmp)
              if (kax /= 0) then
                phase = phato(korb,itmp)*phato(lorb,itmp)
                term = tcof*phase*cfrom(iax,ibx)
                res1 = res1-v1(kax,jbx)*term
                res2 = res2-v2(kax,jbx)*term
              end if
              ! 2. alpha l -> k
              itmp = iafrm(lorb,iax)
              lax = iato(korb,itmp)
              if (lax /= 0) then
                phase = phato(lorb,itmp)*phato(korb,itmp)
                term = tcof*phase*cfrom(lax,ibx)
                res1 = res1+v1(iax,jbx)*term
                res2 = res2+v2(iax,jbx)*term
              end if
            end do
          end if
        end do
      end if
    else if ((.not. projcas) .and. (.not. sc)) then
      ! 1) Alpha excitation
      do ia=1,n1a
        iaxtmp = i1alf(ia,iorb)
        jax = iato(jorb,iaxtmp)
        if (jax /= 0) then
          iax = iato(iorb,iaxtmp)
          tcof = phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
          ! I -> J
          do ixa=ixapr(iax),ixapr(iax+1)-1
            ibx = iapr(ixa)
            ! 2. alpha k -> l
            itmp = iafrm(korb,jax)
            kax = iato(lorb,itmp)
            if (kax /= 0) then
              phase = phato(korb,itmp)*phato(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1+v1(kax,ibx)*term
              res2 = res2+v2(kax,ibx)*term
            end if
            ! 2. beta  k -> l
            itmp = ibfrm(korb,ibx)
            kbx = ibto(lorb,itmp)
            if (kbx /= 0) then
              phase = phbto(korb,itmp)*phbto(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1+v1(jax,kbx)*term
              res2 = res2+v2(jax,kbx)*term
            end if
          end do
        end if
      end do

      if (absym(3)) then
        res1 = Two*res1
        res2 = Two*res2
      else
        ! 2) Beta excitation
        do ib=1,n1b
          ibxtmp = i1bet(ib,iorb)
          jbx = ibto(jorb,ibxtmp)
          if (jbx /= 0) then
            ibx = ibto(iorb,ibxtmp)
            tcof = phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
            ! I -> J
            do ixb=ixbpr(ibx),ixbpr(ibx+1)-1
              iax = ibpr(ixb)
              ! 2. beta  k -> l
              itmp = ibfrm(korb,jbx)
              kbx = ibto(lorb,itmp)
              if (kbx /= 0) then
                phase = phbto(korb,itmp)*phbto(lorb,itmp)
                term = tcof*phase*cfrom(iax,ibx)
                res1 = res1+v1(iax,kbx)*term
                res2 = res2+v2(iax,kbx)*term
              end if
              ! 2. alpha k -> l
              itmp = iafrm(korb,iax)
              kax = iato(lorb,itmp)
              if (kax /= 0) then
                phase = phato(korb,itmp)*phato(lorb,itmp)
                term = tcof*phase*cfrom(iax,ibx)
                res1 = res1+v1(kax,jbx)*term
                res2 = res2+v2(kax,jbx)*term
              end if
            end do
          end if
        end do
      end if
    else if (projcas .and. sc) then
      ! 1) Alpha excitation
      do ia=1,n1a
        iaxtmp = i1alf(ia,iorb)
        jax = iato(jorb,iaxtmp)
        if (jax /= 0) then
          iax = iato(iorb,iaxtmp)
          tcof = phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
          ! I -> J
          ibx = ndb-iax+1
          ! 2. alpha k -> l
          itmp = iafrm(korb,jax)
          kax = iato(lorb,itmp)
          if (kax /= 0) then
            phase = phato(korb,itmp)*phato(lorb,itmp)
            term = tcof*phase*cfrom(iax,ibx)
            res1 = res1+v1(kax,ibx)*term
            res2 = res2+v2(kax,ibx)*term
          end if
          ! 2. alpha l -> k
          itmp = iafrm(lorb,iax)
          lax = iato(korb,itmp)
          if (lax /= 0) then
            phase = phato(lorb,itmp)*phato(korb,itmp)
            term = tcof*phase*cfrom(lax,ibx)
            res1 = res1-v1(jax,ibx)*term
            res2 = res2-v2(jax,ibx)*term
          end if
          ! 2. beta  k -> l
          itmp = ibfrm(korb,ibx)
          kbx = ibto(lorb,itmp)
          if (kbx /= 0) then
            phase = phbto(korb,itmp)*phbto(lorb,itmp)
            term = tcof*phase*cfrom(iax,ibx)
            res1 = res1+v1(jax,kbx)*term
            res2 = res2+v2(jax,kbx)*term
          end if
          ! 2. beta  l -> k
          itmp = ibfrm(lorb,ibx)
          lbx = ibto(korb,itmp)
          if (lbx /= 0) then
            phase = phbto(lorb,itmp)*phbto(korb,itmp)
            term = tcof*phase*cfrom(iax,lbx)
            res1 = res1-v1(jax,ibx)*term
            res2 = res2-v2(jax,ibx)*term
          end if
          ! J -> I
          ibx = ndb-jax+1
          ! 2. alpha k -> l
          itmp = iafrm(korb,jax)
          kax = iato(lorb,itmp)
          if (kax /= 0) then
            phase = phato(korb,itmp)*phato(lorb,itmp)
            term = tcof*phase*cfrom(iax,ibx)
            res1 = res1-v1(kax,ibx)*term
            res2 = res2-v2(kax,ibx)*term
          end if
          ! 2. alpha l -> k
          itmp = iafrm(lorb,iax)
          lax = iato(korb,itmp)
          if (lax /= 0) then
            phase = phato(lorb,itmp)*phato(korb,itmp)
            term = tcof*phase*cfrom(lax,ibx)
            res1 = res1+v1(jax,ibx)*term
            res2 = res2+v2(jax,ibx)*term
          end if
          ! 2. beta  k -> l
          itmp = ibfrm(korb,ibx)
          kbx = ibto(lorb,itmp)
          if (kbx /= 0) then
            phase = phbto(korb,itmp)*phbto(lorb,itmp)
            term = tcof*phase*cfrom(iax,ibx)
            res1 = res1-v1(jax,kbx)*term
            res2 = res2-v2(jax,kbx)*term
          end if
          ! 2. beta  l -> k
          itmp = ibfrm(lorb,ibx)
          lbx = ibto(korb,itmp)
          if (lbx /= 0) then
            phase = phbto(lorb,itmp)*phbto(korb,itmp)
            term = tcof*phase*cfrom(iax,lbx)
            res1 = res1+v1(jax,ibx)*term
            res2 = res2+v2(jax,ibx)*term
          end if
        end if
      end do

      if (absym(3)) then
        res1 = Two*res1
        res2 = Two*res2
      else
        ! 2) Beta excitation
        do ib=1,n1b
          ibxtmp = i1bet(ib,iorb)
          jbx = ibto(jorb,ibxtmp)
          if (jbx /= 0) then
            ibx = ibto(iorb,ibxtmp)
            tcof = phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
            ! I -> J
            iax = nda-ibx+1
            ! 2. beta  k -> l
            itmp = ibfrm(korb,jbx)
            kbx = ibto(lorb,itmp)
            if (kbx /= 0) then
              phase = phbto(korb,itmp)*phbto(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1+v1(iax,kbx)*term
              res2 = res2+v2(iax,kbx)*term
            end if
            ! 2. beta  l -> k
            itmp = ibfrm(lorb,ibx)
            lbx = ibto(korb,itmp)
            if (lbx /= 0) then
              phase = phbto(lorb,itmp)*phbto(korb,itmp)
              term = tcof*phase*cfrom(iax,lbx)
              res1 = res1-v1(iax,jbx)*term
              res2 = res2-v2(iax,jbx)*term
            end if
            ! 2. alpha k -> l
            itmp = iafrm(korb,iax)
            kax = iato(lorb,itmp)
            if (kax /= 0) then
              phase = phato(korb,itmp)*phato(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1+v1(kax,jbx)*term
              res2 = res2+v2(kax,jbx)*term
            end if
            ! 2. alpha l -> k
            itmp = iafrm(lorb,iax)
            lax = iato(korb,itmp)
            if (lax /= 0) then
              phase = phato(lorb,itmp)*phato(korb,itmp)
              term = tcof*phase*cfrom(lax,ibx)
              res1 = res1-v1(iax,jbx)*term
              res2 = res2-v2(iax,jbx)*term
            end if
            ! J -> I
            iax = nda-jbx+1
            ! 2. beta  k -> l
            itmp = ibfrm(korb,jbx)
            kbx = ibto(lorb,itmp)
            if (kbx /= 0) then
              phase = phbto(korb,itmp)*phbto(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1-v1(iax,kbx)*term
              res2 = res2-v2(iax,kbx)*term
            end if
            ! 2. beta  l -> k
            itmp = ibfrm(lorb,ibx)
            lbx = ibto(korb,itmp)
            if (lbx /= 0) then
              phase = phbto(lorb,itmp)*phbto(korb,itmp)
              term = tcof*phase*cfrom(iax,lbx)
              res1 = res1+v1(iax,jbx)*term
              res2 = res2+v2(iax,jbx)*term
            end if
            ! 2. alpha k -> l
            itmp = iafrm(korb,iax)
            kax = iato(lorb,itmp)
            if (kax /= 0) then
              phase = phato(korb,itmp)*phato(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1-v1(kax,jbx)*term
              res2 = res2-v2(kax,jbx)*term
            end if
            ! 2. alpha l -> k
            itmp = iafrm(lorb,iax)
            lax = iato(korb,itmp)
            if (lax /= 0) then
              phase = phato(lorb,itmp)*phato(korb,itmp)
              term = tcof*phase*cfrom(lax,ibx)
              res1 = res1+v1(iax,jbx)*term
              res2 = res2+v2(iax,jbx)*term
            end if
          end if
        end do
      end if
    else if ((.not. projcas) .and. sc) then
      ! 1) Alpha excitation
      do ia=1,n1a
        iaxtmp = i1alf(ia,iorb)
        jax = iato(jorb,iaxtmp)
        if (jax /= 0) then
          iax = iato(iorb,iaxtmp)
          tcof = phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
          ! I -> J
          ibx = ndb-iax+1
          ! 2. alpha k -> l
          itmp = iafrm(korb,jax)
          kax = iato(lorb,itmp)
          if (kax /= 0) then
            phase = phato(korb,itmp)*phato(lorb,itmp)
            term = tcof*phase*cfrom(iax,ibx)
            res1 = res1+v1(kax,ibx)*term
            res2 = res2+v2(kax,ibx)*term
          end if
          ! 2. beta  k -> l
          itmp = ibfrm(korb,ibx)
          kbx = ibto(lorb,itmp)
          if (kbx /= 0) then
            phase = phbto(korb,itmp)*phbto(lorb,itmp)
            term = tcof*phase*cfrom(iax,ibx)
            res1 = res1+v1(jax,kbx)*term
            res2 = res2+v2(jax,kbx)*term
          end if
        end if
      end do

      if (absym(3)) then
        res1 = Two*res1
        res2 = Two*res2
      else
        ! 2) Beta excitation
        do ib=1,n1b
          ibxtmp = i1bet(ib,iorb)
          jbx = ibto(jorb,ibxtmp)
          if (jbx /= 0) then
            ibx = ibto(iorb,ibxtmp)
            tcof = phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
            ! I -> J
            iax = nda-ibx+1
            ! 2. beta  k -> l
            itmp = ibfrm(korb,jbx)
            kbx = ibto(lorb,itmp)
            if (kbx /= 0) then
              phase = phbto(korb,itmp)*phbto(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1+v1(iax,kbx)*term
              res2 = res2+v2(iax,kbx)*term
            end if
            ! 2. alpha k -> l
            itmp = iafrm(korb,iax)
            kax = iato(lorb,itmp)
            if (kax /= 0) then
              phase = phato(korb,itmp)*phato(lorb,itmp)
              term = tcof*phase*cfrom(iax,ibx)
              res1 = res1+v1(kax,jbx)*term
              res2 = res2+v2(kax,jbx)*term
            end if
          end if
        end do
      end if
    end if
    hesst(ip2,ip1) = res1
    if (ip1 /= ip2) then
      ! E_ji E_lk = E_lk E_ji - \delta_kj E_li + \delta_il E_jk
      hesst(ip1,ip2) = res1
      if (korb == jorb) hesst(ip1,ip2) = hesst(ip1,ip2)-gx(lorb,iorb)
      if (iorb == lorb) hesst(ip1,ip2) = hesst(ip1,ip2)+gx(jorb,korb)
    end if
    if ((iorb /= jorb) .and. (korb /= lorb)) then
      t1 = res1
      t2 = res2
      if ((jorb == korb) .and. (iorb /= lorb)) then
        t1 = t1-gx(lorb,iorb)
        li = lorb+(iorb-1)*(norb-1)
        if (lorb > iorb) li = li-1
        t2 = t2-grad2(li)
      end if
      ji = jorb+(iorb-1)*(norb-1)
      if (jorb > iorb) ji = ji-1
      lk = lorb+(korb-1)*(norb-1)
      if (lorb > korb) lk = lk-1
      hessorb(ji,lk) = hessorb(ji,lk)+oaa2*t1+aa1*t2
      hessorb(lk,ji) = hessorb(ji,lk)
    end if
  end do
end do

return

end subroutine dev2b_2_cvb
