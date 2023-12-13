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

subroutine dev2a_2_cvb(v1,v2,cfrom,hessorb,oaa2,aa1)
! Calculates V1 EijEkl CFROM and V2 EijEkl CFROM

use casvb_global, only: absym, i1alf, i1bet, iafrm, iapr, iato, ibfrm, ibpr, ibto, ixapr, ixbpr, n1a, n1b, nda, ndb, norb, nprorb, &
                        phato, phbto, projcas, sc
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: v1(nda,ndb), v2(nda,ndb), cfrom(nda,ndb), oaa2, aa1
real(kind=wp), intent(inout) :: hessorb(nprorb,nprorb)
integer(kind=iwp) :: i1, i2, i3, i4, ia, iax, iaxtmp, ib, ibx, ibxtmp, iorb, iprm1, iprm2, itmp, ixa, ixb, jax, jbx, jorb, kax, &
                     kbx, korb, lax, lbx, lorb
real(kind=wp) :: phase, res1, res2, tcof, term

do iprm1=1,nprorb
  i1 = (iprm1-1)/(norb-1)+1
  i2 = iprm1-(i1-1)*(norb-1)
  if (i2 >= i1) i2 = i2+1
  do iprm2=1,iprm1
    i3 = (iprm2-1)/(norb-1)+1
    i4 = iprm2-(i3-1)*(norb-1)
    if (i4 >= i3) i4 = i4+1
    if ((i1 /= i4) .or. (i2 == i3)) then
      iorb = i3
      jorb = i4
      korb = i1
      lorb = i2
    else
      iorb = i1
      jorb = i2
      korb = i3
      lorb = i4
    end if
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
    hessorb(iprm1,iprm2) = hessorb(iprm1,iprm2)+oaa2*res1+aa1*res2
    hessorb(iprm2,iprm1) = hessorb(iprm1,iprm2)
  end do
end do

return

end subroutine dev2a_2_cvb
