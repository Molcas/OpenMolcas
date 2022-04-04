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

subroutine exppsb(symp,symq,symr,syms,valn,jn,kn,ln)
! this routine realizes expansion of symmetry block
! symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL), provided such integrals exist
! It finds corresponding (IJ|KL) and expands it to open
! NORB(symp) TEMP files with a structure
! indq,indr,inds,value, each TEMP for one p
!
! N.B. This process can be accelerated, if exppbs would be
! divided into exppsb1-8, each for given typ, since this
! routine is common for all types.
!
! types of (ij|kl) NI,J,K,L defined in III
!
! 1 - si=sk, si=sj, sk=sl
! 2 - si=sk, si=sj, sk>sl
! 3 - si=sk, si>sj, sk=sl
! 4 - si=sk, si>sj, sk>sl
! 5 - si>sk, si=sj, sk=sl
! 6 - si>sk, si=sj, sk>sl
! 7 - si>sk, si>sj, sk=sl
! 8 - si>sk, si>sj, sk>sl

use ccsort_global, only: fullprint, idis, LUINTM, mbas, NORB, np, nq, nr, ns, nshow, nsize, typ
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: symp, symq, symr, syms
real(kind=wp), intent(out) :: valn(nsize,mbas)
integer(kind=iwp), intent(out) :: jn(nsize,mbas), kn(nsize,mbas), ln(nsize,mbas)
#include "tratoc.fh"
integer(kind=iwp) :: i1, idis13, ilow, ind(4), indtemp, iold, iup, j1, jlow, jold, jup, k1, kold, kup, l1, lold, lup, m3, nhelp1, &
                     nhelp2, ni, nj, nk, nl, nsi, nsj, nsk, nsl, typp, yes234, yes5, yes678
real(kind=wp) :: val1
real(kind=wp), allocatable :: TWO(:)

!I get address
idis13 = idis(symp,symq,symr)

!II preparing nshow vector

nshow(1:norb(symp)) = 0

!III.1 define order of indices

ni = np(symp,symq,symr)
nj = nq(symp,symq,symr)
nk = nr(symp,symq,symr)
nl = ns(symp,symq,symr)

!III.2 def yes1-8

typp = typ(symp,symq,symr)

!:1 combination (ij|kl) -> (ij|kl)
!   used in types: 1,2,3,4,5,6,7,8 (all)
!yes1 = 1

!:2 combination (ij|kl) -> (ji|kl)
!:3 combination (ij|kl) -> (ij|lk)
!:4 combination (ij|kl) -> (ji|lk)
!   used in types: 1,5 since 2,3,6,7 never appear
if ((typp == 1) .or. (typp == 5)) then
  yes234 = 1
else
  yes234 = 0
end if

!:5 combination (ij|kl) -> (kl|ij)
!   used in types: 1,2,3,4
if ((typp >= 1) .and. (typp <= 4)) then
  yes5 = 1
else
  yes5 = 0
end if

!:6 combination (ij|kl) -> (lk|ij)
!:7 combination (ij|kl) -> (kl|ji)
!:8 combination (ij|kl) -> (lk|ji)
!   used in types: 1 (since 2,3 never appeard)
if (typp == 1) then
  yes678 = 1
else
  yes678 = 0
end if

! define NSI,NSJ,NSK,NSL
ind(ni) = symp
ind(nj) = symq
ind(nk) = symr
ind(nl) = syms
NSI = ind(1)
NSJ = ind(2)
NSK = ind(3)
NSL = ind(4)

call mma_allocate(TWO,nTraBuf,label='TWO')

indtemp = nTraBuf+1
KUP = NORB(NSK)
do KOLD=1,KUP
  if (fullprint >= 3) write(u6,*) ' * K ind ',KOLD

  LUP = NORB(NSL)
  if (NSK == NSL) LUP = KOLD
  do LOLD=1,LUP
    if (fullprint >= 3) write(u6,*) ' ** L ind ',LOLD

    ILOW = 1
    if (NSI == NSK) ILOW = KOLD
    IUP = NORB(NSI)
    do IOLD=ILOW,IUP
      if (fullprint >= 3) write(u6,*) ' *** I ind ',IOLD

      JLOW = 1
      if ((NSI == NSK) .and. (IOLD == KOLD)) JLOW = LOLD
      JUP = NORB(NSJ)
      if (NSI == NSJ) JUP = IOLD
      do JOLD=JLOW,JUP
        if (fullprint >= 3) write(u6,*) ' **** J ind ',JOLD

        ! read block of integrals if necessary

        if (indtemp == (nTraBuf+1)) then
          indtemp = 1
          ! read block
          call dDAFILE(LUINTM,2,TWO,nTraBuf,IDIS13)
        end if

        ! write integrals to appropriate positions

        val1 = TWO(indtemp)

        !:1 combination (ij|kl) -> (ij|kl)
        !   since yes1 is always 1, if structure is skipped
        ind(1) = iold
        ind(2) = jold
        ind(3) = kold
        ind(4) = lold
        j1 = ind(nj)
        l1 = ind(nl)
        if ((symq /= syms) .or. (l1 <= j1)) then
          i1 = ind(ni)
          k1 = ind(nk)

          m3 = nshow(i1)+1
          jn(m3,i1) = j1
          kn(m3,i1) = k1
          ln(m3,i1) = l1
          valn(m3,i1) = val1
          nshow(i1) = m3

          if (m3 == nsize) then
            call zasun(i1,nsize,valn,jn,kn,ln)
            nshow(i1) = 0
          end if
        end if

        if (yes234 == 1) then

          !:2 combination (ij|kl) -> (ji|kl)
          ind(1) = jold
          ind(2) = iold
          ind(3) = kold
          ind(4) = lold
          j1 = ind(nj)
          l1 = ind(nl)
          if ((symq /= syms) .or. (l1 <= j1)) then
            i1 = ind(ni)
            k1 = ind(nk)

            m3 = nshow(i1)+1
            jn(m3,i1) = j1
            kn(m3,i1) = k1
            ln(m3,i1) = l1
            valn(m3,i1) = val1
            nshow(i1) = m3

            if (m3 == nsize) then
              call zasun(i1,nsize,valn,jn,kn,ln)
              nshow(i1) = 0
            end if
          end if

          !:3 combination (ij|kl) -> (ij|lk)
          ind(1) = iold
          ind(2) = jold
          ind(3) = lold
          ind(4) = kold
          j1 = ind(nj)
          l1 = ind(nl)
          if ((symq /= syms) .or. (l1 <= j1)) then
            i1 = ind(ni)
            k1 = ind(nk)

            m3 = nshow(i1)+1
            jn(m3,i1) = j1
            kn(m3,i1) = k1
            ln(m3,i1) = l1
            valn(m3,i1) = val1
            nshow(i1) = m3

            if (m3 == nsize) then
              call zasun(i1,nsize,valn,jn,kn,ln)
              nshow(i1) = 0
            end if
          end if

          !:4 combination (ij|kl) -> (ji|lk)
          ind(1) = jold
          ind(2) = iold
          ind(3) = lold
          ind(4) = kold
          j1 = ind(nj)
          l1 = ind(nl)
          if ((symq /= syms) .or. (l1 <= j1)) then
            i1 = ind(ni)
            k1 = ind(nk)

            m3 = nshow(i1)+1
            jn(m3,i1) = j1
            kn(m3,i1) = k1
            ln(m3,i1) = l1
            valn(m3,i1) = val1
            nshow(i1) = m3

            if (m3 == nsize) then
              call zasun(i1,nsize,valn,jn,kn,ln)
              nshow(i1) = 0
            end if
          end if

        end if

        !:5 combination (ij|kl) -> (kl|ij)
        if (yes5 == 1) then
          ind(1) = kold
          ind(2) = lold
          ind(3) = iold
          ind(4) = jold
          j1 = ind(nj)
          l1 = ind(nl)
          if ((symq /= syms) .or. (l1 <= j1)) then
            i1 = ind(ni)
            k1 = ind(nk)

            m3 = nshow(i1)+1
            jn(m3,i1) = j1
            kn(m3,i1) = k1
            ln(m3,i1) = l1
            valn(m3,i1) = val1
            nshow(i1) = m3

            if (m3 == nsize) then
              call zasun(i1,nsize,valn,jn,kn,ln)
              nshow(i1) = 0
            end if
          end if
        end if

        if (yes678 == 1) then

          !:6 combination (ij|kl) -> (lk|ij)
          ind(1) = lold
          ind(2) = kold
          ind(3) = iold
          ind(4) = jold
          j1 = ind(nj)
          l1 = ind(nl)
          if ((symq /= syms) .or. (l1 <= j1)) then
            i1 = ind(ni)
            k1 = ind(nk)

            m3 = nshow(i1)+1
            jn(m3,i1) = j1
            kn(m3,i1) = k1
            ln(m3,i1) = l1
            valn(m3,i1) = val1
            nshow(i1) = m3

            if (m3 == nsize) then
              call zasun(i1,nsize,valn,jn,kn,ln)
              nshow(i1) = 0
            end if
          end if

          !:7 combination (ij|kl) -> (kl|ji)
          ind(1) = kold
          ind(2) = lold
          ind(3) = jold
          ind(4) = iold
          j1 = ind(nj)
          l1 = ind(nl)
          if ((symq /= syms) .or. (l1 <= j1)) then
            i1 = ind(ni)
            k1 = ind(nk)

            m3 = nshow(i1)+1
            jn(m3,i1) = j1
            kn(m3,i1) = k1
            ln(m3,i1) = l1
            valn(m3,i1) = val1
            nshow(i1) = m3

            if (m3 == nsize) then
              call zasun(i1,nsize,valn,jn,kn,ln)
              nshow(i1) = 0
            end if
          end if

          !:8 combination (ij|kl) -> (lk|ji)
          ind(1) = lold
          ind(2) = kold
          ind(3) = jold
          ind(4) = iold
          j1 = ind(nj)
          l1 = ind(nl)
          if ((symq /= syms) .or. (l1 <= j1)) then
            i1 = ind(ni)
            k1 = ind(nk)

            m3 = nshow(i1)+1
            jn(m3,i1) = j1
            kn(m3,i1) = k1
            ln(m3,i1) = l1
            valn(m3,i1) = val1
            nshow(i1) = m3

            if (m3 == nsize) then
              call zasun(i1,nsize,valn,jn,kn,ln)
              nshow(i1) = 0
            end if
          end if

        end if

        indtemp = indtemp+1

      end do
    end do
  end do
end do

call mma_deallocate(TWO)

!IV write the rest integrals if needed

do nhelp1=1,norb(symp)
  nhelp2 = nshow(nhelp1)
  if (nhelp2 > 0) call zasun(nhelp1,nhelp2,valn,jn,kn,ln)
end do

return

end subroutine exppsb
