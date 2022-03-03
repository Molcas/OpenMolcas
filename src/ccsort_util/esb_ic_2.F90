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

subroutine esb_ic_2(symp,symq,Vic,dimp,dimq,pqind)
! this routine realizes expansion of symmetry block
! symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL),
! (for case symp=symr,symq=syms)
! provided such integrals exist
! It finds corresponding (IJ|KL) and expands it to
! matrix vic (pr,qs)

use ccsort_global, only: idis, LUINTM, mbas, NORB, np, nq, nr, ns, typ
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: symp, symq, dimp, dimq
real(kind=wp), intent(_OUT_) :: Vic(dimp*(dimp+1)/2,dimq*(dimq+1)/2)
integer(kind=iwp), intent(out) :: pqind(mbas,mbas)
#include "tratoc.fh"
integer(kind=iwp) :: i, i1, idis13, ilow, ind(4), indtemp, iold, iup, j, j1, jlow, jold, jup, k1, kold, kup, l1, lold, lup, maxx, &
                     ni, nj, nk, nl, nsi, nsj, nsk, nsl, typp, yes234, yes5, yes678
real(kind=wp) :: val1
real(kind=wp), allocatable :: TWO(:)

!I calc pqind

if (dimp >= dimq) then
  maxx = dimp
else
  maxx = dimq
end if

do i=1,maxx
  do j=1,maxx
    if (i >= j) then
      pqind(i,j) = i*(i-1)/2+j
    else
      pqind(i,j) = j*(j-1)/2+i
    end if
  end do
end do

!II get address
idis13 = idis(symp,symq,symp)

!III.1 define order of indices

ni = np(symp,symq,symp)
nj = nq(symp,symq,symp)
nk = nr(symp,symq,symp)
nl = ns(symp,symq,symp)

!III.2 def yes1-8

typp = typ(symp,symq,symp)

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
ind(nk) = symp
ind(nl) = symq
NSI = ind(1)
NSJ = ind(2)
NSK = ind(3)
NSL = ind(4)

call mma_allocate(TWO,nTraBuf,label='TWO')

indtemp = nTraBuf+1
KUP = NORB(NSK)
do KOLD=1,KUP

  LUP = NORB(NSL)
  if (NSK == NSL) LUP = KOLD
  do LOLD=1,LUP

    ILOW = 1
    if (NSI == NSK) ILOW = KOLD
    IUP = NORB(NSI)
    do IOLD=ILOW,IUP

      JLOW = 1
      if ((NSI == NSK) .and. (IOLD == KOLD)) JLOW = LOLD
      JUP = NORB(NSJ)
      if (NSI == NSJ) JUP = IOLD
      do JOLD=JLOW,JUP

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
        i1 = ind(ni)
        k1 = ind(nk)

        if ((i1 <= dimp) .and. (j1 <= dimq) .and. (k1 <= dimp) .and. (l1 <= dimq)) Vic(pqind(i1,k1),pqind(j1,l1)) = val1

        if (yes234 == 1) then

          !:2 combination (ij|kl) -> (ji|kl)
          ind(1) = jold
          ind(2) = iold
          ind(3) = kold
          ind(4) = lold
          j1 = ind(nj)
          l1 = ind(nl)
          i1 = ind(ni)
          k1 = ind(nk)

          if ((i1 <= dimp) .and. (j1 <= dimq) .and. (k1 <= dimp) .and. (l1 <= dimq)) Vic(pqind(i1,k1),pqind(j1,l1)) = val1

          !:3 combination (ij|kl) -> (ij|lk)
          ind(1) = iold
          ind(2) = jold
          ind(3) = lold
          ind(4) = kold
          j1 = ind(nj)
          l1 = ind(nl)
          i1 = ind(ni)
          k1 = ind(nk)

          if ((i1 <= dimp) .and. (j1 <= dimq) .and. (k1 <= dimp) .and. (l1 <= dimq)) Vic(pqind(i1,k1),pqind(j1,l1)) = val1

          !:4 combination (ij|kl) -> (ji|lk)
          ind(1) = jold
          ind(2) = iold
          ind(3) = lold
          ind(4) = kold
          j1 = ind(nj)
          l1 = ind(nl)
          i1 = ind(ni)
          k1 = ind(nk)

          if ((i1 <= dimp) .and. (j1 <= dimq) .and. (k1 <= dimp) .and. (l1 <= dimq)) Vic(pqind(i1,k1),pqind(j1,l1)) = val1

        end if

        !:5 combination (ij|kl) -> (kl|ij)
        if (yes5 == 1) then
          ind(1) = kold
          ind(2) = lold
          ind(3) = iold
          ind(4) = jold
          j1 = ind(nj)
          l1 = ind(nl)
          i1 = ind(ni)
          k1 = ind(nk)

          if ((i1 <= dimp) .and. (j1 <= dimq) .and. (k1 <= dimp) .and. (l1 <= dimq)) Vic(pqind(i1,k1),pqind(j1,l1)) = val1

        end if

        if (yes678 == 1) then

          !:6 combination (ij|kl) -> (lk|ij)
          ind(1) = lold
          ind(2) = kold
          ind(3) = iold
          ind(4) = jold
          j1 = ind(nj)
          l1 = ind(nl)
          i1 = ind(ni)
          k1 = ind(nk)

          if ((i1 <= dimp) .and. (j1 <= dimq) .and. (k1 <= dimp) .and. (l1 <= dimq)) Vic(pqind(i1,k1),pqind(j1,l1)) = val1

          !:7 combination (ij|kl) -> (kl|ji)
          ind(1) = kold
          ind(2) = lold
          ind(3) = jold
          ind(4) = iold
          j1 = ind(nj)
          l1 = ind(nl)
          i1 = ind(ni)
          k1 = ind(nk)

          if ((i1 <= dimp) .and. (j1 <= dimq) .and. (k1 <= dimp) .and. (l1 <= dimq)) Vic(pqind(i1,k1),pqind(j1,l1)) = val1

          !:8 combination (ij|kl) -> (lk|ji)
          ind(1) = lold
          ind(2) = kold
          ind(3) = jold
          ind(4) = iold
          j1 = ind(nj)
          l1 = ind(nl)
          i1 = ind(ni)
          k1 = ind(nk)

          if ((i1 <= dimp) .and. (j1 <= dimq) .and. (k1 <= dimp) .and. (l1 <= dimq)) Vic(pqind(i1,k1),pqind(j1,l1)) = val1

        end if

        indtemp = indtemp+1

      end do
    end do
  end do
end do

call mma_deallocate(TWO)

return

end subroutine esb_ic_2
