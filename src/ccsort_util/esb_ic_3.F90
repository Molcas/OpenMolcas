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

subroutine esb_ic_3(symp,Vic,dimp,pqind)
! this routine realizes expansion of symmetry block
! symp,symq,symr,syms <p,q|r,s>  <-> (IJ|KL),
! (for case symp=symr=symq=syms)
! provided such integrals exist
! It finds corresponding (IJ|KL) and expands it to
! matrix vic (prqs)

use ccsort_global, only: fullprint, idis, LUINTM, mbas, NORB, np, nq, nr, ns
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: symp, dimp
real(kind=wp), intent(_OUT_) :: Vic((dimp*(dimp+1)/2)*((dimp*(dimp+1)/2)+1)/2)
integer(kind=iwp), intent(out) :: pqind(mbas,mbas)
#include "tratoc.fh"
integer(kind=iwp) :: i, i1, idis13, ik, ikjl, ilow, ind(4), indtemp, iold, iup, j, j1, jl, jlow, jold, jup, k1, kold, kup, l1, &
                     lold, lup, maxx, ni, nj, nk, nl, nsi, nsj, nsk, nsl
real(kind=wp) :: val1
real(kind=wp), allocatable :: TWO(:)

!I calc pqind

!LD if (dimp >= dimp) then
maxx = dimp
!LD else
!LD   maxx = dimq
!LD end if

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
idis13 = idis(symp,symp,symp)

!III.1 define order of indices

ni = np(symp,symp,symp)
nj = nq(symp,symp,symp)
nk = nr(symp,symp,symp)
nl = ns(symp,symp,symp)

! define NSI,NSJ,NSK,NSL
ind(ni) = symp
ind(nj) = symp
ind(nk) = symp
ind(nl) = symp
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
        i1 = ind(ni)
        k1 = ind(nk)

        !:2 def iklj
        ik = pqind(i1,k1)
        jl = pqind(j1,l1)
        if (ik >= jl) then
          ikjl = ik*(ik-1)/2+jl
        else
          ikjl = jl*(jl-1)/2+ik
        end if

        Vic(ikjl) = val1

        indtemp = indtemp+1

      end do
    end do
  end do
end do

call mma_deallocate(TWO)

return

end subroutine esb_ic_3
