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

subroutine exchctl(exch,nneq,neqv,neq,nexch,nmax,lmax,npair,i_pair,MxRank1,MxRank2,imaxrank,Jex,JAex,JAex9,JDMex,JITOexR,JITOexI, &
                   eso,SM,MM,coord,rot,rlg,riso,tpar,upar,lant,itype,Dipol,AnisoLines1,AnisoLines3,AnisoLines9,KE,KEOPT, &
                   DM_exchange,JITO_exchange,W,Z,S,M,iPrint,mem)
! this Subroutine is a control Subroutine for the exchange interaction,
! diagonalization of total hamiltonian and computation of matrix elements
! of magnetic and spin moment
! neqv : max of neq(nneq)
! lant : (takes values from 1-7 for Gd-Yb respectively)
! mem  : memory allocated so far

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, cZero, cOne
use Definitions, only: wp, iwp, u6, CtoB, ItoB, RtoB

implicit none
integer(kind=iwp), intent(in) :: exch, nneq, neqv, neq(nneq), nexch(nneq), nmax, lmax, npair, i_pair(npair,2), MxRank1, MxRank2, &
                                 imaxrank(npair,2), lant, KEOPT, iPrint, mem
real(kind=wp), intent(in) :: Jex(npair), JAex(npair,3), JAex9(npair,3,3), JDMex(npair,3), &
                             JITOexR(nPair,MxRank1,-MxRank1:MxRank1,MxRank2,-MxRank2:MxRank2), &
                             JITOexI(nPair,MxRank1,-MxRank1:MxRank1,MxRank2,-MxRank2:MxRank2), eso(nneq,nmax), coord(nneq,3), &
                             rot(nneq,neqv,3,3), rlg(nneq,neqv,3,3), riso(3,3,nneq), tpar, upar
complex(kind=wp), intent(inout) :: SM(nneq,3,nmax,nmax), MM(nneq,3,nmax,nmax)
character, intent(in) :: itype(nneq)
logical(kind=iwp), intent(in) :: Dipol, AnisoLines1, AnisoLines3, AnisoLines9, KE, DM_exchange, JITO_exchange
real(kind=wp), intent(out) :: W(exch)
complex(kind=wp), intent(out) :: Z(exch,exch), S(3,exch,exch), M(3,exch,exch)
integer(kind=iwp) :: i, i1, i2, is1, isite, j, j1, j2, js1, l, lb1, lb2, lp, mem_local, n1, n2, nb, nb1, nb2, NmaxPop, nmaxR, nsta
real(kind=wp) :: dist, JAex_tmp(3), mg1(3,3), mg2(3,3), p1(3), p2(3), r1(3,3), r2(3,3), vect(3)
integer(kind=iwp), allocatable :: ibas(:,:), ibasR(:,:), icoord(:), intc(:), intcR(:), nexchR(:), nind(:,:)
real(kind=wp), allocatable :: eso_tmp(:,:), rotR(:,:,:,:), rotR_tmp(:,:,:,:), wdip(:), wdmo(:), wito(:), wkex(:), wlin(:), &
                              wlin1(:), wlin3(:), wlin9(:), WR(:)
complex(kind=wp), allocatable :: HDIP(:,:,:,:,:), HDIP_tmp(:,:,:,:,:), HDMO(:,:,:,:,:), HDMO_tmp(:,:,:,:,:), HITO(:,:,:,:,:), &
                                 HITO_tmp(:,:,:,:,:), HKEX(:,:,:,:,:), HKEXR(:,:,:,:,:), HKEXR_tmp(:,:,:,:,:), HLIN1(:,:,:,:,:), &
                                 HLIN1_tmp(:,:,:,:,:), HLIN3(:,:,:,:,:), HLIN3_tmp(:,:,:,:,:), HLIN9(:,:,:,:,:), &
                                 HLIN9_tmp(:,:,:,:,:), HTMP(:,:,:,:), M1(:,:,:), M2(:,:,:), MM_tmp1(:,:,:), MM_tmp2(:,:,:), &
                                 MM1(:,:,:), MM2(:,:,:), MMR(:,:,:,:), MR(:,:,:), S1(:,:,:), S2(:,:,:), SM_tmp1(:,:,:), &
                                 SM_tmp2(:,:,:), SM1(:,:,:), SM2(:,:,:), SMR(:,:,:,:), SR(:,:,:), tmp(:,:), tmp2(:,:), &
                                 tmp3(:,:,:,:,:), ZR(:,:)
integer(kind=iwp), parameter :: exchR = 8
#ifdef _DEBUGPRINT_
#  define _DBG_ .true.
integer(kind=iwp) :: k, k1, k2, q1, q2
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: dbg = _DBG_
integer(kind=iwp), external :: norder
real(kind=wp), external :: dnrm2_

#ifndef _DEBUGPRINT_
#include "macros.fh"
unused_var(riso)
unused_var(mem)
#endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef _DEBUGPRINT_
write(u6,'(A)') 'Enter EXCHCTL'
write(u6,'(A,  i8)') 'exch   = ',exch
write(u6,'(A,  i8)') 'nneq   = ',nneq
write(u6,'(A,  i8)') 'neqv   = ',neqv
write(u6,'(A,  i8)') 'nmax   = ',nmax
write(u6,'(A,  i8)') 'lmax   = ',lmax
write(u6,'(A,  i8)') 'nPair  = ',nPair
write(u6,'(A,  i8)') 'nPair  = ',nPair
write(u6,'(A,  i8)') 'MxRank1= ',MxRank1
write(u6,'(A,  i8)') 'MxRank2= ',MxRank2
write(u6,'(A,10i4)') 'neq()  = ',(neq(i),i=1,nneq)
write(u6,'(A,10i4)') 'nexch()= ',(nexch(i),i=1,nneq)
write(u6,'(A,  L2)') 'AnisoLines1  = ',AnisoLines1
write(u6,'(A,  L2)') 'AnisoLines3  = ',AnisoLines3
write(u6,'(A,  L2)') 'AnisoLines9  = ',AnisoLines9
write(u6,'(A,  L2)') 'DM_exchange  = ',DM_exchange
write(u6,'(A,  L2)') 'Dipol        = ',Dipol
write(u6,'(A,  L2)') 'JITO_exchange= ',JITO_exchange
if (AnisoLines1) then
  do i=1,nPair
    write(u6,'(A,2I4,F10.5)') 'LIN1',i_pair(i,1),i_pair(i,2),Jex(i)
  end do
end if

if (AnisoLines3) then
  do i=1,nPair
    write(u6,'(A,2I4,3F10.5)') 'LIN3',i_pair(i,1),i_pair(i,2),(JAex(i,j),j=1,3)
  end do
end if

if (AnisoLines9) then
  do i=1,nPair
    write(u6,'(A,2I4,9F10.5)') 'LIN9',i_pair(i,1),i_pair(i,2),((JAex9(i,j,k),j=1,3),k=1,3)
  end do
end if

if (DM_exchange) then
  do i=1,nPair
    write(u6,'(A,2I4,3F10.5)') 'DMEX',i_pair(i,1),i_pair(i,2),(JDMex(i,j),j=1,3)
  end do
end if

if (Dipol) then
  write(u6,'(A)') 'COORD(i):'
  do i=1,nneq
    write(u6,'(i3,3F14.8)') i,(coord(i,j),j=1,3)
  end do
end if

if (JITO_exchange) then
  do i=1,nPair
    write(u6,'(A,4I4,3F10.5)') 'ITOJ',i_pair(i,1),i_pair(i,2),imaxrank(i,1),imaxrank(i,2)
    do k1=1,imaxrank(i,1),2
      do q1=-k1,k1
        do k2=1,imaxrank(i,2),2
          do q2=-k2,k2
            write(u6,'(4I3,2x,2ES21.14)') k1,q1,k2,q2,JITOexR(i,k1,q1,k2,q2),JITOexI(i,k1,q1,k2,q2)
          end do
        end do
      end do
    end do
  end do ! ipair
end if

write(u6,'(A)') 'ESO(i):'
do i=1,nneq
  write(u6,'(i3,90F14.8)') i,(eso(i,j),j=1,nmax)
end do
write(u6,'(90A)') 'itype()=',(itype(i),' ',i=1,nneq)

write(u6,'(A,  i8)') 'neqv   = ',neqv
do i=1,nneq
  write(u6,'(A,  i8)') 'riso( site=',i,'):'
  do j=1,3
    write(u6,'(3ES22.14)') (riso(j,k,i),k=1,3)
  end do
end do
#endif
!-----------------------------------------------------------------------
! allocate memory for this function:
mem_local = 0
if (lmax >= 0) then
  ! exchange energy spectrum
  call mma_allocate(intc,lmax,'intc')
  call mma_allocate(icoord,lmax,'icoord')
  call mma_allocate(nind,lmax,2,'nind')
  mem_local = mem_local+size(intc)*ItoB
  mem_local = mem_local+size(icoord)*ItoB
  mem_local = mem_local+size(nind)*ItoB
  if (exch >= 0) then
    call mma_allocate(ibas,exch,lmax,'ibas')
    mem_local = mem_local+size(ibas)*ItoB
  end if
end if
if (exch >= 0) then
  call mma_allocate(wlin,exch,'wlin ')
  call mma_allocate(wlin1,exch,'wlin1')
  call mma_allocate(wlin3,exch,'wlin3')
  call mma_allocate(wlin9,exch,'wlin9')
  call mma_allocate(wdip,exch,'wdip ')
  call mma_allocate(wkex,exch,'wkex ')
  call mma_allocate(wdmo,exch,'wdmo ')
  call mma_allocate(wito,exch,'wito ')
  mem_local = mem_local+size(wlin)*RtoB
  mem_local = mem_local+size(wlin1)*RtoB
  mem_local = mem_local+size(wlin3)*RtoB
  mem_local = mem_local+size(wlin9)*RtoB
  mem_local = mem_local+size(wdip)*RtoB
  mem_local = mem_local+size(wkex)*RtoB
  mem_local = mem_local+size(wdmo)*RtoB
  mem_local = mem_local+size(wito)*RtoB
end if
if (nmax >= 0) then
  call mma_allocate(S1,3,nmax,nmax,' S1')
  call mma_allocate(M1,3,nmax,nmax,' M1')
  call mma_allocate(S2,3,nmax,nmax,' S2')
  call mma_allocate(M2,3,nmax,nmax,' M2')
  call mma_allocate(SM1,3,nmax,nmax,'SM1')
  call mma_allocate(SM2,3,nmax,nmax,'SM2')
  call mma_allocate(MM1,3,nmax,nmax,'MM1')
  call mma_allocate(MM2,3,nmax,nmax,'MM2')
  mem_local = mem_local+size(S1)*CtoB
  mem_local = mem_local+size(M1)*CtoB
  mem_local = mem_local+size(S2)*CtoB
  mem_local = mem_local+size(M2)*CtoB
  mem_local = mem_local+size(SM1)*CtoB
  mem_local = mem_local+size(SM2)*CtoB
  mem_local = mem_local+size(MM1)*CtoB
  mem_local = mem_local+size(MM2)*CtoB

  if (npair >= 0) then
    call mma_allocate(HLIN1,npair,nmax,nmax,nmax,nmax,'HLIN1')
    call mma_allocate(HLIN3,npair,nmax,nmax,nmax,nmax,'HLIN3')
    call mma_allocate(HLIN9,npair,nmax,nmax,nmax,nmax,'HLIN9')
    call mma_allocate(HDIP,npair,nmax,nmax,nmax,nmax,'HDIP')
    call mma_allocate(HKEX,npair,nmax,nmax,nmax,nmax,'HKEX')
    call mma_allocate(HDMO,npair,nmax,nmax,nmax,nmax,'HDMO')
    call mma_allocate(HITO,npair,nmax,nmax,nmax,nmax,'HITO')
    HLIN1(:,:,:,:,:) = cZero
    HLIN3(:,:,:,:,:) = cZero
    HLIN9(:,:,:,:,:) = cZero
    HDIP(:,:,:,:,:) = cZero
    HKEX(:,:,:,:,:) = cZero
    HDMO(:,:,:,:,:) = cZero
    HITO(:,:,:,:,:) = cZero
    mem_local = mem_local+size(HLIN1)*CtoB
    mem_local = mem_local+size(HLIN3)*CtoB
    mem_local = mem_local+size(HLIN9)*CtoB
    mem_local = mem_local+size(HDIP)*CtoB
    mem_local = mem_local+size(HKEX)*CtoB
    mem_local = mem_local+size(HDMO)*CtoB
    mem_local = mem_local+size(HITO)*CtoB
  end if
end if

if (exch >= 0) then
  call mma_allocate(tmp,exch,exch,'tmp')
  mem_local = mem_local+size(tmp)*CtoB
end if

if (nneq >= 0) then
  call mma_allocate(nexchR,nneq,'nexchR')
  mem_local = mem_local+size(nexchR)*ItoB

  call mma_allocate(SMR,nneq,3,2,2,'SMR')
  call mma_allocate(MMR,nneq,3,2,2,'MMR')
  mem_local = mem_local+size(SMR)*CtoB
  mem_local = mem_local+size(MMR)*CtoB

  if (neqv >= 0) then
    call mma_allocate(rotR,nneq,neqv,3,3,'rotR')
    mem_local = mem_local+size(rotR)*RtoB
  end if
end if

if (exchR >= 0) then
  if (lmax >= 0) then
    call mma_allocate(ibasR,nneq,lmax,'ibasR')
    mem_local = mem_local+size(ibasR)*ItoB
  end if
  call mma_allocate(WR,exchR,'WR')
  call mma_allocate(ZR,exchR,exchR,'ZR')
  call mma_allocate(MR,3,exchR,exchR,'MR')
  call mma_allocate(SR,3,exchR,exchR,'SR')
  mem_local = mem_local+size(WR)*RtoB
  mem_local = mem_local+size(ZR)*CtoB
  mem_local = mem_local+size(MR)*CtoB
  mem_local = mem_local+size(SR)*CtoB
end if

if (npair >= 0) then
  call mma_allocate(HKEXR,npair,2,2,2,2,'HKEXR')
  mem_local = mem_local+size(HKEXR)*CtoB
end if

if (lmax >= 0) then
  call mma_allocate(intcR,lmax,'intcR')
  mem_local = mem_local+size(intcR)*ItoB
end if
#ifdef _DEBUGPRINT_
write(u6,*) 'EXCHCTL:  memory allocated (local):'
write(u6,*) 'mem_local=',mem_local
write(u6,*) 'EXCHCTL:  memory allocated (total):'
write(u6,*) 'mem_total=',mem+mem_local
#endif
!-----------------------------------------------------------------------
l = 0
do i=1,nneq
  do j=1,neq(i)
    l = l+1
    nind(l,1) = i
    nind(l,2) = j
  end do
end do
nind(l+1:,:) = 0
intc(1) = 1
if (lmax > 1) then
  do i=2,lmax
    isite = nind(i-1,1)
    intc(i) = intc(i-1)*nexch(isite)
  end do
end if
do nb=1,exch
  nb1 = nb-1
  do i=1,lmax
    ibas(nb,lmax-i+1) = nb1/intc(lmax-i+1)
    nb1 = nb1-ibas(nb,lmax-i+1)*intc(lmax-i+1)
  end do
end do
#ifdef _DEBUGPRINT_
write(u6,'(34x,A,1x,20i3)') 'site Nr.',(i,i=1,lmax)
do nb=1,exch
  write(u6,'(A,i5,A,20i3)') 'COUPLING: basis set:  ibas(',nb,' ,isite) = ',(ibas(nb,i)+1,i=1,lmax)
end do
#endif

!-----------------------------------------------------------------------
! Lines model of magnetic couping  -- 1 parameter
if (AnisoLines1) then
  if (nPair > 0) then
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      call mma_allocate(SM_tmp1,3,n1,n1,label='SM_tmp1')
      call mma_allocate(SM_tmp2,3,n2,n2,label='SM_tmp2')
      call mma_allocate(MM_tmp1,3,n1,n1,label='MM_tmp1')
      call mma_allocate(MM_tmp2,3,n2,n2,label='MM_tmp2')
      SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
      SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
      MM_tmp1(:,:,:) = MM(i1,:,1:n1,1:n1)
      MM_tmp2(:,:,:) = MM(i2,:,1:n2,1:n2)

      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        r1(:,:) = rot(i1,j1,:,:)
        call prep_mom_exchange(n1,r1,SM_tmp1,MM_tmp1,mg1,.true.)
        SM(i1,:,1:n1,1:n1) = SM_tmp1(:,:,:)
        MM(i1,:,1:n1,1:n1) = MM_tmp1(:,:,:)
      end if
      if (itype(i2) == 'A') then
        r2(:,:) = rot(i2,j2,:,:)
        call prep_mom_exchange(n2,r2,SM_tmp2,MM_tmp2,mg2,.true.)
        SM(i2,:,1:n2,1:n2) = SM_tmp2(:,:,:)
        MM(i2,:,1:n2,1:n2) = MM_tmp2(:,:,:)
      end if

      call prMom('SM(i1) bf Lines1',SM_tmp1,n1)
      call prMom('SM(i2) bf Lines1',SM_tmp2,n2)

      ! build the Lines exchange matrix:
      call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
      call Lines_Exchange(Jex(lp),n1,n2,SM_tmp1,SM_tmp2,HTMP)
      HLIN1(lp,1:n1,1:n1,1:n2,1:n2) = HTMP(:,:,:,:)
      call mma_deallocate(HTMP)

#     ifdef _DEBUGPRINT_
      call prMom('SM(i1) af Lines1',SM_tmp1,n1)
      call prMom('SM(i2) af Lines1',SM_tmp2,n2)
      write(u6,'(A,i2)') 'Exchange matrix for pair = ',lp
      write(u6,'(A,i2)') 'in local pseudospin basis'
      do i=1,n1
        do j=1,n1
          do k=1,n2
            do l=1,n2
              write(u6,'(4(a,i2),a,2ES24.14)') 'HLIN1(',i,',',j,',',k,',',l,') = ',HLIN1(lp,i,j,k,l)
            end do
          end do
        end do
      end do
#     endif

      call mma_deallocate(SM_tmp1)
      call mma_deallocate(SM_tmp2)
      call mma_deallocate(MM_tmp1)
      call mma_deallocate(MM_tmp2)

    end do
  end if ! nPair
end if
!-----------------------------------------------------------------------

! Anisotropic Lines model of magnetic couping -- 3 parameters)
! Jxx, Jyy, Jzz
if (AnisoLines3) then
  if (nPair > 0) then
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      call mma_allocate(SM_tmp1,3,n1,n1,label='SM_tmp1')
      call mma_allocate(SM_tmp2,3,n2,n2,label='SM_tmp2')
      call mma_allocate(MM_tmp1,3,n1,n1,label='MM_tmp1')
      call mma_allocate(MM_tmp2,3,n2,n2,label='MM_tmp2')
      SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
      SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
      MM_tmp1(:,:,:) = MM(i1,:,1:n1,1:n1)
      MM_tmp2(:,:,:) = MM(i2,:,1:n2,1:n2)

      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        r1(:,:) = rot(i1,j1,:,:)
        call prep_mom_exchange(n1,r1,SM_tmp1,MM_tmp1,mg1,dbg)
        SM(i1,:,1:n1,1:n1) = SM_tmp1(:,:,:)
        MM(i1,:,1:n1,1:n1) = MM_tmp1(:,:,:)
      end if
      if (itype(i2) == 'A') then
        r2(:,:) = rot(i2,j2,:,:)
        call prep_mom_exchange(n2,r2,SM_tmp2,MM_tmp2,mg2,dbg)
        SM(i2,:,1:n2,1:n2) = SM_tmp2(:,:,:)
        MM(i2,:,1:n2,1:n2) = MM_tmp2(:,:,:)
      end if

      call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
      JAex_tmp(:) = JAex(lp,:)
      call Aniso_Lines_Exchange3(JAex_tmp,n1,n2,SM_tmp1,SM_tmp2,HTMP)
      HLIN3(lp,1:n1,1:n1,1:n2,1:n2) = HTMP(:,:,:,:)
      call mma_deallocate(HTMP)

#     ifdef _DEBUGPRINT_
      write(u6,'(A,i2)') 'Exchange matrix for pair = ',lp
      do i=1,n1
        do j=1,n1
          do k=1,n2
            do l=1,n2
              write(u6,'(4(a,i2),a,2ES24.14)') 'HLIN3(',i,',',j,',',k,',',l,') = ',HLIN3(lp,i,j,k,l)
            end do
          end do
        end do
      end do
#     endif

      call mma_deallocate(SM_tmp1)
      call mma_deallocate(SM_tmp2)
      call mma_deallocate(MM_tmp1)
      call mma_deallocate(MM_tmp2)

    end do ! lp
  end if ! nPair
end if
!-----------------------------------------------------------------------

! Anisotropic Lines model of magnetic couping -- 9 parameters)
! Jxx, Jxy, Jxz, Jyx, Jyy, Jyz, Jzx, Jzy, Jzz
if (AnisoLines9) then
  if (nPair > 0) then
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      call mma_allocate(SM_tmp1,3,n1,n1,label='SM_tmp1')
      call mma_allocate(SM_tmp2,3,n2,n2,label='SM_tmp2')
      call mma_allocate(MM_tmp1,3,n1,n1,label='MM_tmp1')
      call mma_allocate(MM_tmp2,3,n2,n2,label='MM_tmp2')
      SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
      SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
      MM_tmp1(:,:,:) = MM(i1,:,1:n1,1:n1)
      MM_tmp2(:,:,:) = MM(i2,:,1:n2,1:n2)

      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        r1(:,:) = rot(i1,j1,:,:)
        call prep_mom_exchange(n1,r1,SM_tmp1,MM_tmp1,mg1,dbg)
        SM(i1,:,1:n1,1:n1) = SM_tmp1(:,:,:)
        MM(i1,:,1:n1,1:n1) = MM_tmp1(:,:,:)
      end if
      if (itype(i2) == 'A') then
        r2(:,:) = rot(i2,j2,:,:)
        call prep_mom_exchange(n2,r2,SM_tmp2,MM_tmp2,mg2,dbg)
        SM(i2,:,1:n2,1:n2) = SM_tmp2(:,:,:)
        MM(i2,:,1:n2,1:n2) = MM_tmp2(:,:,:)
      end if

      ! rotate the input J matrix by the mg1 and mg2
#     ifdef _DEBUGPRINT_
      call prMom('SM(i1) bf Lines9',SM_tmp1,n1)
      call prMom('SM(i2) bf Lines9',SM_tmp2,n2)
#     endif

      call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
      call Aniso_Lines_Exchange9(JAex9(lp,:,:),n1,n2,SM_tmp1,SM_tmp2,HTMP)
      HLIN9(lp,1:n1,1:n1,1:n2,1:n2) = HTMP(:,:,:,:)
      call mma_deallocate(HTMP)

#     ifdef _DEBUGPRINT_
      call prMom('SM(i1) af Lines9',SM_tmp1,n1)
      call prMom('SM(i2) af Lines9',SM_tmp2,n2)
      write(u6,'(A,i2)') 'Exchange matrix for pair = ',lp
      do i=1,n1
        do j=1,n1
          do k=1,n2
            do l=1,n2
              write(u6,'(4(a,i2),a,2ES24.14)') 'HLIN9(',i,',',j,',',k,',',l,') = ',HLIN9(lp,i,j,k,l)
            end do
          end do
        end do
      end do
#     endif

      call mma_deallocate(SM_tmp1)
      call mma_deallocate(SM_tmp2)
      call mma_deallocate(MM_tmp1)
      call mma_deallocate(MM_tmp2)

    end do
  end if ! nPair
end if

!-----------------------------------------------------------------------
! dipolar couping
if (Dipol) then
  if (nPair > 0) then
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      p1(:) = coord(i1,:)
      p2(:) = coord(i2,:)
      r1(:,:) = rlg(i1,j1,:,:)
      r2(:,:) = rlg(i2,j2,:,:)
      call dirvect(p1,r1,p2,r2,vect,dist)
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i3,3ES20.12,2x,ES20.12)') 'EXCHCTL: DIPOL: lp, vect, R:',lp,(vect(i),i=1,3),dist
#     endif

      call mma_allocate(SM_tmp1,3,n1,n1,label='SM_tmp1')
      call mma_allocate(SM_tmp2,3,n2,n2,label='SM_tmp2')
      call mma_allocate(MM_tmp1,3,n1,n1,label='MM_tmp1')
      call mma_allocate(MM_tmp2,3,n2,n2,label='MM_tmp2')
      SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
      SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
      MM_tmp1(:,:,:) = MM(i1,:,1:n1,1:n1)
      MM_tmp2(:,:,:) = MM(i2,:,1:n2,1:n2)

      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        r1(:,:) = rot(i1,j1,:,:)
        call prep_mom_exchange(n1,r1,SM_tmp1,MM_tmp1,mg1,dbg)
        SM(i1,:,1:n1,1:n1) = SM_tmp1(:,:,:)
        MM(i1,:,1:n1,1:n1) = MM_tmp1(:,:,:)
      end if
      if (itype(i2) == 'A') then
        r2(:,:) = rot(i2,j2,:,:)
        call prep_mom_exchange(n2,r2,SM_tmp2,MM_tmp2,mg2,dbg)
        SM(i2,:,1:n2,1:n2) = SM_tmp2(:,:,:)
        MM(i2,:,1:n2,1:n2) = MM_tmp2(:,:,:)
      end if

      call mma_allocate(HTMP,n1,n1,n2,n2,label='HTMP')
      call Dipol_Exchange(n1,n2,vect,dist,MM_tmp1,MM_tmp2,HTMP)
      HDIP(lp,1:n1,1:n1,1:n2,1:n2) = HTMP(:,:,:,:)
      call mma_deallocate(HTMP)

#     ifdef _DEBUGPRINT_
      write(u6,'(A,i2)') 'Exchange matrix for pair = ',lp
      do i=1,n1
        do j=1,n1
          do k=1,n2
            do l=1,n2
              write(u6,'(4(a,i2),a,2ES24.14)') 'HDIP (',i,',',j,',',k,',',l,') = ',HDIP(lp,i,j,k,l)
            end do
          end do
        end do
      end do
#     endif

      call mma_deallocate(SM_tmp1)
      call mma_deallocate(SM_tmp2)
      call mma_deallocate(MM_tmp1)
      call mma_deallocate(MM_tmp2)

    end do ! lp
  end if
end if
!-----------------------------------------------------------------------

! Dzyaloshinsky-Morya antisymmetric couping
if (DM_exchange) then
  if (nPair > 0) then
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      call mma_allocate(SM_tmp1,3,n1,n1,label='SM_tmp1')
      call mma_allocate(SM_tmp2,3,n2,n2,label='SM_tmp2')
      call mma_allocate(MM_tmp1,3,n1,n1,label='MM_tmp1')
      call mma_allocate(MM_tmp2,3,n2,n2,label='MM_tmp2')
      SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
      SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
      MM_tmp1(:,:,:) = MM(i1,:,1:n1,1:n1)
      MM_tmp2(:,:,:) = MM(i2,:,1:n2,1:n2)

      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        r1(:,:) = rot(i1,j1,:,:)
        call prep_mom_exchange(n1,r1,SM_tmp1,MM_tmp1,mg1,dbg)
        SM(i1,:,1:n1,1:n1) = SM_tmp1(:,:,:)
        MM(i1,:,1:n1,1:n1) = MM_tmp1(:,:,:)
      end if
      if (itype(i2) == 'A') then
        r2(:,:) = rot(i2,j2,:,:)
        call prep_mom_exchange(n2,r2,SM_tmp2,MM_tmp2,mg2,dbg)
        SM(i2,:,1:n2,1:n2) = SM_tmp2(:,:,:)
        MM(i2,:,1:n2,1:n2) = MM_tmp2(:,:,:)
      end if

      call Dzyaloshinsky_Morya_Exchange(JDMex(lp,:),n1,n2,SM_tmp1,SM_tmp2,HDMO(lp,1:n1,1:n1,1:n2,1:n2))
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i2)') 'Exchange matrix for pair = ',lp
      do i=1,n1
        do j=1,n1
          do k=1,n2
            do l=1,n2
              write(u6,'(4(a,i2),a,2ES24.14)') 'HDMO (',i,',',j,',',k,',',l,') = ',HDMO(lp,i,j,k,l)
            end do
          end do
        end do
      end do
#     endif

      call mma_deallocate(SM_tmp1)
      call mma_deallocate(SM_tmp2)
      call mma_deallocate(MM_tmp1)
      call mma_deallocate(MM_tmp2)

    end do
  end if ! nPair
end if
!-----------------------------------------------------------------------

! JITO exchange interaction
if (JITO_exchange) then
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'EXCHCTL:  Entering  JITO_exchange'
# endif
  if (nPair > 0) then
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      call mma_allocate(SM_tmp1,3,n1,n1,label='SM_tmp1')
      call mma_allocate(SM_tmp2,3,n2,n2,label='SM_tmp2')
      call mma_allocate(MM_tmp1,3,n1,n1,label='MM_tmp1')
      call mma_allocate(MM_tmp2,3,n2,n2,label='MM_tmp2')
      SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
      SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
      MM_tmp1(:,:,:) = MM(i1,:,1:n1,1:n1)
      MM_tmp2(:,:,:) = MM(i2,:,1:n2,1:n2)

      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        r1(:,:) = rot(i1,j1,:,:)
        call prep_mom_exchange(n1,r1,SM_tmp1,MM_tmp1,mg1,dbg)
        SM(i1,:,1:n1,1:n1) = SM_tmp1(:,:,:)
        MM(i1,:,1:n1,1:n1) = MM_tmp1(:,:,:)
      end if
      if (itype(i2) == 'A') then
        r2(:,:) = rot(i2,j2,:,:)
        call prep_mom_exchange(n2,r2,SM_tmp2,MM_tmp2,mg2,dbg)
        SM(i2,:,1:n2,1:n2) = SM_tmp2(:,:,:)
        MM(i2,:,1:n2,1:n2) = MM_tmp2(:,:,:)
      end if

      ! using Naoya's ITO:, in general coordinate system
      call JITO_Exchange_Int(MxRank1,MxRank2,imaxrank(lp,:),n1,n2,JITOexR(lp,:,:,:,:),JITOexI(lp,:,:,:,:), &
                             HITO(lp,1:n1,1:n1,1:n2,1:n2))
#     ifdef _DEBUGPRINT_
      write(u6,'(A,i2)') 'Exchange matrix for pair = ',lp
      do i=1,n1
        do j=1,n1
          do k=1,n2
            do l=1,n2
              write(u6,'(4(a,i2),a,2ES24.14)') 'HITO (',i,',',j,',',k,',',l,') = ',HITO(lp,i,j,k,l)
            end do
          end do
        end do
      end do
#     endif

      call mma_deallocate(SM_tmp1)
      call mma_deallocate(SM_tmp2)
      call mma_deallocate(MM_tmp1)
      call mma_deallocate(MM_tmp2)

    end do
  end if ! nPair
# ifdef _DEBUGPRINT_
  write(u6,'(A)') 'EXCHCTL:  Exiting JITO_exchange'
# endif
end if
!-----------------------------------------------------------------------

if (KE) then
  if (nPair > 0) then
    HKEX(:,:,:,:,:) = cZero
    HKEXR(:,:,:,:,:) = cZero
    MMR(:,:,:,:) = cZero
    SMR(:,:,:,:) = cZero
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      call mma_allocate(SM_tmp1,3,n1,n1,label='SM_tmp1')
      call mma_allocate(SM_tmp2,3,n2,n2,label='SM_tmp2')
      call mma_allocate(MM_tmp1,3,n1,n1,label='MM_tmp1')
      call mma_allocate(MM_tmp2,3,n2,n2,label='MM_tmp2')
      SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
      SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
      MM_tmp1(:,:,:) = MM(i1,:,1:n1,1:n1)
      MM_tmp2(:,:,:) = MM(i2,:,1:n2,1:n2)

      r1(:,:) = rot(i1,j1,:,:)
      call mma_allocate(tmp3,3,n1,n1,1,1,label='tmp3')
      call rotmom2(MM_tmp1,n1,r1,tmp3)
      M1(:,1:n1,1:n1) = tmp3(:,:,:,1,1)
      call rotmom2(SM_tmp1,n1,r1,tmp3)
      S1(:,1:n1,1:n1) = tmp3(:,:,:,1,1)
      call mma_deallocate(tmp3)
      r2(:,:) = rot(i2,j2,:,:)
      call mma_allocate(tmp3,3,n2,n2,1,1,label='tmp3')
      call rotmom2(MM_tmp2,n2,r2,tmp3)
      M2(:,1:n2,1:n2) = tmp3(:,:,:,1,1)
      call rotmom2(SM_tmp2,n2,r2,tmp3)
      S2(:,1:n2,1:n2) = tmp3(:,:,:,1,1)
      call mma_deallocate(tmp3)

      call mma_deallocate(SM_tmp1)
      call mma_deallocate(SM_tmp2)
      call mma_deallocate(MM_tmp1)
      call mma_deallocate(MM_tmp2)

      ! KEOPT=1 ! FULL
      ! KEOPT=2 ! Full, 1/U
      ! KEOPT=3 ! FULL  + reduced form
      ! KEOPT=4 ! Full, 1/U + reduced form
      MM1(:,:,:) = cZero
      SM1(:,:,:) = cZero
      MM2(:,:,:) = cZero
      SM2(:,:,:) = cZero
      !FIXME: this call to KE_Exchange does not match at all its definition, please fix (and avoid passing non-contiguous arrays)
      call WarningMessage(2,'There is surely a bug here')
#     include "macros.fh"
      unused_var(tpar)
      unused_var(upar)
      unused_var(lant)
      !call KE_Exchange(n1,n2,M1(:,1:n1,1:n1),S1(:,1:n1,1:n1),M2(:,1:n2,1:n2),S2(:,1:n2,1:n2),eso(i1,1:n1),eso(i2,1:n2),tpar,upar, &
      !                 lant,KEOPT,HKEX(lp,1:n1,1:n1,1:n2,1:n2),MM1(:,1:n1,1:n1),SM1(:,1:n1,1:n1),MM2(:,1:n2,1:n2), &
      !                 SM2(1:3,1:n2,1:n2)

      if ((KEOPT == 1) .or. (KEOPT == 2)) then
        HKEXR(lp,1,1,:,:) = HKEX(lp,1,1,1:n2,1:n2)
        HKEXR(lp,1,2,:,:) = HKEX(lp,1,2,1:n2,1:n2)
        HKEXR(lp,2,1,:,:) = HKEX(lp,2,1,1:n2,1:n2)
        HKEXR(lp,2,2,:,:) = HKEX(lp,2,2,1:n2,1:n2)
        MMR(i1,:,:,:) = MM1(:,1:2,1:2)
        MMR(i2,:,:,:) = MM2(:,1:2,1:2)
        SMR(i1,:,:,:) = SM1(:,1:2,1:2)
        SMR(i2,:,:,:) = SM2(:,1:2,1:2)

      else if ((KEOPT == 3) .or. (KEOPT == 4)) then
        HKEXR(lp,1,1,:,:) = HKEX(lp,1,1,1:n2,1:n2)
        HKEXR(lp,1,2,:,:) = HKEX(lp,1,n1,1:n2,1:n2)
        HKEXR(lp,2,1,:,:) = HKEX(lp,n1,1,1:n2,1:n2)
        HKEXR(lp,2,2,:,:) = HKEX(lp,n1,n1,1:n2,1:n2)

        MMR(i1,:,1,1) = MM1(:,1,1)
        MMR(i1,:,1,2) = MM1(:,1,n1)
        MMR(i1,:,2,1) = MM1(:,n1,1)
        MMR(i1,:,2,2) = MM1(:,n1,n1)

        MMR(i2,:,1,1) = MM2(:,1,1)
        MMR(i2,:,1,2) = MM2(:,1,n2)
        MMR(i2,:,2,1) = MM2(:,n2,1)
        MMR(i2,:,2,2) = MM2(:,n2,n2)

        SMR(i1,:,1,1) = SM1(:,1,1)
        SMR(i1,:,1,2) = SM1(:,1,n1)
        SMR(i1,:,2,1) = SM1(:,n1,1)
        SMR(i1,:,2,2) = SM1(:,n1,n1)

        SMR(i2,:,1,1) = SM2(:,1,1)
        SMR(i2,:,1,2) = SM2(:,1,n2)
        SMR(i2,:,2,1) = SM2(:,n2,1)
        SMR(i2,:,2,2) = SM2(:,n2,n2)

#       ifdef _DEBUGPRINT_
        write(u6,'(A)') 'site Ln'
        do i=1,2
          do j=1,2
            write(u6,'(3(A,i1),A,3(2F20.14,2x))') 'MMR(',i1,',L,',i,',',j,')= ',(MMR(i1,l,i,j),l=1,3)
          end do
        end do
        write(u6,'(/)')
        do i=1,2
          do j=1,2
            write(u6,'(3(A,i1),A,3(2F20.14,2x))') 'SMR(',i1,',L,',i,',',j,')= ',(SMR(i1,l,i,j),l=1,3)
          end do
        end do
        write(u6,'(/)')
        write(u6,'(A)') 'site Radical'
        do i=1,2
          do j=1,2
            write(u6,'(3(A,i1),A,3(2F20.14,2x))') 'MMR(',i2,',L,',i,',',j,')= ',(MMR(i2,l,i,j),l=1,3)
          end do
        end do
        write(u6,'(/)')
        do i=1,2
          do j=1,2
            write(u6,'(3(A,i1),A,3(2F20.14,2x))') 'SMR(',i2,',L,',i,',',j,')= ',(SMR(i2,l,i,j),l=1,3)
          end do
        end do
#       endif

      end if ! KEOPT

    end do ! lp

    ! in case of KEOPT=3 or KEOPT=4 Then we need to compute the spectrum and the properties
    ! in the reduced form, where nexch(i1)=2 ( ground doublet only).
    ! exchnew=8:
    if (KEOPT <= 4) then
      nmaxR = 2
      nexchR(:) = 2
      ibasR(:,:) = 0
      intcR(1) = 1
      if (lmax > 1) then
        do i=2,lmax
          isite = nind(i-1,1)
          intcR(i) = intcR(i-1)*nexchR(isite)
        end do
      end if
      do nb=1,exchR
        nb1 = nb-1
        do i=1,lmax
          ibasR(nb,lmax-i+1) = nb1/intcR(lmax-i+1)
          nb1 = nb1-ibasR(nb,lmax-i+1)*intcR(lmax-i+1)
        end do
      end do
      !-----------------------------------------------------------------
      HLIN1(:,:,:,:,:) = cZero
      HLIN3(:,:,:,:,:) = cZero
      HLIN9(:,:,:,:,:) = cZero
      HDIP(:,:,:,:,:) = cZero
      WLIN(:) = Zero
      WLIN1(:) = Zero
      WLIN3(:) = Zero
      WLIN9(:) = Zero
      WDIP(:) = Zero
      WKEX(:) = Zero
      call mma_allocate(eso_tmp,nneq,nmaxR,label='eso_tmp')
      call mma_allocate(HLIN1_tmp,npair,nmaxR,nmaxR,nmaxR,nmaxR,label='HLIN1_tmp')
      call mma_allocate(HLIN3_tmp,npair,nmaxR,nmaxR,nmaxR,nmaxR,label='HLIN3_tmp')
      call mma_allocate(HLIN9_tmp,npair,nmaxR,nmaxR,nmaxR,nmaxR,label='HLIN9_tmp')
      call mma_allocate(HDIP_tmp,npair,nmaxR,nmaxR,nmaxR,nmaxR,label='HDIP_tmp')
      call mma_allocate(HKEXR_tmp,npair,nmaxR,nmaxR,nmaxR,nmaxR,label='HKEXR_tmp')
      call mma_allocate(HDMO_tmp,npair,nmaxR,nmaxR,nmaxR,nmaxR,label='HDMO_tmp')
      call mma_allocate(HITO_tmp,npair,nmaxR,nmaxR,nmaxR,nmaxR,label='HITO_tmp')
      eso_tmp(:,:) = eso(:,1:nmaxR)
      HLIN1_tmp(:,:,:,:,:) = HLIN1(:,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR)
      HLIN3_tmp(:,:,:,:,:) = HLIN3(:,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR)
      HLIN9_tmp(:,:,:,:,:) = HLIN9(:,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR)
      HDIP_tmp(:,:,:,:,:) = HDIP(:,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR)
      HKEXR_tmp(:,:,:,:,:) = HKEXR(:,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR)
      HDMO_tmp(:,:,:,:,:) = HDMO(:,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR)
      HITO_tmp(:,:,:,:,:) = HITO(:,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR)
      ! print the Exchange Hamiltonian:
      call pa_prham(exchR,npair,i_pair,nneq,neq,nexchR,nmaxR,lmax,eso_tmp,HLIN1_tmp,HLIN3_tmp,HLIN9_tmp,HDIP_tmp,HKEXR_tmp, &
                    HDMO_tmp,HITO_tmp,Dipol,AnisoLines1,AnisoLines3,AnisoLines9,KE,DM_exchange,.false.)
      ! diagonalize the Hamiltonian:
      call pa_diagham(exchR,npair,i_pair,nneq,neq,nexchR,nmaxR,lmax,eso_tmp,HLIN1_tmp,HLIN3_tmp,HLIN9_tmp,HDIP_tmp,HKEXR_tmp, &
                      HDMO_tmp,HITO_tmp,Dipol,.false.,AnisoLines1,AnisoLines3,AnisoLines9,KE,.false.,WLIN1(1:exchR), &
                      WLIN3(1:exchR),WLIN9(1:exchR),WLIN(1:exchR),WDIP(1:exchR),WKEX(1:exchR),WDMO(1:exchR),WITO(1:exchR), &
                      WR,ZR)
      ! print the resulting eigenstates:
      call pa_preigen(exchR,lmax,ibasR,Dipol,AnisoLines1,AnisoLines3,AnisoLines9,KE,.false.,WLIN(1:exchR),WDIP(1:exchR), &
                      WKEX(1:exchR),WITO(1:exchR),WR,ZR,0)
      ! compute the moments:
      rotR(:,:,:,:) = Zero
      rotR(1,1,1,1) = One
      rotR(1,1,2,2) = One
      rotR(1,1,3,3) = One
      rotR(1,2,1,1) = One
      rotR(1,2,2,2) = One
      rotR(1,2,3,3) = One
      rotR(2,1,1,1) = One
      rotR(2,1,2,2) = One
      rotR(2,1,3,3) = One
      MR(:,:,:) = cZero
      SR(:,:,:) = cZero
      call mma_allocate(tmp2,EXCHR,EXCHR,label='tmp2')
      do L=1,3
        do isite=1,lmax
          do nb1=1,exchR
            icoord(:) = ibasR(nb1,:)
            i1 = nind(isite,1)
            j1 = nind(isite,2)
            is1 = ibasR(nb1,isite)+1

            do js1=1,nexchR(i1)
              icoord(isite) = js1-1
              nb2 = norder(icoord,intcR,lmax)
              MR(l,nb1,nb2) = MR(l,nb1,nb2)+rotR(i1,j1,l,1)*MMR(i1,1,is1,js1)+rotR(i1,j1,l,2)*MMR(i1,2,is1,js1)+ &
                              rotR(i1,j1,l,3)*MMR(i1,3,is1,js1)
              SR(l,nb1,nb2) = SR(l,nb1,nb2)+rotR(i1,j1,l,1)*SMR(i1,1,is1,js1)+rotR(i1,j1,l,2)*SMR(i1,2,is1,js1)+ &
                              rotR(i1,j1,l,3)*SMR(i1,3,is1,js1)

            end do  ! js1
          end do  ! nb1
        end do  ! isite
        tmp2(:,:) = MR(L,:,:)
        call ZGEMM_('C','N',EXCHR,EXCHR,EXCHR,cOne,ZR,EXCHR,tmp2,EXCHR,cZero,TMP,EXCHR)
        call ZGEMM_('N','N',EXCHR,EXCHR,EXCHR,cOne,TMP,EXCHR,ZR,EXCHR,cZero,tmp2,EXCHR)
        MR(L,:,:) = tmp2(:,:)
        tmp2(:,:) = SR(L,:,:)
        call ZGEMM_('C','N',EXCHR,EXCHR,EXCHR,cOne,ZR,EXCHR,tmp2,EXCHR,cZero,TMP,EXCHR)
        call ZGEMM_('N','N',EXCHR,EXCHR,EXCHR,cOne,TMP,EXCHR,ZR,EXCHR,cZero,tmp2,EXCHR)
        SR(L,:,:) = tmp2(:,:)
        do is1=1,8
          do js1=1,8
            write(u6,'(3(A,i1),A,2F20.14)') 'MR(',l,',',is1,',',js1,') = ',MR(l,is1,js1)
          end do
        end do
      end do  ! L
      call mma_deallocate(tmp2)
      ! print the localized moments on sites:
      nsta = exchR
      ! assuming max 10 equivalent magnetic sites
      call mma_allocate(MM_tmp1,3,nsta,nsta,label='MM_tmp1')
      call mma_allocate(SM_tmp1,3,nsta,nsta,label='SM_tmp1')
      call mma_allocate(rotR_tmp,nneq,10,3,3,'rotR_tmp')
      call mma_allocate(tmp2,nsta,nsta,label='tmp2')
      call mma_allocate(tmp3,nneq,3,2,2,2,label='tmp3')
      MM_tmp1(:,:,:) = MR(:,1:nsta,1:nsta)
      SM_tmp1(:,:,:) = SR(:,1:nsta,1:nsta)
      rotR_tmp(:,:,:,:) = rotR(:,1:10,:,:)
      tmp2(:,:) = zR(1:nsta,1:nsta)
      tmp3(:,:,:,:,1) = MMR(:,:,1:nmaxR,1:nmaxR)
      tmp3(:,:,:,:,2) = SMR(:,:,1:nmaxR,1:nmaxR)
      call momloc2(nsta,nmaxR,nneq,neq,neqv,rotR_tmp,lmax,nexchR,wR(1:nsta),tmp2,MM_tmp1,SM_tmp1,tmp3(:,:,:,:,1),tmp3(:,:,:,:,2))
      call mma_deallocate(MM_tmp1)
      call mma_deallocate(SM_tmp1)
      call mma_deallocate(rotR_tmp)
      call mma_deallocate(tmp2)
      call mma_deallocate(tmp3)

      call WarningMessage(2,'Wrong code in poly_aniso/exchctl.f')
      ! FIXME: This call is missing 3 arguments
      !call barrier(exchR,MR,WR,1,2)
      call Abend()
    end if !KEOPT

  end if ! npair>0, index lp
end if !KE

!-----------------------------------------------------------------------
! ALL exchange couplings for all pairs are now known.
! printout the Hamiltonians:
#ifdef _DEBUGPRINT_
write(u6,'(A,i5)') 'exch  = ',exch
write(u6,'(A,i5)') 'npair = ',npair
write(u6,'(A,i5)') 'nneq  = ',nneq
write(u6,'(A,i5)') 'nmax  = ',nmax
write(u6,'(A,i5)') 'lmax  = ',lmax
write(u6,'(A,i5)') 'iPrint= ',iPrint
write(u6,*) 'AnisoLines1 = ',AnisoLines1
write(u6,*) 'AnisoLines3 = ',AnisoLines3
write(u6,*) 'AnisoLines9 = ',AnisoLines9
write(u6,*) 'Dipol = ',Dipol
write(u6,*) 'JITO  = ',JITO_exchange
write(u6,*) 'KE    = ',KE
do i=1,npair
  write(u6,'(A,i2,A,i3,3x,A,i2,A,i3)') 'i_pair(',i,',1) = ',i_pair(i,1),'i_pair(',i,',2) = ',i_pair(i,2)
end do
do i=1,nneq
  write(u6,'(A,i2,A,i3)') '   neq(',i,') = ',neq(i)
end do
do i=1,nneq
  write(u6,'(A,i2,A,i3)') ' nexch(',i,') = ',neq(i)
end do
do i=1,nneq
  write(u6,'(A,i2,A,100F10.3)') '   eso(',i,') = ',(eso(i,j),j=1,nexch(i))
end do
call xFlush(u6)
#endif

! print the Exchange Hamiltonian:
#ifndef _DEBUGPRINT_
if (iPrint > 2) then
#endif
  call pa_prham(exch,npair,i_pair,nneq,neq,nexch,nmax,lmax,eso,HLIN1,HLIN3,HLIN9,HDIP,HKEX,HDMO,HITO,Dipol,AnisoLines1, &
                AnisoLines3,AnisoLines9,KE,DM_exchange,JITO_exchange)
#ifndef _DEBUGPRINT_
end if
#endif

! diagonalize the Hamiltonian:
call pa_diagham(exch,npair,i_pair,nneq,neq,nexch,nmax,lmax,eso,HLIN1,HLIN3,HLIN9,HDIP,HKEX,HDMO,HITO,Dipol,DM_exchange, &
                AnisoLines1,AnisoLines3,AnisoLines9,KE,JITO_exchange,WLIN1,WLIN3,WLIN9,WLIN,WDIP,WKEX,WDMO,WITO,W,Z)

! printout the resulting eigenstates:
call pa_preigen(exch,lmax,ibas,Dipol,AnisoLines1,AnisoLines3,AnisoLines9,KE,JITO_exchange,WLIN,WDIP,WKEX,WITO,W,Z,iPrint)
!Z =  exchange eigenstates:
NmaxPop = 500
if (NmaxPop > exch) NmaxPop = exch
call PopAnalysis(nneq,neq,exch,nexch,nmax,lmax,NmaxPop,Z)
w(:) = w(:)-w(1)

! some verification
if (dnrm2_(exch,WLIN,1) > 1.0e-13_wp) call Add_Info('EXCHCTL::  WLIN',[dnrm2_(exch,WLIN,1)],1,8)
if (dnrm2_(exch,WDIP,1) > 1.0e-13_wp) call Add_Info('EXCHCTL::  WDIP',[dnrm2_(exch,WDIP,1)],1,8)
if (dnrm2_(exch,WKEX,1) > 1.0e-13_wp) call Add_Info('EXCHCTL::  WKEX',[dnrm2_(exch,WKEX,1)],1,8)
if (dnrm2_(exch,W,1) > 1.0e-13_wp) call Add_Info('EXCHCTL::     W',[dnrm2_(exch,W,1)],1,8)
! compute the moments:
M(:,:,:) = cZero
S(:,:,:) = cZero
#ifdef _DEBUGPRINT_
write(u6,'(A)') 'Magnetic moments before the build of coupled M and S matrices'
if (nPair > 0) then
  do lp=1,npair
    lb1 = i_pair(lp,1)
    lb2 = i_pair(lp,2)
    i1 = nind(lb1,1) ! indices of non-equivalent sites
    i2 = nind(lb2,1) ! indices of non-equivalent sites
    n1 = nexch(i1)
    n2 = nexch(i2)
    call mma_allocate(SM_tmp1,3,n1,n1,label='SM_tmp1')
    call mma_allocate(SM_tmp2,3,n2,n2,label='SM_tmp2')
    SM_tmp1(:,:,:) = SM(i1,:,1:n1,1:n1)
    SM_tmp2(:,:,:) = SM(i2,:,1:n2,1:n2)
    call prMom('EXCHCTL,Before M ans S, SM(i1):',SM_tmp1,n1)
    call prMom('EXCHCTL,Before M ans S, SM(i2):',SM_tmp2,n2)
    call prMom('EXCHCTL,Before M ans S, MM(i1):',MM_tmp1,n1)
    call prMom('EXCHCTL,Before M ans S, MM(i2):',MM_tmp2,n2)
    call mma_deallocate(SM_tmp1)
    call mma_deallocate(SM_tmp2)
  end do
end if
#endif

call mma_allocate(tmp2,EXCH,EXCH,label='tmp2')
do L=1,3
  do isite=1,lmax
    do nb1=1,exch
      icoord(:) = ibas(nb1,:)
      i1 = nind(isite,1)
      j1 = nind(isite,2)
      is1 = ibas(nb1,isite)+1

      do js1=1,nexch(i1)
        icoord(isite) = js1-1
        nb2 = norder(icoord,intc,lmax)
        M(l,nb1,nb2) = M(l,nb1,nb2)+rot(i1,j1,l,1)*MM(i1,1,is1,js1)+rot(i1,j1,l,2)*MM(i1,2,is1,js1)+rot(i1,j1,l,3)*MM(i1,3,is1,js1)
        S(l,nb1,nb2) = S(l,nb1,nb2)+rot(i1,j1,l,1)*SM(i1,1,is1,js1)+rot(i1,j1,l,2)*SM(i1,2,is1,js1)+rot(i1,j1,l,3)*SM(i1,3,is1,js1)

      end do  ! js1
    end do  ! nb1
  end do  ! isite

  ! magnetic moment
  tmp2(:,:) = M(L,:,:)
  call zgemm_('C','N',EXCH,EXCH,EXCH,cOne,Z,EXCH,tmp2,EXCH,cZero,TMP,EXCH)
  call zgemm_('N','N',EXCH,EXCH,EXCH,cOne,TMP,EXCH,Z,EXCH,cZero,tmp2,EXCH)
  M(L,:,:) = tmp2(:,:)
  ! spin moment
  tmp2(:,:) = S(L,:,:)
  call zgemm_('C','N',EXCH,EXCH,EXCH,cOne,Z,EXCH,tmp2,EXCH,cZero,TMP,EXCH)
  call zgemm_('N','N',EXCH,EXCH,EXCH,cOne,TMP,EXCH,Z,EXCH,cZero,tmp2,EXCH)
  S(L,:,:) = tmp2(:,:)
end do  ! L
call mma_deallocate(tmp2)

! ITO decomposition of the exchange and dipolar interactions:
if (npair > 0) &
  call pr_ito_int(npair,i_pair,lmax,nexch,nneq,neqv,itype,neq,nmax,eso,MM,SM,rot,Dipol,AnisoLines1,AnisoLines3,AnisoLines9, &
                  DM_exchange,JITO_exchange,HLIN1,HLIN3,HLIN9,HDIP,HDMO,HITO)

! projection on the Ising Hamiltonian
! accessible format
!write(u6,'(/)')
!write(u6,'(A)') 'Complete decomposition of the exchange interaction'
!write(u6,'(A)')

!-----------------------------------------------------------------------
! deallocate memory for this function:
if (lmax >= 0) then
  ! exchange energy spectrum
  call mma_deallocate(intc)
  call mma_deallocate(icoord)
  call mma_deallocate(nind)
  if (exch >= 0) call mma_deallocate(ibas)
end if
if (exch >= 0) then
  call mma_deallocate(wlin)
  call mma_deallocate(wlin1)
  call mma_deallocate(wlin3)
  call mma_deallocate(wlin9)
  call mma_deallocate(wdip)
  call mma_deallocate(wkex)
  call mma_deallocate(wdmo)
  call mma_deallocate(wito)
end if
if (nmax >= 0) then
  call mma_deallocate(S1)
  call mma_deallocate(M1)
  call mma_deallocate(S2)
  call mma_deallocate(M2)
  call mma_deallocate(SM1)
  call mma_deallocate(SM2)
  call mma_deallocate(MM1)
  call mma_deallocate(MM2)
  if (npair >= 0) then
    call mma_deallocate(HLIN1)
    call mma_deallocate(HLIN3)
    call mma_deallocate(HLIN9)
    call mma_deallocate(HDIP)
    call mma_deallocate(HKEX)
    call mma_deallocate(HDMO)
    call mma_deallocate(HITO)
  end if
end if

if (exch >= 0) call mma_deallocate(tmp)

if (nneq >= 0) then
  call mma_deallocate(nexchR)
  call mma_deallocate(SMR)
  call mma_deallocate(MMR)
  if (neqv >= 0) call mma_deallocate(rotR)
end if

if (exchR >= 0) then
  if (lmax >= 0) call mma_deallocate(ibasR)
  call mma_deallocate(WR)
  call mma_deallocate(ZR)
  call mma_deallocate(MR)
  call mma_deallocate(SR)
end if

if (npair >= 0) call mma_deallocate(HKEXR)

if (lmax >= 0) call mma_deallocate(intcR)

! results of projection of the exchange interaction on the Ising Hamiltonian:

!open(76,file='eigenstates.txt')
!do i=1,EXCH
!  do j=1,EXCH
!    write(76,'(2ES24.14)') Z(i,j)
!  end do
!end do
!close(76)
!
!open(77,file='eigenstates_full.txt')
!do i=1,EXCH
!  do j=1,EXCH
!    read(77,'(2ES24.14)') ZZR(i,j),ZZI(i,j)
!  end do
!end do
!close(77)
!ZF(:,:) = cmplx(ZZR(:,:),ZZI(:,:),kind=wp)
!call ZGEMM_('C','N',EXCH,EXCH,EXCH,cOne,Z,EXCH,ZF,EXCH,cZero,OVLP,EXCH)
!do i=1,1
!  do j=1,EXCH
!    if (ABS(OVLP(j,i)) > 1.0e-7_wp) &
!      write(u6,'(A,i3,A,i3,A,2F18.14,3x,A,F18.14,3x,A,F18.14)') 'OVLP-1(',j,',',i,') = ',OVLP(j,i),'ABS = ',ABS(OVLP(j,i)), &
!                                                               'ABS^2 = ',ABS(OVLP(j,i))*ABS(OVLP(j,i))
!  end do
!end do
!call ZGEMM_('C','N',EXCH,EXCH,EXCH,cOne,ZF,EXCH,Z,EXCH,cZero,OVLP,EXCH)
!do i=1,1
!  do j=1,EXCH
!    if (ABS(OVLP(j,i)) > 1.0e-7_wp) &
!      write(u6,'(A,i3,A,i3,A,2F18.14,3x,A,F18.14,3x,A,F18.14)') 'OVLP-2(',i,',',j,') = ',OVLP(i,j),'ABS = ',ABS(OVLP(i,j)), &
!                                                                'ABS^2 = ',ABS(OVLP(i,j))*ABS(OVLP(i,j))
!    end if
!  end do
!end do
!199 continue

return

end subroutine exchctl
