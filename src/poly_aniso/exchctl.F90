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

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
! global variables:
integer, intent(in) :: nneq
integer, intent(in) :: neqv ! max of neq(nneq)
integer, intent(in) :: neq(nneq)
integer, intent(in) :: nexch(nneq)
integer, intent(in) :: nmax
integer, intent(in) :: lmax
integer, intent(in) :: npair
integer, intent(in) :: i_pair(npair,2)
integer, intent(in) :: exch
integer, intent(in) :: lant ! (takes values from 1-7 for Gd-Yb respectively)
integer, intent(in) :: iPrint
integer, intent(in) :: mem ! memory allocated so far
integer, intent(in) :: MxRank1, MxRank2
integer, intent(in) :: imaxrank(npair,2)
character, intent(in) :: itype(nneq)
real(kind=8), intent(in) :: eso(nneq,nmax)
real(kind=8), intent(in) :: Jex(npair)
real(kind=8), intent(in) :: JAex(npair,3)
real(kind=8), intent(in) :: JDMex(npair,3)
real(kind=8), intent(in) :: JAex9(npair,3,3)
real(kind=8), intent(in) :: JITOexR(nPair,MxRank1,-MxRank1:MxRank1, MxRank2,-MxRank2:MxRank2)
real(kind=8), intent(in) :: JITOexI(nPair,MxRank1,-MxRank1:MxRank1, MxRank2,-MxRank2:MxRank2)
real(kind=8), intent(in) :: coord(nneq,3)
real(kind=8), intent(in) :: rot(nneq,neqv,3,3)
real(kind=8), intent(in) :: rlg(nneq,neqv,3,3)
real(kind=8), intent(in) :: riso(nneq,3,3)
real(kind=8), intent(in) :: tpar
real(kind=8), intent(in) :: upar
complex(kind=8), intent(inout) :: SM(nneq,3,nmax,nmax)
complex(kind=8), intent(inout) :: MM(nneq,3,nmax,nmax)
logical, intent(in) :: AnisoLines1
logical, intent(in) :: AnisoLines3
logical, intent(in) :: AnisoLines9
logical, intent(in) :: Dipol
logical, intent(in) :: KE
logical, intent(in) :: DM_exchange
logical, intent(in) :: JITO_exchange
real(kind=8), intent(out) :: W(exch)
complex(kind=8), intent(out) :: Z(exch,exch)
complex(kind=8), intent(out) :: S(3,exch,exch)
complex(kind=8), intent(out) :: M(3,exch,exch)
!------------------------------------------------------------------
! local variables
integer :: i, j, l, lp, lb1, lb2, nb, nb1, nb2, isite, is1, js1, i1, i2, j1, j2, k, ibuf, q1, q2, k1, k2, n1, n2
integer :: norder
integer :: NmaxPop
integer, allocatable :: intc(:)   !  intc(lmax)
integer, allocatable :: ibas(:,:) !  ibas(exch,lmax)
integer, allocatable :: icoord(:) !  icoord(lmax)
integer, allocatable :: nind(:,:) !  nind(lmax,2)
real(kind=8) :: vect(3)
real(kind=8) :: dist
real(kind=8), allocatable :: wlin(:) ! wlin(exch)
real(kind=8), allocatable :: wlin1(:)! wlin1(exch)
real(kind=8), allocatable :: wlin3(:)! wlin3(exch)
real(kind=8), allocatable :: wlin9(:)! wlin9(exch)
real(kind=8), allocatable :: wdip(:) ! wdip(exch)
real(kind=8), allocatable :: wkex(:) ! wkex(exch)
real(kind=8), allocatable :: wdmo(:) ! wdmo(exch)
real(kind=8), allocatable :: wito(:) ! wito(exch)
complex(kind=8), allocatable :: S1(:,:,:) ! S1(3,nmax,nmax)
complex(kind=8), allocatable :: M1(:,:,:) ! M1(3,nmax,nmax)
complex(kind=8), allocatable :: S2(:,:,:) ! S2(3,nmax,nmax)
complex(kind=8), allocatable :: M2(:,:,:) ! M2(3,nmax,nmax)
complex(kind=8), allocatable :: ZA1(:,:), ZA2(:,:)
complex(kind=8), allocatable :: SM1(:,:,:) ! SM1(3,nmax,nmax)
complex(kind=8), allocatable :: MM1(:,:,:) ! MM1(3,nmax,nmax)
complex(kind=8), allocatable :: SM2(:,:,:) ! SM2(3,nmax,nmax)
complex(kind=8), allocatable :: MM2(:,:,:) ! MM2(3,nmax,nmax)
complex(kind=8), allocatable :: HLIN1(:,:,:,:,:) ! HLIN1(npair,nmax,nmax,nmax,nmax)
complex(kind=8), allocatable :: HLIN3(:,:,:,:,:) ! HLIN3(npair,nmax,nmax,nmax,nmax)
complex(kind=8), allocatable :: HLIN9(:,:,:,:,:) ! HLIN9(npair,nmax,nmax,nmax,nmax)
complex(kind=8), allocatable :: HDIP(:,:,:,:,:) ! HDIP(npair,nmax,nmax,nmax,nmax)
complex(kind=8), allocatable :: HKEX(:,:,:,:,:) ! HKEX(npair,nmax,nmax,nmax,nmax)
complex(kind=8), allocatable :: HDMO(:,:,:,:,:) ! HDMO(npair,nmax,nmax,nmax,nmax)
complex(kind=8), allocatable :: HITO(:,:,:,:,:) ! HITO(npair,nmax,nmax,nmax,nmax)
complex(kind=8), allocatable :: tmp(:,:) ! tmp(exch,exch)
! two options for KE:
integer :: KEOPT
integer, parameter :: exchR = 8
integer :: nmaxR
integer :: nsta
integer, allocatable :: nexchR(:)  ! nexchR(nneq)
integer, allocatable :: ibasR(:,:) ! ibasR(exchR,lmax)
integer, allocatable :: intcR(:)   ! intcR(lmax)
real(kind=8), allocatable :: WR(:)  ! WR(exchR)
real(kind=8), allocatable :: rotR(:,:,:,:) ! rotR(nneq,neqv,3,3)
complex(kind=8), allocatable :: ZR(:,:) ! ZR(exchR,exchR)
complex(kind=8), allocatable :: HKEXR(:,:,:,:,:) ! HKEXR(npair,2,2,2,2)
complex(kind=8), allocatable :: MR(:,:,:) ! MR(3,exchR,exchR)
complex(kind=8), allocatable :: SR(:,:,:) ! SR(3,exchR,exchR)
complex(kind=8), allocatable :: SMR(:,:,:,:) ! SMR(nneq,3,2,2)
complex(kind=8), allocatable :: MMR(:,:,:,:) ! MMR(nneq,3,2,2)
!complex(kind=8) :: JAllDip(npair,nmax,-nmax:nmax,nmax,-nmax:nmax)
!complex(kind=8) :: JAllEx(npair,nmax,-nmax:nmax,nmax,-nmax:nmax)
!real(kind=8) :: J1Dip(npair,3,3)
!real(kind=8) :: J1Ex(npair,3,3)
real(kind=8) :: mg1(3,3), mg2(3,3)
integer :: CtoB, RtoB, ItoB, mem_local
logical :: DBG !, testlines
real(kind=8) :: dnrm2_
external :: norder, dnrm2_  !,ilaenv

DBG = .false.
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
if (DBG) then
  write(6,'(A)') 'Enter EXCHCTL'
  write(6,'(A,  i8)') 'exch   = ',exch
  write(6,'(A,  i8)') 'nneq   = ',nneq
  write(6,'(A,  i8)') 'neqv   = ',neqv
  write(6,'(A,  i8)') 'nmax   = ',nmax
  write(6,'(A,  i8)') 'lmax   = ',lmax
  write(6,'(A,  i8)') 'nPair  = ',nPair
  write(6,'(A,  i8)') 'nPair  = ',nPair
  write(6,'(A,  i8)') 'MxRank1= ',MxRank1
  write(6,'(A,  i8)') 'MxRank2= ',MxRank2
  write(6,'(A,10i4)') 'neq()  = ',(neq(i),i=1,nneq)
  write(6,'(A,10i4)') 'nexch()= ',(nexch(i),i=1,nneq)
  write(6,'(A,  L2)') 'AnisoLines1  = ',AnisoLines1
  write(6,'(A,  L2)') 'AnisoLines3  = ',AnisoLines3
  write(6,'(A,  L2)') 'AnisoLines9  = ',AnisoLines9
  write(6,'(A,  L2)') 'DM_exchange  = ',DM_exchange
  write(6,'(A,  L2)') 'Dipol        = ',Dipol
  write(6,'(A,  L2)') 'JITO_exchange= ',JITO_exchange
  if (AnisoLines1) then
    do i=1,nPair
      write(6,'(A,2I4,F10.5)') 'LIN1',i_pair(i,1),i_pair(i,2),Jex(i)
    end do
  end if

  if (AnisoLines3) then
    do i=1,nPair
      write(6,'(A,2I4,3F10.5)') 'LIN3',i_pair(i,1),i_pair(i,2),(JAex(i,j),j=1,3)
    end do
  end if

  if (AnisoLines9) then
    do i=1,nPair
      write(6,'(A,2I4,9F10.5)') 'LIN9',i_pair(i,1),i_pair(i,2),((JAex9(i,j,k),j=1,3),k=1,3)
    end do
  end if

  if (DM_exchange) then
    do i=1,nPair
      write(6,'(A,2I4,3F10.5)') 'DMEX',i_pair(i,1),i_pair(i,2),(JDMex(i,j),j=1,3)
    end do
  end if

  if (Dipol) then
    write(6,'(A)') 'COORD(i):'
    do i=1,nneq
      write(6,'(i3,3F14.8)') i,(coord(i,j),j=1,3)
    end do
  end if

  if (JITO_exchange) then
    do i=1,nPair
      write(6,'(A,4I4,3F10.5)') 'ITOJ',i_pair(i,1),i_pair(i,2),imaxrank(i,1),imaxrank(i,2)
      do k1=1,imaxrank(i,1),2
        do q1=-k1,k1
          do k2=1,imaxrank(i,2),2
            do q2=-k2,k2
              write(6,'(4I3,2x,2ES21.14)') k1,q1,k2,q2,JITOexR(i,k1,q1,k2,q2),JITOexI(i,k1,q1,k2,q2)
            end do
          end do
        end do
      end do
    end do ! ipair
  end if

  write(6,'(A)') 'ESO(i):'
  do i=1,nneq
    write(6,'(i3,90F14.8)') i,(eso(i,j),j=1,nmax)
  end do
  write(6,'(90A)') 'itype()=',(itype(i),' ',i=1,nneq)

  write(6,'(A,  i8)') 'neqv   = ',neqv
  do i=1,nneq
    write(6,'(A,  i8)') 'riso( site=',i,'):'
    do j=1,3
      write(6,'(3ES22.14)') (riso(i,j,k),k=1,3)
    end do
  end do
end if ! dbg
!-----------------------------------------------------------------------
! allocate memory for this function:
ItoB = 8
RtoB = 8
CtoB = 16
mem_local = 0
if (lmax >= 0) then
  ! exchange energy spectrum
  call mma_allocate(intc,lmax,'intc')
  call mma_allocate(icoord,lmax,'icoord')
  call mma_allocate(nind,lmax,2,'nind')
  call icopy(lmax,[0],0,intc,1)
  call icopy(lmax,[0],0,icoord,1)
  call icopy(2*lmax,[0],0,nind,1)
  mem_local = mem_local+4*lmax*ItoB
  if (exch >= 0) then
    call mma_allocate(ibas,exch,lmax,'ibas')
    call icopy(exch*lmax,[0],0,ibas,1)
    mem_local = mem_local+exch*lmax*ItoB
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
  call dcopy_(exch,[0.0_wp],0,wlin,1)
  call dcopy_(exch,[0.0_wp],0,wlin1,1)
  call dcopy_(exch,[0.0_wp],0,wlin3,1)
  call dcopy_(exch,[0.0_wp],0,wlin9,1)
  call dcopy_(exch,[0.0_wp],0,wdip,1)
  call dcopy_(exch,[0.0_wp],0,wkex,1)
  call dcopy_(exch,[0.0_wp],0,wdmo,1)
  call dcopy_(exch,[0.0_wp],0,wito,1)
  mem_local = mem_local+8*exch*RtoB
end if
if (nmax >= 0) then
  call mma_allocate(S1,3,nmax,nmax,' S1')
  call mma_allocate(M1,3,nmax,nmax,' M1')
  call mma_allocate(S2,3,nmax,nmax,' S2')
  call mma_allocate(M2,3,nmax,nmax,' M2')
  call mma_allocate(ZA1,nmax,nmax,' Z1')
  call mma_allocate(ZA2,nmax,nmax,' Z2')
  call mma_allocate(SM1,3,nmax,nmax,'SM1')
  call mma_allocate(SM2,3,nmax,nmax,'SM2')
  call mma_allocate(MM1,3,nmax,nmax,'MM1')
  call mma_allocate(MM2,3,nmax,nmax,'MM2')
  call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,S1,1)
  call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,M1,1)
  call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,S2,1)
  call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,M2,1)
  call zcopy_(nmax*nmax,[(0.0_wp,0.0_wp)],0,ZA1,1)
  call zcopy_(nmax*nmax,[(0.0_wp,0.0_wp)],0,ZA2,1)
  call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,SM1,1)
  call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,SM2,1)
  call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,MM1,1)
  call zcopy_(3*nmax*nmax,[(0.0_wp,0.0_wp)],0,MM2,1)
  mem_local = mem_local+8*3*nmax*nmax*CtoB

  if (npair >= 0) then
    ibuf = npair*nmax*nmax*nmax*nmax
    call mma_allocate(HLIN1,npair,nmax,nmax,nmax,nmax,'HLIN1')
    call mma_allocate(HLIN3,npair,nmax,nmax,nmax,nmax,'HLIN3')
    call mma_allocate(HLIN9,npair,nmax,nmax,nmax,nmax,'HLIN9')
    call mma_allocate(HDIP,npair,nmax,nmax,nmax,nmax,'HDIP')
    call mma_allocate(HKEX,npair,nmax,nmax,nmax,nmax,'HKEX')
    call mma_allocate(HDMO,npair,nmax,nmax,nmax,nmax,'HDMO')
    call mma_allocate(HITO,npair,nmax,nmax,nmax,nmax,'HITO')
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN1,1)
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN3,1)
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN9,1)
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HDIP,1)
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HKEX,1)
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HDMO,1)
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HITO,1)
    mem_local = mem_local+7*ibuf*CtoB
  end if
end if

if (exch >= 0) then
  call mma_allocate(tmp,exch,exch,'tmp')
  call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,tmp,1)
  mem_local = mem_local+exch*exch*CtoB
end if

if (nneq >= 0) then
  call mma_allocate(nexchR,nneq,'nexchR')
  call icopy(nneq,[0],0,nexchR,1)
  mem_local = mem_local+nneq*ItoB

  call mma_allocate(SMR,nneq,3,2,2,'SMR')
  call mma_allocate(MMR,nneq,3,2,2,'MMR')
  call zcopy_(nneq*3*2*2,[(0.0_wp,0.0_wp)],0,SMR,1)
  call zcopy_(nneq*3*2*2,[(0.0_wp,0.0_wp)],0,MMR,1)
  mem_local = mem_local+2*nneq*3*2*2*CtoB

  if (neqv >= 0) then
    call mma_allocate(rotR,nneq,neqv,3,3,'rotR')
    call dcopy_(nneq*neqv*3*3,[0.0_wp],0,rotR,1)
    mem_local = mem_local+nneq*neqv*3*3*RtoB
  end if
end if

if (exchR >= 0) then
  if (lmax >= 0) then
    call mma_allocate(ibasR,nneq,lmax,'ibasR')
    call icopy(nneq*lmax,[0],0,ibasR,1)
    mem_local = mem_local+nneq*lmax*ItoB
  end if
  call mma_allocate(WR,exchR,'WR')
  call dcopy_(exchR,[0.0_wp],0,WR,1)
  mem_local = mem_local+exchR*RtoB

  call mma_allocate(ZR,exchR,exchR,'ZR')
  call mma_allocate(MR,3,exchR,exchR,'MR')
  call mma_allocate(SR,3,exchR,exchR,'SR')
  call zcopy_(exchR*exchR,[(0.0_wp,0.0_wp)],0,ZR,1)
  call zcopy_(3*exchR*exchR,[(0.0_wp,0.0_wp)],0,MR,1)
  call zcopy_(3*exchR*exchR,[(0.0_wp,0.0_wp)],0,SR,1)
  mem_local = mem_local+7*exchR*exchR*CtoB
end if

if (npair >= 0) then
  call mma_allocate(HKEXR,npair,2,2,2,2,'HKEXR')
  call zcopy_(npair*2*2*2*2,[(0.0_wp,0.0_wp)],0,HKEXR,1)
  mem_local = mem_local+npair*2*2*2*2*CtoB
end if

if (lmax >= 0) then
  call mma_allocate(intcR,lmax,'intcR')
  call icopy(lmax,[0],0,intcR,1)
  mem_local = mem_local+lmax*ItoB
end if
if (dbg) write(6,*) 'EXCHCTL:  memory allocated (local):'
if (dbg) write(6,*) 'mem_local=',mem_local
if (dbg) write(6,*) 'EXCHCTL:  memory allocated (total):'
if (dbg) write(6,*) 'mem_total=',mem+mem_local
!-----------------------------------------------------------------------
l = 0
do i=1,nneq
  do j=1,neq(i)
    l = l+1
    nind(l,1) = i
    nind(l,2) = j
  end do
end do
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
if (dbg) then
  write(6,'(34x,A,1x,20i3)') 'site Nr.',(i,i=1,lmax)
  do nb=1,exch
    write(6,'(A,i5,A,20i3)') 'COUPLING: basis set:  ibas(',nb,' ,isite) = ',(ibas(nb,i)+1,i=1,lmax)
  end do
end if ! dbg

!-----------------------------------------------------------------------
! Lines model of magnetic couping  -- 1 parameter
if (AnisoLines1) then
  if (nPair > 0) then
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN1,1)
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        call prep_mom_exchange(n1,rot(i1,j1,1:3,1:3),SM(i1,1:3,1:n1,1:n1),MM(i1,1:3,1:n1,1:n1),mg1,.true.)
      end if
      if (itype(i2) == 'A') then
        call prep_mom_exchange(n2,rot(i2,j2,1:3,1:3),SM(i2,1:3,1:n2,1:n2),MM(i2,1:3,1:n2,1:n2),mg2,.true.)
      end if

      !if (dbg)
      call prMom('SM(i1) bf Lines1',SM(i1,1:3,1:n1,1:n1),n1)
      !if (dbg)
      call prMom('SM(i2) bf Lines1',SM(i2,1:3,1:n2,1:n2),n2)

      ! build the Lines exchange matrix:
      call Lines_Exchange(Jex(lp),n1,n2,SM(i1,1:3,1:n1,1:n1),SM(i2,1:3,1:n2,1:n2),HLIN1(lp,1:n1,1:n1,1:n2,1:n2))

      if (dbg) then
        call prMom('SM(i1) af Lines1',SM(i1,1:3,1:n1,1:n1),n1)
        call prMom('SM(i2) af Lines1',SM(i2,1:3,1:n2,1:n2),n2)
        write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
        write(6,'(A,i2)') 'in local pseudospin basis'
        do i=1,n1
          do j=1,n1
            do k=1,n2
              do l=1,n2
                write(6,'(4(a,i2),a,2ES24.14)') 'HLIN1(',i,',',j,',',k,',',l,') = ',HLIN1(lp,i,j,k,l)
              end do
            end do
          end do
        end do
      end if! dbg

    end do
  end if ! nPair
end if
!-----------------------------------------------------------------------

! Anisotropic Lines model of magnetic couping -- 3 parameters)
! Jxx, Jyy, Jzz
if (AnisoLines3) then
  if (nPair > 0) then
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN3,1)
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        call prep_mom_exchange(n1,rot(i1,j1,1:3,1:3),SM(i1,1:3,1:n1,1:n1),MM(i1,1:3,1:n1,1:n1),mg1,dbg)
      end if

      if (itype(i2) == 'A') then
        call prep_mom_exchange(n2,rot(i2,j2,1:3,1:3),SM(i2,1:3,1:n2,1:n2),MM(i2,1:3,1:n2,1:n2),mg2,dbg)
      end if

      call Aniso_Lines_Exchange3(JAex(lp,1:3),n1,n2,SM(i1,1:3,1:n1,1:n1),SM(i2,1:3,1:n2,1:n2),HLIN3(lp,1:n1,1:n1,1:n2,1:n2))

      if (dbg) then
        write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
        do i=1,n1
          do j=1,n1
            do k=1,n2
              do l=1,n2
                write(6,'(4(a,i2),a,2ES24.14)') 'HLIN3(',i,',',j,',',k,',',l,') = ',HLIN3(lp,i,j,k,l)
              end do
            end do
          end do
        end do
      end if ! dbg

    end do ! lp
  end if ! nPair
end if
!-----------------------------------------------------------------------

! Anisotropic Lines model of magnetic couping -- 9 parameters)
! Jxx, Jxy, Jxz, Jyx, Jyy, Jyz, Jzx, Jzy, Jzz
if (AnisoLines9) then
  if (nPair > 0) then
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HLIN9,1)
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        call prep_mom_exchange(n1,rot(i1,j1,1:3,1:3),SM(i1,1:3,1:n1,1:n1),MM(i1,1:3,1:n1,1:n1),mg1,dbg)
      end if
      if (itype(i2) == 'A') then
        call prep_mom_exchange(n2,rot(i2,j2,1:3,1:3),SM(i2,1:3,1:n2,1:n2),MM(i2,1:3,1:n2,1:n2),mg2,dbg)
      end if
      ! rotate the input J matrix by the mg1 and mg2
      if (dbg) call prMom('SM(i1) bf Lines9',SM(i1,1:3,1:n1,1:n1),n1)
      if (dbg) call prMom('SM(i2) bf Lines9',SM(i2,1:3,1:n2,1:n2),n2)

      call Aniso_Lines_Exchange9(JAex9(lp,1:3,1:3),n1,n2,SM(i1,1:3,1:n1,1:n1),SM(i2,1:3,1:n2,1:n2),HLIN9(lp,1:n1,1:n1,1:n2,1:n2))

      if (dbg) then
        call prMom('SM(i1) af Lines9',SM(i1,1:3,1:n1,1:n1),n1)
        call prMom('SM(i2) af Lines9',SM(i2,1:3,1:n2,1:n2),n2)
        write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
        do i=1,n1
          do j=1,n1
            do k=1,n2
              do l=1,n2
                write(6,'(4(a,i2),a,2ES24.14)') 'HLIN9(',i,',',j,',',k,',',l,') = ',HLIN9(lp,i,j,k,l)
              end do
            end do
          end do
        end do
      end if! dbg

    end do
  end if ! nPair
end if

!-----------------------------------------------------------------------
! dipolar couping
if (Dipol) then
  if (nPair > 0) then
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HDIP,1)
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)

      vect(1:3) = 0.0_wp
      dist = 0.0_wp
      call dirvect(coord(i1,1:3),rlg(i1,j1,1:3,1:3),coord(i2,1:3),rlg(i2,j2,1:3,1:3),vect(1:3),dist)
      if (DBG) write(6,'(A,i3,3ES20.12,2x,ES20.12)') 'EXCHCTL: DIPOL: lp, vect, R:',lp,(vect(i),i=1,3),dist
      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        call prep_mom_exchange(n1,rot(i1,j1,1:3,1:3),SM(i1,1:3,1:n1,1:n1),MM(i1,1:3,1:n1,1:n1),mg1,dbg)
      end if
      if (itype(i2) == 'A') then
        call prep_mom_exchange(n2,rot(i2,j2,1:3,1:3),SM(i2,1:3,1:n2,1:n2),MM(i2,1:3,1:n2,1:n2),mg2,dbg)
      end if

      call Dipol_Exchange(n1,n2,vect(1:3),dist,MM(i1,1:3,1:n1,1:n1),MM(i2,1:3,1:n2,1:n2),HDIP(lp,1:n1,1:n1,1:n2,1:n2))

      if (dbg) then
        write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
        do i=1,n1
          do j=1,n1
            do k=1,n2
              do l=1,n2
                write(6,'(4(a,i2),a,2ES24.14)') 'HDIP (',i,',',j,',',k,',',l,') = ',HDIP(lp,i,j,k,l)
              end do
            end do
          end do
        end do
      end if ! dbg

    end do ! lp
  end if
end if
!-----------------------------------------------------------------------

! Dzyaloshinsky-Morya antisymmetric couping
if (DM_exchange) then
  if (nPair > 0) then
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HDMO,1)
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)
      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        call prep_mom_exchange(n1,rot(i1,j1,1:3,1:3),SM(i1,1:3,1:n1,1:n1),MM(i1,1:3,1:n1,1:n1),mg1,dbg)
      end if
      if (itype(i2) == 'A') then
        call prep_mom_exchange(n2,rot(i2,j2,1:3,1:3),SM(i2,1:3,1:n2,1:n2),MM(i2,1:3,1:n2,1:n2),mg2,dbg)
      end if

      call Dzyaloshinsky_Morya_Exchange(JDMex(lp,1:3),n1,n2,SM(i1,1:3,1:n1,1:n1),SM(i2,1:3,1:n2,1:n2),HDMO(lp,1:n1,1:n1,1:n2,1:n2))
      if (dbg) then
        write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
        do i=1,n1
          do j=1,n1
            do k=1,n2
              do l=1,n2
                write(6,'(4(a,i2),a,2ES24.14)') 'HDMO (',i,',',j,',',k,',',l,') = ',HDMO(lp,i,j,k,l)
              end do
            end do
          end do
        end do
      end if! dbg

    end do
  end if ! nPair
end if
!-----------------------------------------------------------------------

! JITO exchange interaction
if (JITO_exchange) then
  if (dbg) write(6,'(A)') 'EXCHCTL:  Entering  JITO_exchange'
  if (nPair > 0) then
    call zcopy_(ibuf,[(0.0_wp,0.0_wp)],0,HITO,1)
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)
      ! find local pseudospin and rotate the spin and magnetic moment
      ! to the local pseudospin basis
      if (itype(i1) == 'A') then
        call prep_mom_exchange(n1,rot(i1,j1,1:3,1:3),SM(i1,1:3,1:n1,1:n1),MM(i1,1:3,1:n1,1:n1),mg1,dbg)
      end if
      if (itype(i2) == 'A') then
        call prep_mom_exchange(n2,rot(i2,j2,1:3,1:3),SM(i2,1:3,1:n2,1:n2),MM(i2,1:3,1:n2,1:n2),mg2,dbg)
      end if

      ! using Naoya's ITO:, in general coordinate system
      call JITO_Exchange_Int(MxRank1,MxRank2,imaxrank(lp,1:2),n1,n2, &
                             JITOexR(lp,1:MxRank1,-MxRank1:MxRank1,1:MxRank2,-MxRank2:MxRank2), &
                             JITOexI(lp,1:MxRank1,-MxRank1:MxRank1,1:MxRank2,-MxRank2:MxRank2),HITO(lp,1:n1,1:n1,1:n2,1:n2))
      if (dbg) then
        write(6,'(A,i2)') 'Exchange matrix for pair = ',lp
        do i=1,n1
          do j=1,n1
            do k=1,n2
              do l=1,n2
                write(6,'(4(a,i2),a,2ES24.14)') 'HITO (',i,',',j,',',k,',',l,') = ',HITO(lp,i,j,k,l)
              end do
            end do
          end do
        end do
      end if! dbg

    end do
  end if ! nPair
  if (dbg) write(6,'(A)') 'EXCHCTL:  Exiting JITO_exchange'
end if
!-----------------------------------------------------------------------

if (KE) then
  if (nPair > 0) then
    HKEX = (0.0_wp,0.0_wp)
    HKEX = (0.0_wp,0.0_wp)
    HKEXR = (0.0_wp,0.0_wp)
    MMR = (0.0_wp,0.0_wp)
    SMR = (0.0_wp,0.0_wp)
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      j1 = nind(lb1,2) ! indices of equivalent sites
      j2 = nind(lb2,2) ! indices of equivalent sites

      n1 = nexch(i1)
      n2 = nexch(i2)
      call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,S1,1)
      call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,S2,1)
      call zcopy_(3*n1*n1,[(0.0_wp,0.0_wp)],0,M1,1)
      call zcopy_(3*n2*n2,[(0.0_wp,0.0_wp)],0,M2,1)
      call rotmom2(MM(i1,1:3,1:n1,1:n1),n1,rot(i1,j1,1:3,1:3),M1(1:3,1:n1,1:n1))
      call rotmom2(SM(i1,1:3,1:n1,1:n1),n1,rot(i1,j1,1:3,1:3),S1(1:3,1:n1,1:n1))
      call rotmom2(MM(i2,1:3,1:n2,1:n2),n2,rot(i2,j2,1:3,1:3),M2(1:3,1:n2,1:n2))
      call rotmom2(SM(i2,1:3,1:n2,1:n2),n2,rot(i2,j2,1:3,1:3),S2(1:3,1:n2,1:n2))
      ! KEOPT=1 ! FULL
      ! KEOPT=2 ! Full, 1/U
      ! KEOPT=3 ! FULL  + reduced form
      ! KEOPT=4 ! Full, 1/U + reduced form
      MM1 = (0.0_wp,0.0_wp)
      SM1 = (0.0_wp,0.0_wp)
      MM2 = (0.0_wp,0.0_wp)
      SM2 = (0.0_wp,0.0_wp)
      !FIXME: this call to KE_Exchange does not match at all its definition, please fix
      call WarningMessage(2,'There is surely a bug here')
      if (.false.) call Unused_real(tpar)
      if (.false.) call Unused_real(upar)
      if (.false.) call Unused_integer(lant)
      !call KE_Exchange(n1,n2,M1(1:3,1:n1,1:n1),S1(1:3,1:n1,1:n1),M2(1:3,1:n2,1:n2),S2(1:3,1:n2,1:n2),eso(i1,1:n1),eso(i2,1:n2), &
      !                 tpar,upar,lant,KEOPT,HKEX(lp,1:n1,1:n1,1:n2,1:n2),MM1(1:3,1:n1,1:n1),SM1(1:3,1:n1,1:n1),MM2(1:3,1:n2,1:n2), &
      !                 SM2(1:3,1:n2,1:n2)

      if ((KEOPT == 1) .or. (KEOPT == 2)) then
        do is1=1,n2
          do js1=1,n2
            HKEXR(lp,1,1,is1,js1) = HKEX(lp,1,1,is1,js1)
            HKEXR(lp,1,2,is1,js1) = HKEX(lp,1,2,is1,js1)
            HKEXR(lp,2,1,is1,js1) = HKEX(lp,2,1,is1,js1)
            HKEXR(lp,2,2,is1,js1) = HKEX(lp,2,2,is1,js1)
          end do
        end do
        do l=1,3
          do is1=1,2
            do js1=1,2
              MMR(i1,l,is1,js1) = MM1(l,is1,js1)
              MMR(i2,l,is1,js1) = MM2(l,is1,js1)
              SMR(i1,l,is1,js1) = SM1(l,is1,js1)
              SMR(i2,l,is1,js1) = SM2(l,is1,js1)
            end do
          end do
        end do

      else if ((KEOPT == 3) .or. (KEOPT == 4)) then
        do is1=1,n2
          do js1=1,n2
            HKEXR(lp,1,1,is1,js1) = HKEX(lp,1,1,is1,js1)
            HKEXR(lp,1,2,is1,js1) = HKEX(lp,1,n1,is1,js1)
            HKEXR(lp,2,1,is1,js1) = HKEX(lp,n1,1,is1,js1)
            HKEXR(lp,2,2,is1,js1) = HKEX(lp,n1,n1,is1,js1)
          end do
        end do

        do l=1,3
          MMR(i1,l,1,1) = MM1(l,1,1)
          MMR(i1,l,1,2) = MM1(l,1,n1)
          MMR(i1,l,2,1) = MM1(l,n1,1)
          MMR(i1,l,2,2) = MM1(l,n1,n1)

          MMR(i2,l,1,1) = MM2(l,1,1)
          MMR(i2,l,1,2) = MM2(l,1,n2)
          MMR(i2,l,2,1) = MM2(l,n2,1)
          MMR(i2,l,2,2) = MM2(l,n2,n2)

          SMR(i1,l,1,1) = SM1(l,1,1)
          SMR(i1,l,1,2) = SM1(l,1,n1)
          SMR(i1,l,2,1) = SM1(l,n1,1)
          SMR(i1,l,2,2) = SM1(l,n1,n1)

          SMR(i2,l,1,1) = SM2(l,1,1)
          SMR(i2,l,1,2) = SM2(l,1,n2)
          SMR(i2,l,2,1) = SM2(l,n2,1)
          SMR(i2,l,2,2) = SM2(l,n2,n2)
        end do

        if (DBG) then
          write(6,'(A)') 'site Ln'
          do i=1,2
            do j=1,2
              write(6,'(3(A,i1),A,3(2F20.14,2x))') 'MMR(',i1,',L,',i,',',j,')= ',(MMR(i1,l,i,j),l=1,3)
            end do
          end do
          write(6,'(/)')
          do i=1,2
            do j=1,2
              write(6,'(3(A,i1),A,3(2F20.14,2x))') 'SMR(',i1,',L,',i,',',j,')= ',(SMR(i1,l,i,j),l=1,3)
            end do
          end do
          write(6,'(/)')
          write(6,'(A)') 'site Radical'
          do i=1,2
            do j=1,2
              write(6,'(3(A,i1),A,3(2F20.14,2x))') 'MMR(',i2,',L,',i,',',j,')= ',(MMR(i2,l,i,j),l=1,3)
            end do
          end do
          write(6,'(/)')
          do i=1,2
            do j=1,2
              write(6,'(3(A,i1),A,3(2F20.14,2x))') 'SMR(',i2,',L,',i,',',j,')= ',(SMR(i2,l,i,j),l=1,3)
            end do
          end do
        end if ! DBG

      end if ! KEOPT
    end do ! lp

    ! in case of KEOPT=3 or KEOPT=4 Then we need to compute the spectrum and the properties
    ! in the reduced form, where nexch(i1)=2 ( ground doublet only).
    ! exchnew=8:
    if (KEOPT <= 4) then
      nmaxR = 2
      nexchR(1) = 2
      nexchR(2) = 2
      do i=1,nneq
        nexchR(i) = 2
      end do
      intcR(:) = 0
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
      HLIN1 = (0.0_wp,0.0_wp)
      HLIN3 = (0.0_wp,0.0_wp)
      HLIN9 = (0.0_wp,0.0_wp)
      HDIP = (0.0_wp,0.0_wp)
      WLIN = 0.0_wp
      WLIN1 = 0.0_wp
      WLIN3 = 0.0_wp
      WLIN9 = 0.0_wp
      WDIP = 0.0_wp
      WKEX = 0.0_wp
      WR = 0.0_wp
      ZR = (0.0_wp,0.0_wp)
      ! print the Exchange Hamiltonian:
      call pa_prham(exchR,npair,i_pair,nneq,neq,nexchR,nmaxR,lmax,eso(1:nneq,1:nmaxR), &
                    HLIN1(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),HLIN3(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR), &
                    HLIN9(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),HDIP(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR), &
                    HKEXR(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),HDMO(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR), &
                    HITO(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),Dipol,AnisoLines1,AnisoLines3,AnisoLines9,KE,DM_exchange,.false.)
      ! diagonalize the Hamiltonian:
      call pa_diagham(exchR,npair,i_pair,nneq,neq,nexchR,nmaxR,lmax,eso(1:nneq,1:nmaxR), &
                      HLIN1(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),HLIN3(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR), &
                      HLIN9(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),HDIP(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR), &
                      HKEXR(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),HDMO(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR), &
                      HITO(1:npair,1:nmaxR,1:nmaxR,1:nmaxR,1:nmaxR),Dipol,.false.,AnisoLines1,AnisoLines3,AnisoLines9,KE,.false., &
                      WLIN1(1:exchR),WLIN3(1:exchR),WLIN9(1:exchR),WLIN(1:exchR),WDIP(1:exchR),WKEX(1:exchR),WDMO(1:exchR), &
                      WITO(1:exchR),WR(1:exchR),ZR(1:exchR,1:exchR))
      ! print the resulting eigenstates:
      call pa_preigen(exchR,lmax,ibasR,Dipol,AnisoLines1,AnisoLines3,AnisoLines9,KE,.false.,.false.,WLIN(1:exchR),WDIP(1:exchR), &
                      WKEX(1:exchR),WDMO(1:exchR),WITO(1:exchR),WR(1:exchR),ZR(1:exchR,1:exchR),0)
      ! compute the moments:
      rotR = 0.0_wp
      rotR(1,1,1,1) = 1.0_wp
      rotR(1,1,2,2) = 1.0_wp
      rotR(1,1,3,3) = 1.0_wp
      rotR(1,2,1,1) = 1.0_wp
      rotR(1,2,2,2) = 1.0_wp
      rotR(1,2,3,3) = 1.0_wp
      rotR(2,1,1,1) = 1.0_wp
      rotR(2,1,2,2) = 1.0_wp
      rotR(2,1,3,3) = 1.0_wp
      MR = (0.0_wp,0.0_wp)
      SR = (0.0_wp,0.0_wp)
      do L=1,3
        do isite=1,lmax
          do nb1=1,exchR
            do lp=1,lmax
              icoord(lp) = ibasR(nb1,lp)
            end do
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
        TMP(:,:) = (0.0_wp,0.0_wp)
        call ZGEMM_('C','N',EXCHR,EXCHR,EXCHR,(1.0_wp,0.0_wp),ZR,EXCHR,MR(L,:,:),EXCHR,(0.0_wp,0.0_wp),TMP,EXCHR)
        MR(L,:,:) = (0.0_wp,0.0_wp)
        call ZGEMM_('N','N',EXCHR,EXCHR,EXCHR,(1.0_wp,0.0_wp),TMP,EXCHR,ZR,EXCHR,(0.0_wp,0.0_wp),MR(L,:,:),EXCHR)
        TMP(:,:) = (0.0_wp,0.0_wp)
        call ZGEMM_('C','N',EXCHR,EXCHR,EXCHR,(1.0_wp,0.0_wp),ZR,EXCHR,SR(L,:,:),EXCHR,(0.0_wp,0.0_wp),TMP,EXCHR)
        SR(L,:,:) = (0.0_wp,0.0_wp)
        call ZGEMM_('N','N',EXCHR,EXCHR,EXCHR,(1.0_wp,0.0_wp),TMP,EXCHR,ZR,EXCHR,(0.0_wp,0.0_wp),SR(L,:,:),EXCHR)
        do is1=1,8
          do js1=1,8
            write(6,'(3(A,i1),A,2F20.14)') 'MR(',l,',',is1,',',js1,') = ',MR(l,is1,js1)
          end do
        end do
      end do  ! L
      ! print the localized moments on sites:
      nsta = exchR
      ! assuming max 10 equivalent magnetic sites
      call momloc2(nsta,nmaxR,nneq,neq,neqv,rotR(1:nneq,1:10,:,:),lmax,nexchR,wR(1:nsta),zR(1:nsta,1:nsta),MR(1:3,1:nsta,1:nsta), &
                   SR(1:3,1:nsta,1:nsta),MMR(1:nneq,1:3,1:nmaxR,1:nmaxR),SMR(1:nneq,1:3,1:nmaxR,1:nmaxR))

      call WarningMessage(2,'Wrong code in poly_aniso/exchctl.f')
      ! FIXME: This call is missing 3 arguments
      !call barrier(exchR,MR(1:3,1:exchR,1:exchR),WR(1:exchR),1,2)
      call Abend()
    end if !KEOPT

  end if ! npair>0, index lp
end if !KE

!-----------------------------------------------------------------------
! ALL exchange couplings for all pairs are now known.
! printout the Hamiltonians:
if (DBG) then
  write(6,'(A,i5)') 'exch  = ',exch
  write(6,'(A,i5)') 'npair = ',npair
  write(6,'(A,i5)') 'nneq  = ',nneq
  write(6,'(A,i5)') 'nmax  = ',nmax
  write(6,'(A,i5)') 'lmax  = ',lmax
  write(6,'(A,i5)') 'iPrint= ',iPrint
  write(6,*) 'AnisoLines1 = ',AnisoLines1
  write(6,*) 'AnisoLines3 = ',AnisoLines3
  write(6,*) 'AnisoLines9 = ',AnisoLines9
  write(6,*) 'Dipol = ',Dipol
  write(6,*) 'JITO  = ',JITO_exchange
  write(6,*) 'KE    = ',KE
  do i=1,npair
    write(6,'(A,i2,A,i3,3x,A,i2,A,i3)') 'i_pair(',i,',1) = ',i_pair(i,1),'i_pair(',i,',2) = ',i_pair(i,2)
  end do
  do i=1,nneq
    write(6,'(A,i2,A,i3)') '   neq(',i,') = ',neq(i)
  end do
  do i=1,nneq
    write(6,'(A,i2,A,i3)') ' nexch(',i,') = ',neq(i)
  end do
  do i=1,nneq
    write(6,'(A,i2,A,100F10.3)') '   eso(',i,') = ',(eso(i,j),j=1,nexch(i))
  end do
  call xFlush(6)
end if

if ((iPrint > 2) .or. DBG) then
  ! print the Exchange Hamiltonian:
  call pa_prham(exch,npair,i_pair,nneq,neq,nexch,nmax,lmax,eso,HLIN1,HLIN3,HLIN9,HDIP,HKEX,HDMO,HITO,Dipol,AnisoLines1, &
                AnisoLines3,AnisoLines9,KE,DM_exchange,JITO_exchange)
end if

! diagonalize the Hamiltonian:
call pa_diagham(exch,npair,i_pair,nneq,neq,nexch,nmax,lmax,eso,HLIN1,HLIN3,HLIN9,HDIP,HKEX,HDMO,HITO,Dipol,DM_exchange, &
                AnisoLines1,AnisoLines3,AnisoLines9,KE,JITO_exchange,WLIN1,WLIN3,WLIN9,WLIN,WDIP,WKEX,WDMO,WITO,W,Z)

! printout the resulting eigenstates:
call pa_preigen(exch,lmax,ibas,Dipol,AnisoLines1,AnisoLines3,AnisoLines9,KE,DM_exchange,JITO_exchange,WLIN,WDIP,WKEX,WDMO,WITO,W, &
                Z,iPrint)
!Z =  exchange eigenstates:
NmaxPop = 500
if (NmaxPop > exch) then
  NmaxPop = exch
end if
call PopAnalysis(nneq,neq,exch,nexch,nmax,lmax,NmaxPop,Z)
do i=exch,1,-1
  w(i) = w(i)-w(1)
end do

! some verification
if (dnrm2_(exch,WLIN,1) > 1.0d-13) call Add_Info('EXCHCTL::  WLIN',[dnrm2_(exch,WLIN,1)],1,8)
if (dnrm2_(exch,WDIP,1) > 1.0d-13) call Add_Info('EXCHCTL::  WDIP',[dnrm2_(exch,WDIP,1)],1,8)
if (dnrm2_(exch,WKEX,1) > 1.0d-13) call Add_Info('EXCHCTL::  WKEX',[dnrm2_(exch,WKEX,1)],1,8)
if (dnrm2_(exch,W,1) > 1.0d-13) call Add_Info('EXCHCTL::     W',[dnrm2_(exch,W,1)],1,8)
! compute the moments:
call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,M,1)
call zcopy_(3*exch*exch,[(0.0_wp,0.0_wp)],0,S,1)
if (dbg) then
  write(6,'(A)') 'Magnetic moments before the build of coupled M and S matrices'
  if (nPair > 0) then
    do lp=1,npair
      lb1 = i_pair(lp,1)
      lb2 = i_pair(lp,2)
      i1 = nind(lb1,1) ! indices of non-equivalent sites
      i2 = nind(lb2,1) ! indices of non-equivalent sites
      n1 = nexch(i1)
      n2 = nexch(i2)
      call prMom('EXCHCTL,Before M ans S, SM(i1):',SM(i1,1:3,1:n1,1:n1),n1)
      call prMom('EXCHCTL,Before M ans S, SM(i2):',SM(i2,1:3,1:n2,1:n2),n2)
      call prMom('EXCHCTL,Before M ans S, MM(i1):',MM(i1,1:3,1:n1,1:n1),n1)
      call prMom('EXCHCTL,Before M ans S, MM(i2):',MM(i2,1:3,1:n2,1:n2),n2)
    end do
  end if
end if

do L=1,3
  do isite=1,lmax
    do nb1=1,exch
      do lp=1,lmax
        icoord(lp) = ibas(nb1,lp)
      end do
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
  call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,TMP,1)
  call zgemm_('C','N',EXCH,EXCH,EXCH,(1.0_wp,0.0_wp),Z,EXCH,M(L,:,:),EXCH,(0.0_wp,0.0_wp),TMP,EXCH)
  call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,M(L,:,:),1)
  call zgemm_('N','N',EXCH,EXCH,EXCH,(1.0_wp,0.0_wp),TMP,EXCH,Z,EXCH,(0.0_wp,0.0_wp),M(L,:,:),EXCH)
  call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,TMP,1)
  ! spin moment
  call zgemm_('C','N',EXCH,EXCH,EXCH,(1.0_wp,0.0_wp),Z,EXCH,S(L,:,:),EXCH,(0.0_wp,0.0_wp),TMP,EXCH)
  call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,S(L,:,:),1)
  call zgemm_('N','N',EXCH,EXCH,EXCH,(1.0_wp,0.0_wp),TMP,EXCH,Z,EXCH,(0.0_wp,0.0_wp),S(L,:,:),EXCH)
end do  ! L

if (npair > 0) then
  ! ITO decomposition of the exchange and dipolar interactions:
  call pr_ito_int(npair,i_pair,lmax,nexch,nneq,neqv,itype,neq,nmax,eso(1:nneq,1:nmax),MM(1:nneq,1:3,1:nmax,1:nmax), &
                  SM(1:nneq,1:3,1:nmax,1:nmax),rot,Dipol,AnisoLines1,AnisoLines3,AnisoLines9,DM_exchange,JITO_exchange,HLIN1, &
                  HLIN3,HLIN9,HDIP,HDMO,HITO)
end if

! projection on the Ising Hamiltonian
! accessible format
!write(6,'(/)')
!write(6,'(A)') 'Complete decomposition of the exchange interaction'
!write(6,'(A)')

!-----------------------------------------------------------------------
! deallocate memory for this function:
if (lmax >= 0) then
  ! exchange energy spectrum
  call mma_deallocate(intc)
  call mma_deallocate(icoord)
  call mma_deallocate(nind)
  if (exch >= 0) then
    call mma_deallocate(ibas)
  end if
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
  call mma_deallocate(ZA1)
  call mma_deallocate(ZA2)
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

if (exch >= 0) then
  call mma_deallocate(tmp)
end if

if (nneq >= 0) then
  call mma_deallocate(nexchR)
  call mma_deallocate(SMR)
  call mma_deallocate(MMR)
  if (neqv >= 0) then
    call mma_deallocate(rotR)
  end if
end if

if (exchR >= 0) then
  if (lmax >= 0) then
    call mma_deallocate(ibasR)
  end if
  call mma_deallocate(WR)
  call mma_deallocate(ZR)
  call mma_deallocate(MR)
  call mma_deallocate(SR)
end if

if (npair >= 0) then
  call mma_deallocate(HKEXR)
end if

if (lmax >= 0) then
  call mma_deallocate(intcR)
end if

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
!ZZR = 0.0_wp
!ZZI = 0.0_wp
!do i=1,EXCH
!  do j=1,EXCH
!    read(77,'(2ES24.14)') ZZR(i,j),ZZI(i,j)
!  end do
!end do
!close(77)
!ZF = (0.0_wp,0.0_wp)
!do i=1,EXCH
!  do j=1,EXCH
!    ZF(i,j) = cmplx(ZZR(i,j),ZZI(i,j),8)
!  end do
!end do
!OVLP = (0.0_wp,0.0_wp)
!call ZGEMM_('C','N',EXCH,EXCH,EXCH,(1.0_wp,0.0_wp),Z(1:EXCH,1:EXCH),EXCH,ZF(1:EXCH,1:EXCH),EXCH,(0.0_wp,0.0_wp), &
!            OVLP(1:EXCH,1:EXCH),EXCH )
!do i=1,1
!  do j=1,EXCH
!    if (ABS(OVLP(j,i)) > 1.0e-7_wp) &
!      write(6,'(A,i3,A,i3,A,2F18.14,3x,A,F18.14,3x,A,F18.14)') 'OVLP-1(',j,',',i,') = ',OVLP(j,i),'ABS = ',ABS(OVLP(j,i)), &
!                                                               'ABS^2 = ',ABS(OVLP(j,i))*ABS(OVLP(j,i))
!  end do
!end do
!OVLP = (0.0_wp,0.0_wp)
!call ZGEMM_('C','N',EXCH,EXCH,EXCH,(1.0_wp,0.0_wp),ZF(1:EXCH,1:EXCH),EXCH,Z(1:EXCH,1:EXCH),EXCH,(0.0_wp,0.0_wp), &
!            OVLP(1:EXCH,1:EXCH),EXCH )
!do i=1,1
!  do j=1,EXCH
!    if (ABS(OVLP(j,i)) > 1.0e-7_wp) &
!      write(6,'(A,i3,A,i3,A,2F18.14,3x,A,F18.14,3x,A,F18.14)') 'OVLP-2(',i,',',j,') = ',OVLP(i,j),'ABS = ',ABS(OVLP(i,j)), &
!                                                                'ABS^2 = ',ABS(OVLP(i,j))*ABS(OVLP(i,j))
!    end if
!  end do
!end do
!199 continue

return

end subroutine exchctl
