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
!***************************************************************************
! GA to/from ScaLAPACK (square block scattered decomposition) interface    *
!                                                                          *
! common variables:                                                        *
!     nnodes        - number of processors                                 *
!     iam           - my processor number                                  *
!     nprow / npcol - number of processor rows/cols in virtual proc. grid  *
!     myrow / mycol - cordinates of my processor in virtual proc. grid     *
!                                                                          *
!***************************************************************************
!* 04/12/96  GVT  Changed the code to adapt to a new version of ScaLAPACK  *
!*           Giuseppe Vitillaro peppe@unipg.it                             *
!* 11/12/14  Slightly modified for the DGA interface                       *
!*           Victor Vysotskiy                                              *
!***************************************************************************

#include "compiler_features.h"
#if defined (_SCALAPACK_) && defined (_MOLCAS_MPP_)

subroutine ga_pdsyevx_(g_a,g_b,eval,nb8)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, BLASInt

implicit none
integer(kind=iwp) :: g_a, g_b, nb8
real(kind=wp) :: eval(*)
integer(kind=iwp) :: dimA18, dimA28, dimB18, dimB28, elemA, elemB, iclu, info8, lcwork, liwork, me, ngaps, two4n, typeA, typeB
integer(kind=BLASInt) :: descA(9), descB(9), dimA1, dimA2, dimB1, dimB2, iceil, il, info, iu, lcwork4, lda, ldb, liwork4, m, mpA, &
                         mpB, mq0, n, nb, nn, np0, nqA, nqB, numroc, nz
real(kind=wp) :: abstol, orfac, vl, vu
logical(kind=iwp) :: oactive
character :: jobz, rng, uplo
integer(kind=iwp), allocatable :: adrclustr(:), adrfail(:), adriwork(:)
real(kind=wp), allocatable :: adrA(:), adrB(:), adrcwork(:), adrgaps(:)
real(kind=wp), external :: pdlamch
#include "scalapack_.fh"
#include "global.fh"
#include "mafdecls.fh"

! oactive     : true iff this process participates
! mpA,nqA     : rows/cols of A held by the processor
! mpB,nqB     : rows/cols of B held by the processor
! nb          : block size
! descA,descB : descriptors for scalapack

! processor dependent; machine dependent

if (nb8 == 0) then
  nb = 16
else
  nb = nb8
end if

!*** check environment

me = ga_nodeid()

!*** check GA info for input arrays

call ga_inquire(g_a,typeA,dimA18,dimA28)
call ga_inquire(g_b,typeB,dimB18,dimB28)

dima1 = dima18
dima2 = dima28
dimb1 = dimb18
dimb2 = dimb28
n = dima1
if (nb < 1) nb = 1
if (dimA1 /= dima2) call ga_error('ga_pdsyevx: matrix A not square ',0)
if (dimb1 /= dimb2) call ga_error('ga_pdsyevx: matrix B not square ',0)
if (dimb1 /= n) call ga_error('ga_pdsyevx: size matrix A and B differ ',0)

!*** initialize SL interface

call SLinit2_(n)
oactive = iam < maxproc
call ga_sync
if (oactive) then

!*** find SBS format parameters

  mpA = numroc(dimA1,nb,myrow,0_BLASInt,nprow)
  nqA = numroc(dimA2,nb,mycol,0_BLASInt,npcol)

  mpB = numroc(dimB1,nb,myrow,0_BLASInt,nprow)
  nqB = numroc(dimB2,nb,mycol,0_BLASInt,npcol)

  lda = max(1_BLASInt,mpA)
  ldb = max(1_BLASInt,mpB)

  ! let scalapack check for errors

  ! should check to see if this is a compute node
  ! check to see how this works in the new data server model

  elemA = mpA*nqA
  if (elemA /= 0) call mma_allocate(adrA,elemA,Label='adrA')

  !*** copy g_a to A using the block cyclic scalapack format

  call ga_to_SL2_(g_a,dimA1,dimA2,nb,nb,adrA,lda)

  elemB = mpB*nqB

  if (elemB /= 0) call mma_allocate(adrB,elemB,Label='adrB')

  ngaps = nprow*npcol
  if (ngaps /= 0) call mma_allocate(adrgaps,ngaps,Label='adrgaps')

  iclu = 2*nprow*npcol
  two4n = 2_BLASInt*n
  iclu = max(two4n,iclu)
  if (iclu /= 0) call mma_allocate(adrclustr,iclu,Label='adrclustr')
  call mma_allocate(adrfail,dima18,Label='adrfail')
end if
call ga_sync()
if (oactive) then

  !*** fill SCALAPACK matrix descriptors

  call descinit(descA,dimA1,dimA2,nb,nb,0_BLASInt,0_BLASInt,islctxt2,lda,info)
  info8 = info
  if (info8 /= 0) call ga_error(' ga_pdsyevx: descinit A failed ',-info8)
  call descinit(descB,dimB1,dimB2,nb,nb,0_BLASInt,0_BLASInt,islctxt2,ldb,info)
  info8 = info
  if (info8 /= 0) call ga_error(' ga_pdsyevx: descinit B failed ',-info8)

  jobz = 'V'
  rng = 'I'
  uplo = 'L'
  vl = Zero
  vu = Zero
  il = 0
  iu = 0
  il = 1
  iu = n
  nz = 0

  ! ability to deal with orthonormality ; let's just
  ! have the regular scalapack stuff for the moment

  liwork = 6*max(n,nprow*npcol+1_BLASInt,4_BLASInt)
  liwork = liwork
  if (liwork /= 0) call mma_allocate(adriwork,liwork,Label='adriwork')

  nn = max(n,nb,2_BLASInt)
  np0 = numroc(nn,nb,0_BLASInt,0_BLASInt,nprow)
  mq0 = numroc(nn,nb,0_BLASInt,0_BLASInt,npcol)

  orfac = 1.0e-3_wp

  lcwork = 5*n+max(5*NN,(NP0*MQ0+2*nb*nb))+ICEIL(N,NPROW*NPCOL)*NN+1
  if (lcwork /= 0) call mma_allocate(adrcwork,lcwork,Label='adrcwork')

  abstol = pdlamch(islctxt2,'U')

  liwork4 = liwork
  lcwork4 = lcwork
  call pdsyevx(jobz,rng,uplo,n,adrA,1_BLASInt,1_BLASInt,descA,vl,vu,il,iu,abstol,m,nz,eval,orfac,adrB,1_BLASInt,1_BLASInt,descB, &
               adrcwork,lcwork4,adriwork,liwork4,adrfail,adrclustr,adrgaps,info)

  if (nz /= n) then
    if (info /= 0) then
      if (info > 0) then
        call ga_error(' ga_pdsyevx: argument is illegal ',info)
      else
        call ga_error(' ga_pdsyevx: eigenvectors failed to converge ',info)
      end if
    end if
  end if

  !*** copy solution matrix back to g_c

  call ga_from_SL2_(g_b,dimA1,dimB2,nb,nb,adrB,ldb)

  !*** deallocate work/SL arrays

  if (lcwork /= 0) call mma_deallocate(adrcwork)
  if (liwork /= 0) call mma_deallocate(adriwork)
  if (iclu /= 0) call mma_deallocate(adrclustr)
  if (ngaps /= 0) call mma_deallocate(adrgaps)
  if (elemB /= 0) call mma_deallocate(adrB)
  if (elemA /= 0) call mma_deallocate(adrA)
  call mma_deallocate(adrfail)
end if

call ga_sync()
! broadcast evals
if ((maxproc < nnodes) .or. (dima1 <= nb)) call GA_Brdcst(MT_DBL,eval,dima18,0)

end subroutine ga_pdsyevx_
!=!=
subroutine slinit2_(n)

use Definitions, only: iwp, BLASInt

implicit none
integer(kind=BLASInt) :: n
logical(kind=iwp) :: initialized = .false.
integer(kind=BLASInt), external :: slgetmxproc_
#include "scalapack_.fh"

if (initialized) return

!**** call ga_sync before to enter in BLACS and after
call ga_sync()

call blacs_pinfo(iam,nnodes)

! determine optimal nprocs for eigensolvers based on matrix size n

maxproc = slgetmxproc_(n,nnodes)
call FindGrid_(maxproc,nprow,npcol)

call blacs_get(0_BLASInt,0_BLASInt,islctxt2)
call blacs_gridinit(islctxt2,'R',nprow,npcol)

if (iam < maxproc) then
  call blacs_gridinfo(iSLctxt2,nprow,npcol,myrow,mycol)
else
  nprow = 0
  npcol = 0
  myrow = 0
  mycol = 0
end if
initialized = .true.

call ga_sync()

end subroutine slinit2_
!=!=
function slgetmxproc_(n,nnodes)

use Constants, only: Two, Eight
use Definitions, only: wp, iwp, BLASInt

implicit none
integer(kind=BLASInt) :: slgetmxproc_
integer(kind=BLASInt) :: n
integer(kind=BLASInt) :: i, nnodes, twoi
real(kind=wp) :: nprocs0
integer(kind=iwp), parameter :: nmax = 11
real(kind=wp), parameter :: fact = 7108.0_wp**2/1024.0_wp
!new real(kind=wp), parameter :: fact = 7108.0_wp**2/512.0_wp

! lower bound of 8 procs
nprocs0 = max(n**2/fact,Eight)

! try to get powers of two

do i=nmax,0,-1
  if (nint(nprocs0/(Two**i)) == 1) goto 1
end do
i = 4
1 twoi = 2**i
slgetmxproc_ = min(nnodes,twoi)

end function slgetmxproc_
!=!=
subroutine FindGrid_(nnodes,npr,npc)
!*** determine npr, npc from nnodes
!*** solution is searched in the neighborhood of the square grid

use Definitions, only: wp, BLASInt

implicit none
integer(kind=BLASInt) :: nnodes, npr, npc
integer(kind=BLASInt) :: i

! try to get the 1:4 ratio

npc = 2*int(sqrt(real(nnodes,kind=wp)))
do i=npc,1,-1
  if (mod(nnodes,i) == 0) goto 1
end do
1 continue
npc = i
npr = nnodes/npc
if (npr > npc) then
  i = npc
  npc = npr
  npr = i
end if

end subroutine FindGrid_
!=!=
subroutine ga_to_SL2_(g_a,dim_1,dim_2,nbr,nbc,s_a,lda4)
!*** transforms a GA to SL format
!*** reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'
!
! g_a          - handle for Global Array that is being converted SL format
! dim_1, dim_2 - dimensions of Global Array
! nbr          - number of block rows?
! nbc          - number of block columns?
! s_a          - local array holding SL formatted data
! lda4         - leading dimension of s_a

use Definitions, only: wp, iwp, BLASInt

implicit none
integer(kind=iwp) :: g_a
integer(kind=BLASInt) :: dim_1, dim_2, nbr, nbc, lda4
real(kind=wp) :: s_a(lda4,*)
integer(kind=iwp) :: cbase, col, coll, lda, marg1, marg2, r0i, r0im1, r1i, r1im1, rbase, row, rowl, tcol, trow
integer(kind=BLASInt) :: pcol, prow
logical(kind=iwp) :: putpending
#include "scalapack_.fh"

!***  Synchronize at the beginning

lda = lda4
rbase = 1
cbase = 1

Trow = nbr*nprow
Tcol = nbc*npcol

do col=1,dim_2,nbc
  ! processor column that holds "col"
  pcol = mod(col,Tcol)/nbc
  if (mycol == pcol) then
    r0im1 = -9999
    r1im1 = 0
    putpending = .false.
    marg1 = col+nbc-1
    marg2 = dim_2
    coll = min(marg1,marg2)

    do row=1,dim_1,nbr
      ! processor row that holds "row"
      prow = mod(row,Trow)/nbr
      if (myrow == prow) then
        if (.not. putpending) then
          r0im1 = row
          marg1 = row+nbr-1
          marg2 = dim_1
          r1im1 = min(marg1,marg2)
        end if
        r0i = row
        marg1 = row+nbr-1
        marg2 = dim_1
        r1i = min(marg1,marg2)

        if ((r0i == r1im1+1) .and. (r1i /= dim_1)) then
          r1im1 = r1i
          putpending = .true.
        else
          rowl = r1i
          call ga_get(g_a,r0im1,rowl,col,coll,s_a(rbase,cbase),lda)
          putpending = .false.
          rbase = rbase+rowl-r0im1+1
          r0im1 = -1
          r1im1 = -1
        end if
      end if
    end do
    rbase = 1
    cbase = cbase+nbc
  end if
end do

!**** ... and at the end
end subroutine ga_to_SL2_
!=!=
subroutine ga_from_SL2_(g_a,dim_1,dim_2,nbr,nbc,s_a,lda4)
!*** transforms a matrix from SL to GA format
!*** reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'

use Definitions, only: wp, iwp, BLASInt

implicit none
integer(kind=iwp) :: g_a
integer(kind=BLASInt) :: dim_1, dim_2, nbr, nbc, lda4
real(kind=wp) :: s_a(lda4,*)
integer(kind=iwp) :: cbase, col, coll, lda, marg1, marg2, r0i, r0im1, r1i, r1im1, rbase, row, rowl, tcol, trow
integer(kind=BLASInt) :: pcol, prow
logical(kind=iwp) :: putpending
#include "scalapack_.fh"

!**** Syncronize at the beginning

lda = lda4
rbase = 1
cbase = 1

Trow = nbr*nprow
Tcol = nbc*npcol

do col=1,dim_2,nbc
  ! processor column that holds "col"
  pcol = mod(col,Tcol)/nbc
  if (mycol == pcol) then
    r0im1 = -9999
    r1im1 = 0
    putpending = .false.
    marg1 = col+nbc-1
    marg2 = dim_2
    coll = min(marg1,marg2)
    do row=1,dim_1,nbr
      ! processor row that holds "row"
      prow = mod(row,Trow)/nbr
      if (myrow == prow) then
        if (.not. putpending) then
          r0im1 = row
          marg1 = row+nbr-1
          marg2 = dim_1
          r1im1 = min(marg1,marg2)
        end if
        r0i = row
        marg1 = row+nbr-1
        marg2 = dim_1
        r1i = min(marg1,marg2)

        if ((r0i == r1im1+1) .and. (r1i /= dim_1)) then
          r1im1 = r1i
          putpending = .true.
        else
          rowl = r1i
          call ga_put(g_a,r0im1,rowl,col,coll,s_a(rbase,cbase),lda)
          putpending = .false.
          rbase = rbase+rowl-r0im1+1
          r0im1 = -1
          r1im1 = -1
        end if
      end if
    end do
    rbase = 1
    cbase = cbase+nbc
  end if
end do

!**** ... and at the end
end subroutine ga_from_SL2_

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(ga_pdsyevx)

#endif
