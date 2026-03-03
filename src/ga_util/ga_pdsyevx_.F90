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
!     nprow/ npcol  - number of processor rows/cols in virtual proc. grid  *
!     myrow/mycol   - cordinates of my processor in virtual proc. grid     *
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

implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack_.fh"
integer g_a               ! matrix A
integer g_b               ! matrix B
integer nb8               ! block size
real*8 eval(*)
character*1 jobz, range, uplo
real*8, allocatable :: adrA(:)          ! A
real*8, allocatable :: adrB(:)          ! B
logical oactive           ! true iff this process participates
integer dimA18, dimA28, typeA
integer dimB18, dimB28, typeB
SCALAPACKINT dimA1,dimA2
SCALAPACKINT dimB1,dimB2
SCALAPACKINT mpA,nqA     ! rows/cols of A held by the processor
SCALAPACKINT mpB,nqB     ! rows/cols of B held by the processor
integer me
SCALAPACKINT lda,ldb
integer elemA, elemB
SCALAPACKINT numroc
SCALAPACKINT nb                ! block size
SCALAPACKINT descA(9),descB(9) !descriptor for scalapack
integer ngaps
real*8, allocatable :: adrgaps(:)
integer iclu
integer, allocatable :: adrclustr(:)
integer, allocatable :: adrfail(:)
integer liwork
integer, allocatable :: adriwork(:)
integer lcwork
real*8, allocatable :: adrcwork(:)
SCALAPACKINT lcwork4
SCALAPACKINT liwork4
SCALAPACKINT nn,mq0,np0
real*8 vl, vu, abstol, orfac
SCALAPACKINT il,iu
SCALAPACKINT m,nz
SCALAPACKINT n
SCALAPACKINT info
integer info8
SCALAPACKINT zero4,one4,two4,four4
integer two4n
parameter(zero4=0,one4=1,two4=2,four4=4)
real*8 pdlamch
external pdlamch
SCALAPACKINT iceil

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

  mpA = numroc(dimA1,nb,myrow2,zero4,nprow2)
  nqA = numroc(dimA2,nb,mycol2,zero4,npcol2)

  mpB = numroc(dimB1,nb,myrow2,zero4,nprow2)
  nqB = numroc(dimB2,nb,mycol2,zero4,npcol2)

  lda = max(one4,mpA)
  ldb = max(one4,mpB)

  ! let scalapack check for errors

  ! should check to see if this is a compute node
  ! check to see how this works in the new data server model

  elemA = mpA*nqA
  if (elemA /= 0) call mma_allocate(adrA,elemA,Label='adrA')

  !*** copy g_a to A using the block cyclic scalapack format

  call ga_to_SL2_(g_a,dimA1,dimA2,nb,nb,adrA,lda)

  elemB = mpB*nqB

  if (elemB /= 0) call mma_allocate(adrB,elemB,Label='adrB')

  ngaps = nprow2*npcol2
  if (ngaps /= 0) call mma_allocate(adrgaps,ngaps,Label='adrgaps')

  iclu = 2*nprow2*npcol2
  two4n = two4*n
  iclu = max(two4n,iclu)
  if (iclu /= 0) call mma_allocate(adrclustr,iclu,Label='adrclustr')
  call mma_allocate(adrfail,dima18,Label='adrfail')
end if
call ga_sync()
if (oactive) then

  !*** fill SCALAPACK matrix descriptors

  call descinit(descA,dimA1,dimA2,nb,nb,zero4,zero4,islctxt2,lda,info)
  info8 = info
  if (info8 /= 0) call ga_error(' ga_pdsyevx: descinit A failed ',-info8)
  call descinit(descB,dimB1,dimB2,nb,nb,zero4,zero4,islctxt2,ldb,info)
  info8 = info
  if (info8 /= 0) call ga_error(' ga_pdsyevx: descinit B failed ',-info8)

  jobz = 'V'
  range = 'A'
  range = 'I'
  uplo = 'L'
  vl = 0.d0
  vu = 0.d0
  il = 0
  iu = 0
  il = 1
  iu = n
  nz = 0

  ! ability to deal with orthonormality ; let's just
  ! have the regular scalapack stuff for the moment

  liwork = 6*max(n,nprow2*npcol2+one4,four4)
  liwork = liwork
  if (liwork /= 0) call mma_allocate(adriwork,liwork,Label='adriwork')

  nn = max(n,nb,two4)
  np0 = numroc(nn,nb,zero4,zero4,nprow2)
  mq0 = numroc(nn,nb,zero4,zero4,npcol2)

  orfac = 1.d-3

  lcwork = 5*n+max(5*NN,(NP0*MQ0+2*nb*nb))+ICEIL(N,NPROW2*NPCOL2)*NN+1
  if (lcwork /= 0) call mma_allocate(adrcwork,lcwork,Label='adrcwork')

  abstol = pdlamch(islctxt2,'U')

  liwork4 = liwork
  lcwork4 = lcwork
  call pdsyevx(jobz,range,uplo,n,adrA,one4,one4,descA,vl,vu,il,iu,abstol,m,nz,eval,orfac,adrB,one4,one4,descB,adrcwork,lcwork4, &
               adriwork,liwork4,adrfail,adrclustr,adrgaps,info)

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

implicit none
#include "scalapack_.fh"
#include "global.fh"
SCALAPACKINT n
SCALAPACKINT zero4
parameter(zero4=0)
SCALAPACKINT slgetmxproc_
external slgetmxproc_

if (init2) return

!**** call ga_sync before to enter in BLACS and after
call ga_sync()

call blacs_pinfo(iam,nnodes)

! determine optimal nprocs for eigensolvers based on matrix size n

maxproc = slgetmxproc_(n,nnodes)
call FindGrid_(maxproc,nprow2,npcol2)

call blacs_get(zero4,zero4,islctxt2)
call blacs_gridinit(islctxt2,'R',nprow2,npcol2)

if (iam < maxproc) then
  call blacs_gridinfo(iSLctxt2,nprow2,npcol2,myrow2,mycol2)
else
  nprow2 = 0
  npcol2 = 0
  myrow2 = 0
  mycol2 = 0
end if
init2 = .true.

call ga_sync()

end subroutine slinit2_
!=!=
SCALAPACKINT function slgetmxproc_(n,nnodes)

implicit none
SCALAPACKINT nnodes
SCALAPACKINT n
SCALAPACKINT i
real*8 fact
SCALAPACKINT nmax,twoi
real*8 nprocs0
real*8 otto
parameter(nmax=11,fact=((7108d0*7108d0)/1024d0),otto=8d0)
!new parameter(fact=((7108d0*7108d0)/512d0))

! lower bound of 8 procs
nprocs0 = max((n*n)/fact,otto)

! try to get powers of two

do i=nmax,0,-1
  if (nint(nprocs0/(2d0**i)) == 1) goto 1
end do
i = 4
1 twoi = 2**i
slgetmxproc_ = min(nnodes,twoi)

end function slgetmxproc_
!=!=
subroutine FindGrid_(nnodes,nprow,npcol)
!*** determine nprow, npcol from nnodes
!*** solution is searched in the neighborhood of the square grid

implicit none
SCALAPACKINT nnodes,nprow,npcol,i

! try to get the 1:4 ratio

npcol = 2*int(sqrt(dble(nnodes)))
do i=npcol,1,-1
  if (mod(nnodes,i) == 0) goto 1
end do
1 continue
npcol = i
nprow = nnodes/npcol
if (nprow > npcol) then
  i = npcol
  npcol = nprow
  nprow = i
end if
end subroutine FindGrid_
!=!=
subroutine ga_to_SL2_(g_a,dim1,dim2,nbr,nbc,s_a,lda4)
!*** transforms a GA to SL format
!*** reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'
!
! g_a         - handle for Global Array that is being converted SL format
! dim1, dim2  - dimensions of Global Array
! nbr         - number of block rows?
! nbc         - number of block columns?
! s_a         - local array holding SL formatted data
! lda4        - leading dimension of s_a

implicit none
#include "scalapack_.fh"
#include "global.fh"
SCALAPACKINT nbr,nbc,lda4,dim1,dim2
integer g_a
real*8 s_a(lda4,*)
integer row, col, tcol, trow, rbase, cbase
SCALAPACKINT pcol,prow
integer lda
integer r0i, r1i, r0im1, r1im1, rowl, coll
integer marg1, marg2
logical putpending

!***  Synchronize at the beginning

lda = lda4
rbase = 1
cbase = 1

Trow = nbr*nprow2
Tcol = nbc*npcol2
do col=1,dim2,nbc
  ! processor column that holds "col"
  pcol = mod(col,Tcol)/nbc
  if (mycol2 == pcol) then
    r0im1 = -9999
    r1im1 = 0
    putpending = .false.
    marg1 = col+nbc-1
    marg2 = dim2
    coll = min(marg1,marg2)

    do row=1,dim1,nbr
      ! processor row that holds "row"
      prow = mod(row,Trow)/nbr
      if (myrow2 == prow) then
        if (.not. putpending) then
          r0im1 = row
          marg1 = row+nbr-1
          marg2 = dim1
          r1im1 = min(marg1,marg2)
        end if
        r0i = row
        marg1 = row+nbr-1
        marg2 = dim1
        r1i = min(marg1,marg2)

        if ((r0i == r1im1+1) .and. (r1i /= dim1)) then
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
subroutine ga_from_SL2_(g_a,dim1,dim2,nbr,nbc,s_a,lda4)
!*** transforms a matrix from SL to GA format
!*** reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'

implicit none
#include "scalapack_.fh"
#include "global.fh"
integer g_a, lda
SCALAPACKINT nbr,nbc,lda4,dim1,dim2
real*8 s_a(lda4,*)
integer tcol, trow
integer rbase, cbase
SCALAPACKINT pcol,prow
integer row, col, rowl, coll
integer r0i, r1i, r0im1, r1im1
integer marg1, marg2
logical putpending

lda = lda4

!**** Syncronize at the beginning

rbase = 1
cbase = 1

Trow = nbr*nprow2
Tcol = nbc*npcol2

do col=1,dim2,nbc
  ! processor column that holds "col"
  pcol = mod(col,Tcol)/nbc
  if (mycol2 == pcol) then
    r0im1 = -9999
    r1im1 = 0
    putpending = .false.
    marg1 = col+nbc-1
    marg2 = dim2
    coll = min(marg1,marg2)
    do row=1,dim1,nbr
      ! processor row that holds "row"
      prow = mod(row,Trow)/nbr
      if (myrow2 == prow) then
        if (.not. putpending) then
          r0im1 = row
          marg1 = row+nbr-1
          marg2 = dim1
          r1im1 = min(marg1,marg2)
        end if
        r0i = row
        marg1 = row+nbr-1
        marg2 = dim1
        r1i = min(marg1,marg2)

        if ((r0i == r1im1+1) .and. (r1i /= dim1)) then
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
