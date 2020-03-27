************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
****************************************************************************
* GA to/from ScaLAPACK (square block scattered decomposition) interface    *
*                                                                          *
* common variables:                                                        *
*     nnodes        - number of processors                                 *
*     iam           - my processor number                                  *
*     nprow/ npcol  - number of processor rows/cols in virtual proc. grid  *
*     myrow/mycol   - cordinates of my processor in virtual proc. grid     *
*                                                                          *
c***************************************************************************
c* 04/12/96  GVT  Changed the code to adapt to a new version of ScaLAPACK  *
c*           Giuseppe Vitillaro peppe@unipg.it                             *
c* 11/12/14  Slightly modified for the DGA interface                       *
c*           Victor Vysotskiy                                              *
c***************************************************************************
#if defined(SCALAPACK) && defined (_MOLCAS_MPP_)
      subroutine ga_pdsyevx_(g_a, g_b, eval, nb8)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "scalapack_.fh"

      integer g_a               ! matrix A
      integer g_b               ! matrix B
      integer nb8               ! block size
      real*8 eval(*)

      character*1 jobz, range, uplo
c
      integer adra          ! A
      integer adrb          ! B
c
c
      logical oactive           ! true iff this process participates
      integer dimA18, dimA28, typeA
      integer dimB18, dimB28, typeB
      SCALAPACKINT dimA1, dimA2
      SCALAPACKINT dimB1, dimB2
c
      SCALAPACKINT mpA, nqA     ! rows/cols of A held by the processor
      SCALAPACKINT mpB, nqB     ! rows/cols of B held by the processor
c
      integer me
      SCALAPACKINT lda, ldb
      integer elemA,elemB
      SCALAPACKINT numroc

      SCALAPACKINT nb                ! block size
      SCALAPACKINT descA(9), descB(9) !descriptor for scalapack

c
      integer ngaps, adrgaps
      integer iclu, adrclustr
      integer adrfail
      integer liwork, adriwork
      integer lcwork, adrcwork
      SCALAPACKINT lcwork4
      SCALAPACKINT liwork4
c
      SCALAPACKINT nn,mq0, np0
      real*8 vl, vu, abstol, orfac
      SCALAPACKINT il, iu
      SCALAPACKINT m, nz
      SCALAPACKINT n
      SCALAPACKINT info
      integer info8
      SCALAPACKINT zero4,one4,two4,four4
      integer two4n
      parameter(zero4=0,one4=1,two4=2,four4=4)

      real*8 pdlamch
      external pdlamch
      SCALAPACKINT iceil

#include "WrkSpc.fh"

c
c     processor dependent; machine dependent
c
      if(nb8.eq.0) then
         nb=16
      else
         nb=nb8
      endif
c
c***  check environment
c
      me     = ga_nodeid()
c
c***  check GA info for input arrays
c
      call ga_inquire(g_a, typeA, dimA18, dimA28)
      call ga_inquire(g_b, typeB, dimB18, dimB28)

      dima1=dima18
      dima2=dima28
      dimb1=dimb18
      dimb2=dimb28
      n=dima1
      if(nb.lt.1) nb=1
      if(dimA1.ne.dima2) call ga_error(
     $     'ga_pdsyevx: matrix A not square ',0)
      if(dimb1.ne.dimb2) call ga_error(
     $     'ga_pdsyevx: matrix B not square ',0)
      if(dimb1.ne.n) call ga_error(
     $     'ga_pdsyevx: size matrix A and B differ ',0)


c
c
c***  initialize SL interface
c
      call SLinit2_(n)
      oactive=iam.lt.maxproc
      call ga_sync
      if (oactive) then
c
c***  find SBS format parameters
c
c
         mpA = numroc(dimA1, nb, myrow2, zero4, nprow2)
         nqA = numroc(dimA2, nb, mycol2, zero4, npcol2)
c
         mpB = numroc(dimB1, nb, myrow2, zero4, nprow2)
         nqB = numroc(dimB2, nb, mycol2, zero4, npcol2)
c
c
         lda = max(one4,mpA)
         ldb = max(one4,mpB)
c
c
c     let scalapack check for errors
c
c     should check to see if this is a compute node
c     check to see how this works in the new data server model
c
c
         elemA= mpA*nqA
         if(elemA.ne.0)
     $   Call GetMem('A','ALLO','REAL',adrA,elemA)
c
c***  copy g_a to A using the block cyclic scalapack format
c

         call ga_to_SL2_(g_a, dimA1, dimA2, nb, nb,
     $        WORK(adrA), lda)
c
         elemB= mpB*nqB

         if(elemB.ne.0)
     $   Call GetMem('B','ALLO','REAL',adrB,elemB)

c
         ngaps = nprow2*npcol2
         if(ngaps.ne.0)
     $   Call GetMem('GAP','ALLO','REAL',adrgaps,ngaps)
c
         iclu = 2*nprow2*npcol2
         two4n=two4*n
         iclu = max(two4n,iclu)
         if(iclu.ne.0)
     $   Call GetMem('ICLUS','ALLO','INTE',adrclustr,iclu)
         Call GetMem('IFAIL','ALLO','INTE',adrfail,dima18)
      endif
      call ga_sync()
      if(oactive) then
c
c***  fill SCALAPACK matrix descriptors
c
         call descinit(descA, dimA1, dimA2, nb, nb, zero4, zero4,
     $        islctxt2, lda, info)
         info8=info
         if(info8.ne.0) call ga_error(' ga_pdsyevx: descinit A failed ',
     $        -info8)
         call descinit(descB, dimB1, dimB2, nb, nb, zero4, zero4,
     $        islctxt2, ldb, info)
         info8=info
         if(info8.ne.0) call ga_error(' ga_pdsyevx: descinit B failed ',
     $        -info8)
c
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
c
c     ability to deal with orthonormality ; let's just
c     have the regular scalapack stuff for the moment
c
         liwork = 6*max(n, nprow2*npcol2+one4, four4)
         liwork=liwork
         if(liwork.ne.0)
     $   Call GetMem('IWORK','ALLO','INTE',adriwork,liwork)
c
         nn = max(n, nb, two4)
         np0 = numroc(nn, nb, zero4, zero4, nprow2)
         mq0 = numroc(nn, nb, zero4, zero4, npcol2)
c
c
         orfac = 1.d-3
c
c
         lcwork = 5*n +MAX(5*NN,(NP0*MQ0 + 2*nb*nb))+
     $          ICEIL( N, NPROW2*NPCOL2)*NN+1
         if(lcwork.ne.0)
     $   Call GetMem('CWORK','ALLO','REAL',adrcwork,lcwork)
c
c
         abstol=pdlamch(islctxt2, 'U')
c
c
         liwork4=liwork
         lcwork4=lcwork
         call pdsyevx(jobz, range, uplo,
     $        n, work(adrA), one4, one4, descA,vl,
     $        vu, il, iu, abstol, m, nz, eval, orfac, work(adrB),
     $        one4, one4, descB, work(adrcwork), lcwork4,
     $        iwork(adriwork), liwork4, iwork(adrfail),
     $        iwork(adrclustr), work(adrgaps), info)
c

         if (nz .ne. n ) then
            if ( info .ne. 0 ) then
               if ( info .gt. 0 ) then
         call ga_error(' ga_pdsyevx: argument is illegal ', info)
               else
         call ga_error(' ga_pdsyevx: eigenvectors failed to converge ',
     $                 info)
               endif
            endif
         endif
c
c
c
c***  copy solution matrix back to g_c
c
         call ga_from_SL2_(g_b, dimA1, dimB2,
     $        nb, nb, work(adrB),
     $        ldb)
c
c
c
c***  deallocate work/SL arrays
c
         if ( lcwork .ne. 0 )
     $   Call GetMem('CWORK','FREE','REAL',adrcwork,lcwork)
         if ( liwork .ne. 0 )
     $   Call GetMem('IWORK','FREE','INTE',adriwork,liwork)
         if ( iclu .ne. 0 )
     $   Call GetMem('ICLUS','FREE','INTE',adrclustr,iclu)
         if ( ngaps.ne.0 )
     $   Call GetMem('GAP','FREE','REAL',adrgaps,ngaps)
         if ( elemB .ne. 0 )
     $   Call GetMem('B','FREE','REAL',adrB,elemB)
         if ( elemA .ne. 0 )
     $   Call GetMem('A','FREE','REAL',adrA,elemA)
         Call GetMem('IFAIL','FREE','INTE',adrfail,dimA18)
      endif
c
      call ga_sync()
      if(maxproc.lt.nnodes.or.dima1.le.nb) then
c     broadcast evals
        Call GA_Brdcst(MT_DBL, eval, dima18, 0)
      endif
      return
      end

      subroutine slinit2_(n)
      implicit none
#include "scalapack_.fh"
#include "global.fh"
      SCALAPACKINT n
c
      SCALAPACKINT zero4
      parameter(zero4=0)
      SCALAPACKINT slgetmxproc_
      external slgetmxproc_
c
      if(init2)return
c
c**** call ga_sync before to enter in BLACS and after
      call ga_sync()
c
      call blacs_pinfo(iam, nnodes)
c
c     determine optimal nprocs for eigensolvers based on matrix size n
c
      maxproc=slgetmxproc_(n,nnodes)
      call FindGrid_(maxproc, nprow2, npcol2)
c
      call blacs_get( zero4, zero4, islctxt2 )
      call blacs_gridinit(islctxt2, 'R', nprow2, npcol2)


c
      if(iam.lt.maxproc) then
         call blacs_gridinfo(iSLctxt2, nprow2, npcol2, myrow2, mycol2)
      else
         nprow2=0
         npcol2=0
         myrow2=0
         mycol2=0
      endif
      init2=.true.
c
      call ga_sync()

c
      end
c
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
cnew      parameter(nmax=11,fact=((7108d0*7108d0)/512d0),otto=8d0)
c     lower bound of 8 procs
      nprocs0=max((n*n)/fact,otto)
c
c     try to get powers of two
c
      do i = nmax, 0, -1
         if(nint(nprocs0/(2d0**i)).eq.1) goto 1
      enddo
      i=4
1     twoi=2**i
      slgetmxproc_=min(nnodes,twoi)
      return
      end

      subroutine FindGrid_(nnodes, nprow, npcol)
c
c***  determine nprow, npcol from nnodes
c***  solution is searched in the neighborhood of the square grid
c
      implicit none
      SCALAPACKINT nnodes, nprow, npcol,i
c
c     try to get the 1:4 ratio
c
      npcol = 2*int(sqrt(dble(nnodes)))
      do i = npcol, 1, -1
         if(mod(nnodes,i).eq.0) goto 1
      enddo
1     continue
      npcol = i
      nprow = nnodes/npcol
      if(nprow.gt.npcol) then
         i=npcol
         npcol=nprow
         nprow=i
      endif
      end

      subroutine ga_to_SL2_(g_a, dim1, dim2, nbr, nbc, s_a, lda4)
c
c***  transforms a GA to SL format
c***  reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'
c
c   g_a         - handle for Global Array that is being converted SL format
c   dim1, dim2  - dimensions of Global Array
c   nbr         - number of block rows?
c   nbc         - number of block columns?
c   s_a         - local array holding SL formatted data
c   lda4        - leading dimension of s_a
      implicit none
#include "scalapack_.fh"
#include "global.fh"
      SCALAPACKINT  nbr, nbc, lda4, dim1, dim2
      integer g_a
      real*8  s_a(lda4,*)
      integer row, col, tcol, trow, rbase, cbase
      SCALAPACKINT pcol, prow
      integer lda
      integer r0i,r1i,r0im1,r1im1,rowl,coll
      integer marg1,marg2
      logical putpending

c***  Synchronize at the beginning
c
      lda=lda4
      rbase = 1
      cbase = 1
c
      Trow = nbr * nprow2
      Tcol = nbc * npcol2
      do col  = 1, dim2, nbc
         ! processor column that holds "col"
         pcol = mod(col,Tcol)/nbc
         if(mycol2.eq.pcol) then
            r0im1=-9999
            r1im1=0
            putpending=.false.
            marg1=col+nbc-1
            marg2=dim2
            coll=min(marg1,marg2)

            do row  = 1, dim1, nbr
               ! processor row that holds "row"
               prow = mod(row,Trow)/nbr
               if(myrow2.eq.prow) then
                  if(.not.putpending) then
                     r0im1=row
                     marg1=row+nbr-1
                     marg2=dim1
                     r1im1=min(marg1,marg2)
                  endif
                  r0i=row
                  marg1=row+nbr-1
                  marg2=dim1
                  r1i=min(marg1,marg2)

                  if(r0i.eq.(r1im1+1).and.r1i.ne.dim1) then
                     r1im1=r1i
                     putpending=.true.
                  else
                  rowl=r1i
                  call ga_get(g_a,r0im1,rowl,col,
     $                    coll, s_a(rbase,cbase), lda)
                     putpending=.false.
                     rbase = rbase + rowl - r0im1 +1
                     r0im1=-1
                     r1im1=-1
                  endif
               endif
            enddo
            rbase = 1
            cbase = cbase + nbc
         endif
      enddo

c**** ... and at the end
      end



      subroutine ga_from_SL2_(g_a,dim1,dim2,nbr, nbc, s_a, lda4)
c
c***  transforms a matrix from SL to GA format
c***  reference: Dongarra et al, 'A look at scalable dense lin. alg. libs'
c
      implicit none
#include "scalapack_.fh"
#include "global.fh"
      integer g_a,lda
      SCALAPACKINT  nbr, nbc, lda4, dim1, dim2
      real*8 s_a(lda4,*)
      integer  tcol, trow
      integer rbase, cbase
      SCALAPACKINT pcol, prow
      integer row,col,rowl,coll
      integer r0i,r1i,r0im1,r1im1
      integer marg1,marg2
      logical putpending
c
      lda=lda4
c
c**** Syncronize at the beginning
c
      rbase = 1
      cbase = 1
c
      Trow = nbr * nprow2
      Tcol = nbc * npcol2
c
      do col  = 1, dim2, nbc
         ! processor column that holds "col"
         pcol = mod(col,Tcol)/nbc
         if(mycol2.eq.pcol) then
            r0im1=-9999
            r1im1=0
            putpending=.false.
            marg1=col+nbc-1
            marg2=dim2
            coll=min(marg1,marg2)
            do row  = 1, dim1, nbr
               ! processor row that holds "row"
               prow = mod(row,Trow)/nbr
               if(myrow2.eq.prow) then
                  if(.not.putpending) then
                     r0im1=row
                     marg1=row+nbr-1
                     marg2=dim1
                     r1im1=min(marg1,marg2)
                  endif
                  r0i=row
                  marg1=row+nbr-1
                  marg2=dim1
                  r1i=min(marg1,marg2)

                  if(r0i.eq.(r1im1+1).and.r1i.ne.dim1) then
                     r1im1=r1i
                     putpending=.true.
                  else
                  rowl=r1i
                  call ga_put(g_a,r0im1,rowl,col,
     $                    coll, s_a(rbase,cbase), lda)
                     putpending=.false.
                     rbase = rbase + rowl - r0im1 +1
                     r0im1=-1
                     r1im1=-1
                  endif
               endif
            enddo
            rbase = 1
            cbase = cbase + nbc
         endif
      enddo

c**** ... and at the end
      end
#else
      subroutine dummy_ga_pdsyevx()
      return
      end
#endif
