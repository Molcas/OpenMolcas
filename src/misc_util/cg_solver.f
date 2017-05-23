************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2014, Ignacio Fdez. Galvan                             *
************************************************************************
*-----------------------------------------------------------------------
* <DOC>
*   <NAME>CG\_Solver</NAME>
*   <SYNTAX>Call CG\_Solver(n,nij,A,ija,b,x,info)</Syntax>
*   <ARGUMENTS>
*     \Argument{n}{Size of the system}{Integer}{in}
*     \Argument{nij}{Size of A and ija}{Integer}{in}
*     \Argument{A}{Input matrix in dense or sparse format}{Real*8 (nij)}{in}
*     \Argument{ija}{Index vector of matrix A}{Integer (*)}{in}
*     \Argument{b}{Vector of independent terms}{Real*8 (n)}{in}
*     \Argument{x}{Solution vector}{Real*8 (n)}{inout}
*     \Argument{info}{Status info}{Real*8}{inout}
*   </ARGUMENTS>
*   <PURPOSE>Solve a system of linear equations with a preconditioned conjugate gradients method</PURPOSE>
*   <DEPENDENCIES>Sp\_MV,Sp\_ICD,Sp\_Transpose,Sp\_TriSolve,dcopy,dscal,daxpy,dsymv,dpotrf,dpotrf</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Given a system of linear equations $A x = b$, this routine attempts to solve
*     it with a preconditioned conjugate gradients method.
*     The matrix A should be symmetric positive definite.
*     The method is useful with a sparse matrix, where A*x can be computed efficiently.
*     The preconditioner is an incomplete Cholesky decomposition, which becomes complete
*     with a dense matrix, and then it should converge at the first iteration.
*     If a dense matrix is used, the ija(1) value must be 0.
*     On output, info contains the number of iterations needed for convergence,
*     or -1 if it did not converge.
*   </DESCRIPTION>
* </DOC>
*-----------------------------------------------------------------------
      SUBROUTINE CG_Solver(n,nij,A,ija,b,x,info)
      IMPLICIT NONE
      INTEGER n, ija(*), nij, maxk, recomp, k, info
      REAL*8 A(nij), b(n), x(n)
      REAL*8, DIMENSION(:), ALLOCATABLE :: r, p, Ap, z, y
      REAL*8 rr, alpha, beta
      REAL*8 Thr, RelThr
      INTEGER ipLo, ipUp, ipijLo, ipijUp
      PARAMETER (Thr=1.0d-20)
      real*8 ddot_
      external ddot_
      LOGICAL Sparse
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      CALL mma_allocate(r,n)
      CALL mma_allocate(p,n)
      CALL mma_allocate(Ap,n)
      CALL mma_allocate(z,n)
      CALL mma_allocate(y,n)

c     Same algorithm, with different calls for sparse or dense matrices
      Sparse=(ija(1).GT.0)

      maxk=MAX(10,n*n)
      recomp=MAX(50,INT(n/DBLE(10)))
      call dcopy_(n,b,1,r,1)

      IF (Sparse) THEN
        CALL Allocate_Work(ipLo,nij)
        CALL Allocate_Work(ipUp,nij)
        CALL Allocate_iWork(ipijLo,nij)
        CALL Allocate_iWork(ipijUp,nij)
c
c       Compute the preconditioner
        CALL Sp_ICD(n,A,ija,Work(ipLo),iWork(ipijLo))
        CALL Sp_Transpose(n,Work(ipLo),iWork(ipijLo),
     &                          Work(ipUp),iWork(ipijUp),nij)

c
c       Initial guess: r = A*x-b
        CALL Sp_MV(n,-One,A,ija,x,One,r)
        CALL Sp_TriSolve(n,'L',Work(ipLo),iWork(ipijLo),r,y)
        CALL Sp_TriSolve(n,'U',Work(ipUp),iWork(ipijUp),y,z)
        call dcopy_(n,z,1,p,1)
        rr=DDot_(n,z,1,r,1)
        RelThr=Thr*MAX(rr,One)
        k=1
        DO WHILE ((ABS(rr).GE.RelThr).AND.(k.LE.maxk))
          CALL Sp_MV(n,One,A,ija,p,Zero,Ap)
          alpha=rr/DDot_(n,p,1,Ap,1)
          call daxpy_(n,alpha,p,1,x,1)
          beta=rr
c
c         Recompute or update the residual
          IF (MOD(k,recomp).EQ.0) THEN
            call dcopy_(n,b,1,r,1)
            CALL Sp_MV(n,-One,A,ija,x,One,r)
          ELSE
            call daxpy_(n,-alpha,Ap,1,r,1)
          END IF
          CALL Sp_TriSolve(n,'L',Work(ipLo),iWork(ipijLo),r,y)
          CALL Sp_TriSolve(n,'U',Work(ipUp),iWork(ipijUp),y,z)
          rr=DDot_(n,z,1,r,1)
          call dscal_(n,rr/beta,p,1)
          call daxpy_(n,One,z,1,p,1)
          k=k+1
        END DO
        CALL Free_Work(ipLo)
        CALL Free_Work(ipUp)
        CALL Free_iWork(ipijLo)
        CALL Free_iWork(ipijUp)
      ELSE
        CALL Allocate_Work(ipLo,nij)

c
c       With a dense matrix, the preconditioner could be replaced with
c       something else, otherwise this is just solving the system with
c       a direct method
        call dcopy_(n*n,A,1,Work(ipLo),1)
        call dpotrf_('L',n,Work(ipLo),n,info)

        CALL DSyMV('L',n,-One,A,n,x,1,One,r,1)
        call dcopy_(n,r,1,z,1)
        CALL DPoTrS('L',n,1,Work(ipLo),n,z,n,info)
        call dcopy_(n,z,1,p,1)
        rr=DDot_(n,z,1,r,1)
        RelThr=Thr*MAX(rr,One)
        k=1
        DO WHILE ((ABS(rr).GE.RelThr).AND.(k.LE.maxk))
          Call dGeMV_('N',n,n,One,A,n,p,1,Zero,Ap,1)
          alpha=rr/DDot_(n,p,1,Ap,1)
          call daxpy_(n,alpha,p,1,x,1)
          beta=rr
          IF (MOD(k,recomp).EQ.0) THEN
            call dcopy_(n,b,1,r,1)
            Call dGeMV_('N',n,n,-One,A,n,x,1,One,r,1)
          ELSE
            call daxpy_(n,-alpha,Ap,1,r,1)
          END IF
          call dcopy_(n,r,1,z,1)
          CALL DPoTrS('L',n,1,Work(ipLo),n,z,n,info)
          rr=DDot_(n,z,1,r,1)
          call dscal_(n,rr/beta,p,1)
          call daxpy_(n,One,z,1,p,1)
          k=k+1
        END DO
        CALL Free_Work(ipLo)
      END IF

c
c     Set the return value
      IF (k.LE.maxk) THEN
        info=k
      ELSE
        info=-1
      END IF

      CALL mma_deallocate(r)
      CALL mma_deallocate(p)
      CALL mma_deallocate(Ap)
      CALL mma_deallocate(z)
      CALL mma_deallocate(y)

      END
