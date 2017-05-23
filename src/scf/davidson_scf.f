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
*               2017, Roland Lindh                                     *
************************************************************************
*-----------------------------------------------------------------------
* <DOC>
*   <NAME>Davidson</NAME>
*   <SYNTAX>Call Davidson_SCF(HDiag,g,m,k,Fact,Eig,Vec,MemRsv,iRC)</SYNTAX>
*   <ARGUMENTS>
*     \Argument{HDiag}{Diagonal of the Hessian matrix}{Real*8 m}{in}
*     \Argument{g}{Gradient vector}{Real*8 m}{in}
*     \Argument{m}{Size of diagonal Hessian and gradient}{Integer}{in}
*     \Argument{k}{Number of lowest eigenvalues to compute}{Integer}{inout}
*     \Argument{Fact}{Scaling factor}{Real*8 }{in}
*     \Argument{Eig}{Lowest eigenvalues}{Real*8 (k)}{out}
*     \Argument{Vec}{Lowest eigenvectors}{Real*8 (m+1,k)}{inout}
*     \Argument{MemRsv}{Amount of reserved memory}{Integer}{in}
*     \Argument{iRC}{Return code (0 if converged)}{Integer}{out}
*   </ARGUMENTS>
*   <PURPOSE>Compute the lowest k eigenvalues of a symmetric matrix</PURPOSE>
*   <DEPENDENCIES>Add\_Vector</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY>Roland Lindh</MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Simple application of the Davidson procedure to obtain the lowest k eigenvalues
*     and corresponding eigenvectors of a symmetric matrix.
*     On input, Vec can contain an initial guess for the eigenvectors (from a previous
*     run with smaller k, for example), only the non-zero vectors are used.
*     If k $>$ n on input, it will be set to n on output.
*     This routine is adapted to an augmented Hessian, which is not explicitly expressed
*     rather the original Hessian is implicitly there, via a diagonal and an on-the-fly
*     update when multiplied by a vector, and the gradient is explicit there.
*   </DESCRIPTION>
* </DOC>
*-----------------------------------------------------------------------
*define _DEBUG_SORUPV_
#ifdef _DEBUG_SORUPV_
      SUBROUTINE Davidson_SCF(HDiag,g,m,k,Fact,Eig,Vec,MemRsv,iRC,
     &                        iter_SCF,Update_H)
#else
      SUBROUTINE Davidson_SCF(HDiag,g,m,k,Fact,Eig,Vec,MemRsv,iRC)
#endif
      IMPLICIT NONE
      INTEGER m,n,k,iRC, MemRsv
      REAL*8  HDiag(m),g(m),Eig(k),Vec(m+1,k), Fact
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: Sub, Ab
      REAL*8, DIMENSION(:), ALLOCATABLE :: Eig_old, EVec, Proj, EVal
      INTEGER, DIMENSION(:), ALLOCATABLE :: Index_D
#ifdef _DEBUG_SORUPV_
      Real*8, Allocatable:: Hessian(:,:), Vector(:)
      Logical :: Active=.False.
      Save Active
      Real*8, Allocatable:: Hss(:,:), Gamma(:), Delta(:)
      Save Hss
      Integer iter_scf, inode
      Logical Update_H
#include "file.fh"
#include "llists.fh"
#endif
      REAL*8 Aux,Thr,Thr2,Thr3,Conv,Alpha,tmp
      real*8 ddot_
      INTEGER mk,old_mk,mink,maxk,ig,info,nTmp,iter,maxiter
      INTEGER i,j,ii,jj
      INTEGER ipTmp,ipDum
      INTEGER ipDiag,ipTVec,ipTAV,ipTRes
      LOGICAL Last,Augmented,Reduced
      external ddot_
      PARAMETER (Thr=1.0D-6, maxiter=300, Thr2=1.0D-16, Thr3=1.0D-16)
*
#include "stdalloc.fh"
#include "real.fh"
#include "WrkSpc.fh"
      INTEGER iPrint,iRout
#include "print.fh"

      iRout=216
      iPrint=nPrint(iRout)
      n=m+1

*define _DEBUG_
#ifdef _DEBUG_
      Call NrmClc(HDiag,m,'Davidson_SCF','HDiag')
      Call NrmClc(    g,m,'Davidson_SCF','g')
*     CALL RecPrt('HDiag',' ',HDiag,m,1)
*     CALL RecPrt('g',' ',g,m,1)
      Do i = 1, m
         Write (6,'(2G10.1)') g(i),HDiag(i)
      End Do
#endif

*---- Initialize some parameters
*      mk   = subspace size (initially k)
*      maxk = maximum subspace size (25 if k=1)
*      mink = subspace size to reduce to when the maximum is exceeded (5 if k=1)
*
*     (do not use MIN here, so that it's possible to use an explicit value for k)
      IF (k.GT.n) k=n
      mink=MIN(MAX(k+2,5),n)
      maxk=MIN(5*mink,n)
      mk=k
      iRC=0

*---- Allocate matrices
*      Sub  = Vectors (columns) defining the subspace (maximum maxk vectors)
*      Ab   = A*b vectors (A * Sub)
*      EVal  = Set of computed eigenvalues (maximum maxk elements)
*      EVec  = Set of computed eigenvectors, in the subspace (maximum maxk*maxk)
*
      CALL mma_allocate(Eig_old,k,label="Eig_old")
      Call mma_allocate(Sub,n,maxk,Label='Sub')
      Call mma_allocate(Ab ,n,maxk,Label='Ab ')
      Call mma_allocate(Proj,maxk*maxk,Label='Proj')
      Call mma_allocate(EVal,maxk     ,Label='EVal')
      Call mma_allocate(EVec,maxk*maxk,Label='EVec')
      CALL DZero(Ab,n*maxk)
      CALL DZero(EVal,maxk)
      CALL DZero(EVec,maxk*maxk)

*---- Build an index of sorted diagonal elements in A
*
      Call mma_allocate(Index_D,n,Label='Index_D')
      DO i=1,n
         Index_D(i)=i
      END DO
*
      DO i=1,n
         ii=Index_D(i)
         If (ii.eq.n) Then
            Aux=0.0D0
         Else
            Aux=HDiag(ii)
         End If
         ii=i
         DO j=i,n
            jj=Index_D(j)
            If (jj.eq.n) Then
               Tmp=0.0D0
            Else
               Tmp=HDiag(jj)
            End If
            IF (Tmp .LT. Aux) THEN
               Aux=Tmp
               ii=j
            END IF
         END DO
         IF (ii .NE. i) THEN
            jj=Index_D(ii)
            Index_D(ii)=Index_D(i)
            Index_D(i)=jj
         END IF
      END DO
#ifdef _DEBUG_
*     Write (6,*) 'Index_D=',Index_D
#endif

*---- Setup the initial subspace
*      Read the non-linear-dependent columns from the initial eigenvector matrix
*      Fill up to mk with selected base vectors from the initial matrix
*       (those corresponding to the lowest diagonal elements)
*      The rest is set to zero, just in case
*
      nTmp=0
      CALL Allocate_Work(ipTmp,n)
      DO i=1,k
        call dcopy_(n,Vec(1,i),1,Work(ipTmp),1)
        CALL Add_Vector(n,nTmp,Sub,Work(ipTmp),Thr3)
      END DO
*
      ii=0
      CALL DZero(Work(ipTmp),n)
*
      DO WHILE ((nTmp .LT. mk) .AND. (ii .LT. n))
         ii=ii+1
         jj=Index_D(ii)
*
*        A large value indicated a forbidden rotation.
*        A negative value indicates large rotation to another
*        global minimum. Avoid these!!!
*
         If (HDiag(jj).lt.1.0D20.and.HDiag(jj).gt.-0.10D0) Then
            Work(ipTmp+jj-1)=One
            CALL Add_Vector(n,nTmp,Sub,Work(ipTmp),Thr3)
            Work(ipTmp+jj-1)=Zero
         End If
      END DO
*
*     ig will be a global counter to loop across all n base vectors
      ig=ii
      CALL Free_Work(ipTmp)
      CALL DZero(Sub(1,mk+1),(maxk-mk)*n)

*---- Iterative procedure starts here
*      mk     = subspace size at each iteration
*      old_mk = value of mk at the previous iteration
*
      Augmented=.FALSE.
      Reduced=.FALSE.
      Last=.FALSE.
      old_mk=0
      iter=0
      CALL Allocate_Work(ipDum,1)
      CALL Allocate_Work(ipDiag,n)
      CALL Allocate_Work(ipTVec,n)
      CALL Allocate_Work(ipTAV,n)
      CALL Allocate_Work(ipTRes,n)
      DO WHILE (.NOT. Last)
        iter=iter+1
        IF (iter .GT. 1) call dcopy_(k,Eig,1,Eig_old,1)
#ifdef _DEBUG_
        IF (.NOT. Reduced) THEN
          WRITE(6,'(A)') '---------------'
          WRITE(6,'(A,1X,I5)') 'Iteration',iter
        END IF
        CALL RecPrt('Orthonormalized subspace',' ',Sub,n,mk)
#endif

*----   Compute the matrix product
*        Ab = A * Sub
*        Only the new vectors since the last iterations need to be
*        calculated
*
*       Note that the A-matrix is the augmented Hessian of a rs-rfo
*       approach. The A matrix is not explicitly stored but rather only
*       the associated gradient is. The original Hessian is implicitly
*       there and a vector corresponding the contraction of the updated
*       Hessian and a trial vector can be computed on-the-fly.
*
        Do j=old_mk,mk-1
#ifdef _DEBUG_
           Write (6,*) 'Davidson_SCF: j,Fact=',j,Fact
           Call NrmClc(Sub(1,j+1),n,'Davidson_SCF','Sub(1,j+1)')
           Call RecPrt('Sub',' ',Sub(1,j+1),1,n)
#endif
*
*define _DEBUG_SAVE_ _DEBUG_
#undef _DEBUG_
#ifdef _DEBUG_SORUPV_
           If (iter.eq.1.and.j.eq.old_mk) Then
*
*             The code below is to debug the on-the-fly update plus
*             contraction-with-a-trial-vector code. Here we will
*             keep the explicit Hessian in memory and update it step
*             by step and then compare the result. Did I say that this
*             is memeory demanding. Speaking about this, the procedure
*             fail in deallocating the memory of Hss - how cares.
*
              Call mma_allocate(Vector,m,Label='Vector')
              Call mma_allocate(Hessian,m,m,Label='Hessian')
#ifdef _DEBUG_
*
*             Tweak the code to spit out the updated Hessian from
*             SorUpV
*
              Write (6,*)
              Write (6,*) '*************************************'
              Write (6,*) '*************************************'
              Do i = 1, m
                 Call FZero(Hessian(1,i),m)
                 Call FZero(Vector,m)
                 Vector(i)=1.0D0
*                Call SOrUpV(MemRsv,Vector,HDiag,m,Hessian(1,i),'GRAD')
                 Call SOrUpV2(MemRsv,Vector,HDiag,m,Hessian(1,i),'GRAD')
              End Do
              Call NrmClc(Hessian,m**2,'Davidson_SCF','Hessian')
              Call RecPrt('Hessian',' ',Hessian,m,m)
*             Call RecPrt('HDiag',' ',HDiag,1,m)
#endif
              Call mma_deallocate(Hessian)
*
              If (.NOT.Active) Then
                 Active=.True.
                 Call mma_allocate(Hss,m,m,Label='Hss')
                 Call FZero(Hss,m**2)
                 Call DCopy_(m,HDiag,1,Hss,m+1)
              Else If (Update_H) Then
                 Update_H=.False.
*
*                Do the DFP update
*
                 Call mma_allocate(Gamma,m,Label='Gamma')
                 Call mma_allocate(Delta,m,Label='Delta')
                 Call FZero(Vector,m)
*
                 Call GetNod(iter_SCF-1,LLdGrd,inode)
                 If (inode.eq.0) Then
                    Write (6,*) 'Gamma: inode.eq.0'
                    Call Abend()
                 End If
                 Call iVPtr(LudGd,Gamma,m,inode)
*
                 Call GetNod(iter_SCF-1,LLDelt,inode)
                 If (inode.eq.0) Then
                    Write (6,*) 'Gamma: inode.eq.0'
                    Call Abend()
                 End If
                 Call iVPtr(LuDel,Delta,m,inode)
*
#ifdef _DEBUG_
                 Write (6,*) 'n-1=',iter_SCF-1
                 Call NrmClc(Delta,m,'Davison_SCF','dX(n-1)')
                 Call NrmClc(Gamma,m,'Davison_SCF','dg(n-1)')
#endif
*
*                Do the actual update.
*
*                Call BFGS_HInv(Hss,m,Vector,Delta,Gamma)
                 Call DFP(Hss,m,Vector,Delta,Gamma)
*
                 Call mma_deallocate(Gamma)
                 Call mma_deallocate(Delta)
*
              End If
              Call mma_deallocate(Vector)
#ifdef _DEBUG_
              Write (6,*) '*************************************'
              Write (6,*) '*************************************'
              Write (6,*)
#endif
           End If
#ifdef _DEBUG_
           Call NrmClc(Hss    ,m**2,'Davidson_SCF','Hss    ')
           Call RecPrt('Davidsion_SCF: Hss',' ',Hss,m,m)
#endif
*
           Call DGEMM_('N','N',
     &                 m,1,m,
     &                 1.0D0,Hss,m,
     &                       Sub(1,j+1),m,
     &                 0.0D0,Ab(1,j+1),m)
#else
*
*          Pick up the contribution for the updated Hessian (BFGS update)
*
*          Call SOrUpV(MemRsv,Sub(1,j+1),HDiag,m,Ab(1,j+1),'GRAD')
           Call SOrUpV2(MemRsv,Sub(1,j+1),HDiag,m,Ab(1,j+1),'GRAD')
#endif
*define _DEBUG_ _DEBUG_SAVE_
           Call DScal_(m,1.0D0/Fact,Ab(1,j+1),1)
*
*          Add contribution from the gradient
*
           tmp=Sub(n,j+1)
           Call DaXpY_(m,1.0D0/Sqrt(Fact),g,1,Ab(1,j+1),1)
*
           Ab(n,j+1) = DDot_(m,g,1,Sub(1,j+1),1)
#ifdef _DEBUG_
           Call NrmClc(Ab(1,j+1),n,'Davidson_SCF','Ab(1,j+1)')
           Call RecPrt('Ab',' ',Ab(1,j+1),1,n)
#endif
*
        End Do
*
*----   Compute the matrix to diagonalize (symmetric)
*        Proj = Sub^t * Ab
*        Again, only the new rows/columns are needed
*
        IF (old_mk .EQ. 0) THEN
           CALL DGeMM_('T','N',
     &                 mk,mk,n,
     &                 One,Sub,n,
     &                     Ab,n,
     &                 Zero,Proj,maxk)
        ELSE
          DO i=0,mk-1
            DO j=MAX(old_mk,i),mk-1
              Proj(1+i*maxk+j)=DDot_(n,Sub(1,j+1),1,Ab(1,i+1),1)
            END DO
          END DO
        END IF

*----   Compute the eigenvalues of the projected matrix
*        Make sure the eigenpairs are sorted
*        If the subspace has been reduced, no need to compute new eigenpairs
*
        IF (.NOT. Reduced) THEN
#ifdef _DEBUG_
          WRITE(6,'(2X,A,1X,I5)') 'Solving for subspace size:',mk
#endif
          call dcopy_(maxk*maxk,Proj,1,EVec,1)
          call dsyev_('V','L',mk,EVec,maxk,EVal,
     &                          Work(ipDum),-1,info)
          nTmp=INT(Work(ipDum))
          CALL Allocate_Work(ipTmp,nTmp)
          call dsyev_('V','L',mk,EVec,maxk,EVal,
     &                          Work(ipTmp),nTmp,info)
          CALL Free_Work(ipTmp)
          CALL JacOrd2(EVal,EVec,mk,maxk)
          call dcopy_(k,EVal,1,Eig,1)
#ifdef _DEBUG_
          CALL RecPrt('Current guess',' ',Eig,1,k)
#endif
        END IF
#ifdef _DEBUG_
        IF (iPrint .GE. 99) THEN
          CALL RecPrt('Eigenvalues',' ',EVal,1,mk)
          CALL SubRecPrt('Subspace Eigenvectors',' ',EVec,
     &      maxk,mk,mk)
          WRITE(6,*)
        END IF
#endif
*                                                                      *
************************************************************************
************************************************************************
*----   Check for convergence                                          *
*
*       Converge if the change in the eigenvalues is small
*       (but if a mink size has been reached)
*       Converge if the full system has been solved
*       Stop if the number of iterations exceeds the maximum
*       Stop if no new vectors to add are found
*
*                                                                      *
************************************************************************
*                                                                      *
*       Compute Conv
*
        IF (iter .GT. 1) THEN
*                                                                      *
************************************************************************
*                                                                      *
           Conv=Zero
           DO i=1,k
              IF (ABS(Eig(i)) .GT. Thr2) THEN
                 Conv=MAX(Conv,ABS((Eig(i)-Eig_old(i))/Eig(i)))
              ELSE
                 Conv=MAX(Conv,ABS(Eig(i)-Eig_old(i)))
              END IF
           END DO
#ifdef _DEBUG_
           IF (Augmented)
     &       WRITE(6,'(2X,A,1X,G12.6)') 'Maximum relative eigenvalue '//
     &                                 'change:',Conv
#endif
*                                                                      *
************************************************************************
*                                                                      *
        ELSE
*                                                                      *
************************************************************************
*                                                                      *
           Conv=Ten*Thr
*                                                                      *
************************************************************************
*                                                                      *
        END IF
*                                                                      *
************************************************************************
*                                                                      *
*       Now check for convergence
*
        old_mk=mk
*                                                                      *
************************************************************************
*                                                                      *
        IF (Augmented .AND. (Conv .LE. Thr)
     &      .AND. (mk .GE. mink)) THEN
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
              WRITE(6,'(A)') 'Converged due to small change'
#endif
          Last=.TRUE.
*                                                                      *
************************************************************************
*                                                                      *
        ELSE IF (mk .EQ. n) THEN
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
          WRITE(6,'(A)') 'Complete system solved'
#endif
          Last=.TRUE.
*                                                                      *
************************************************************************
*                                                                      *
        ELSE IF (iter .GE. maxiter) THEN
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
          WRITE(6,'(A)') 'Not converged'
#endif
          Last=.TRUE.
          iRC=1

*----     Reduce the subspace size if it exceeds the maximum (maxk)
*          Sub' = Sub * Vec(1:mink)
*          Sub' should be orthonormal if Sub is orthonormal
*        (A reduction does not consume an iteration)
*       There is also a reduction if the process is stagnated
*
*                                                                      *
************************************************************************
*                                                                      *
        ELSE IF ((MIN(mk+k,n) .GT. maxk) .OR. (iRC .EQ. 2)) THEN
*                                                                      *
************************************************************************
*                                                                      *
          IF (iRC .EQ. 2) iRC=0
#ifdef _DEBUG_
          WRITE(6,'(2X,A,1X,I5)') 'Reducing search space to',mink
#endif
          CALL Allocate_Work(ipTmp,mink*n)
          CALL DGeMM_('N','N',
     &                n,mink,mk,
     &                One,Sub,n,
     &                    EVec,maxk,
     &                Zero,Work(ipTmp),n)
          call dcopy_(mink*n,Work(ipTmp),1,Sub,1)
          CALL Free_Work(ipTmp)

*----     To make sure Sub' is orthonormal, add the vectors one by one
*
          j=0
          i=0
          DO WHILE ((j .LT. mink) .AND. (i .LT. mk))
            i=i+1
            Call Add_Vector(n,j,Sub,Sub(1,i),Thr3)
          END DO

*----     j should be mink, but who knows...
*
#ifdef _DEBUG_
          IF (j .LT. mink) THEN
            WRITE(6,'(2X,A,1X,I5)') 'Fewer vectors found:',j
          END IF
#endif
          CALL DZero(Sub(1,j+1),(maxk-j)*n)
          CALL DZero(Ab(1,j+1),(maxk-j)*n)
          CALL DZero(EVec,maxk*maxk)
          DO i=0,j-1
            EVec(1+i*(maxk+1))=One
          END DO
          mk=j
          old_mk=0
          Augmented=.FALSE.
          Reduced=.TRUE.
          iter=iter-1
*                                                                      *
************************************************************************
*                                                                      *
        ELSE
*                                                                      *
************************************************************************
*                                                                      *
*----     Expand the subspace
*          For each eigenpair i of the first k,
*          check convergence for the residuals r:
*           r = Ab * Vec(i) - Val(i) * Sub * Vec(i)
*          Add a new vector, orthonormalized with the previous vectors,
*          computed from r and the eigenpair
*          (different possible variants)
*
          CALL Allocate_Work(ipTmp,n)
          Conv=Zero
*
          jj=0
          DO i=0,k-1
*            Vector in full space: Sub*Vec(i)
             Call dGeMV_('N',n,mk,One,Sub,n,
     &                               EVec(1+i*maxk),1,
     &                           Zero,Work(ipTVec),1)
*            Product of matrix and vector: Ab*Vec(i)
             Call dGeMV_('N',n,mk,One,Ab,n,
     &                               EVec(1+i*maxk),1,
     &                           Zero,Work(ipTAV),1)
*            Residual: (A-Val(i))*Vec(i) = Ab*Vec(i) - Val(i)*Sub*Vec(i)
             call dcopy_(n,Work(ipTAV),1,Work(ipTRes),1)
             call daxpy_(n,-EVal(1+i),Work(ipTVec),1,Work(ipTRes),1)
             Conv=MAX(Conv,DDot_(n,Work(ipTRes),1,Work(ipTRes),1))

*----        Scale vector, orthonormalize, and add to subspace
*
*            Diagonal matrix to scale the vectors: 1/(A(j,j)-Val(i))
             DO j=0,n-1
*               Aux=A((j+1)*(j+2)/2)-EVal(1+i)
                If (j.eq.n-1) Then
                   Aux=          -Eval(1+i)
                Else
                   Aux=HDiag(j+1)-Eval(1+i)
                End If
                If (HDiag(j+1).lt.1.0D20) Then
                   Work(ipDiag+j)=One/SIGN(MAX(ABS(Aux),Thr2),Aux)
                Else
                   Work(ipDiag+j)=1.0D20
                End If
             END DO
*
*            scale
             DO j=0,n-1
                If (Work(ipDiag+j).lt.1.0D02) Then
                   Work(ipTmp+j)=Work(ipTRes+j)*Work(ipDiag+j)
                Else
                   Work(ipTmp+j)=0.0D0
                End If
             END DO
*
             Alpha=Zero
             DO j=0,n-1
                If (Work(ipDiag+j).lt.1.0D02) Then
                   Alpha=Alpha+Work(ipDiag+j)*Work(ipTVec+j)**2
                End If
             END DO
             Alpha=DDot_(n,Work(ipTVec),1,Work(ipTmp),1)/Alpha
*            subtract
             DO j=0,n-1
                If (Work(ipDiag+j).lt.1.0D02) Then
                   Work(ipTVec+j)=Work(ipTVec+j)*Work(ipDiag+j)
                Else
                   Work(ipTVec+j)=0.0D0
                End If
             END DO
             call daxpy_(n,-Alpha,Work(ipTVec),1,Work(ipTmp),1)
*
             IF (mk+jj .LE. n-1) THEN
                jj=mk+jj
                CALL Add_Vector(n,jj,Sub,Work(ipTmp),Thr3)
                jj=jj-mk
             END IF
          END DO
*
#ifdef _DEBUG_
          WRITE(6,'(2X,A,1X,G12.6)') 'Maximum residual:',Conv
#endif
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
          IF ((Conv .LT. Thr3) .AND. (mk .GE. mink)) THEN
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
#ifdef _DEBUG_
             WRITE(6,'(A)') 'Converged due to small residual'
#endif
             Last=.TRUE.
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
          ELSE
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
             mk=MIN(mk+jj,n)

*----        If no new vector is found to add to the subspace, we are in
*            trouble. Try to find a non-linear-dependent base vector in
*            the original matrix
*
             IF (jj .EQ. 0) THEN
#ifdef _DEBUG_
               WRITE(6,'(A)') 'Process stagnated'
#endif
               IF (mk .LT. maxk) THEN
                  CALL DZero(Work(ipTmp),n)
                  i=0
*
                  DO WHILE ((jj .LT. 1) .AND. (i .LT. n))
                     i=i+1
                     ig=MOD(ig,n)+1
                     ii=Index_D(ig)
*
*                    Avoid excplicitly rotations between fermions of
*                    different types. Avoid rotations which will be
*                    large,
*
                     If (HDiag(ii).lt.1.0D20  .and.
     &                   HDiag(ii).gt.-0.10D0) Then
                        Work(ipTmp+ii-1)=One
                        jj=mk+jj
                        CALL Add_Vector(n,jj,Sub,Work(ipTmp),Thr3)
                        Work(ipTmp+ii-1)=Zero
                        jj=jj-mk
                     End If
*
                  END DO
*
                  mk=MIN(mk+jj,n)
                  IF (jj .GT. 0) Augmented=.TRUE.
               END IF
               IF (jj .EQ. 0) THEN
                  Augmented=.FALSE.
                  iRC=2
               ENDIF
            ELSE
               Augmented=.TRUE.
            END IF
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
          END IF
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
          CALL Free_Work(ipTmp)
          Reduced=.FALSE.
*                                                                      *
************************************************************************
*                                                                      *
        END IF
*                                                                      *
************************************************************************
*                                                                      *
      END DO
*
      CALL Free_Work(ipDum)
      CALL Free_Work(ipDiag)
      CALL Free_Work(ipTVec)
      CALL Free_Work(ipTAV)
      CALL Free_Work(ipTRes)
      Call mma_deallocate(Index_D)

*---- Store the current lowest k eigenvectors (in the full space)
*      Vec' = Sub * Vec(1:k)
*
      CALL DGeMM_('N','N',
     &            n,k,mk,
     &            One,Sub,n,
     &                EVec,maxk,
     &            Zero,Vec,n)

#ifdef _DEBUG_
      Call NrmClc(Vec(1,1),m,'Davidson_SCF','Vec(1-m)')
      Call NrmClc(Vec(n,1),1,'Davidson_SCF','Vec(n)')
#endif
      Call mma_deallocate(Sub)
      Call mma_deallocate(Ab)
      Call mma_deallocate(Proj)
      Call mma_deallocate(EVal)
      Call mma_deallocate(EVec)
      CALL mma_deallocate(Eig_old)

      END
#ifdef _DEBUG_SORUPV_
      Subroutine BFGS_HInv(Hss,m,HdX,dX,dg)
      Implicit Real*8 (a-h,o-z)
      Real*8 Hss(m,m), HdX(m), dX(m), dg(m)
*
#ifdef _DEBUG_
      Call RecPrt('BFGS: H(n-1)',' ',Hss,m,m)
      Call RecPrt('BFGS: dX(n-1)',' ',dX,1,m)
      Call RecPrt('BFGS: dg(n-1)',' ',dg,1,m)
#endif
      Call FZero(HdX,m)
      Call DGEMM_('N','N',
     &            m,1,m,
     &            1.0D0,Hss,m,
     &                  dX,m,
     &            0.0D0,HdX,m)
      alpha=1.0D0/DDot_(m,dX,1,dg,1)
      beta =DDot_(m,HdX,1,dX,1)
      epsilon=(1.0D0 + alpha*beta)*alpha
*define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('BFGS: HdX(n-1)',' ',HdX,1,m)
      Write (6,*) 'BFGS: a,b,e=',alpha,beta,epsilon
#endif
#undef _DEBUG_
*
      Do i = 1, m
         Do j = 1, m
            Hss(i,j) = Hss(i,j)
     &               + epsilon * dg(i) * dg(j)
     &               - alpha   * (dg(i)*HdX(j)
     &                           +HdX(i)*dg(j))
         End Do
      End Do
#ifdef _DEBUG_
      Call RecPrt('BFGS: H(n)',' ',Hss,m,m)
#endif
*
      Return
      End
#endif
