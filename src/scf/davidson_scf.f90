!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2014, Ignacio Fdez. Galvan                             *
!               2017, Roland Lindh                                     *
!***********************************************************************
!  Davidson_SCF
!
!> @brief
!>   Compute the lowest \p k eigenvalues of a symmetric matrix.
!> @author Ignacio Fdez. Galv&aacute;n
!> @modified_by Roland Lindh
!>
!> @details
!> Simple application of the Davidson procedure to obtain the lowest \p k eigenvalues
!> and corresponding eigenvectors of a symmetric matrix.
!> On input, \p Vec can contain an initial guess for the eigenvectors (from a previous
!> run with smaller \p k, for example), only the non-zero vectors are used.
!> This routine is adapted to an augmented Hessian, which is not explicitly expressed,
!> rather the original Hessian is implicitly there, via a diagonal and an on-the-fly
!> update when multiplied by a vector, and the gradient is explicit there.
!>
!> @param[in]     HDiag  Diagonal of the Hessian matrix
!> @param[in]     g      Gradient vector
!> @param[in]     m      Size of diagonal Hessian and gradient
!> @param[in]     k      Number of lowest eigenvalues to compute
!> @param[in]     Fact   Scaling factor
!> @param[out]    Eig    Lowest eigenvalues
!> @param[in,out] Vec    Lowest eigenvectors
!> @param[out]    iRC    Return code (0 if converged)
!***********************************************************************
!#define _DEBUGPRINT_
      SUBROUTINE Davidson_SCF(g,m,k,Fact,Eig,Vec,iRC)
      Use SCF_Arrays, only: HDiag
      use Constants, only: Zero, One, Ten
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT NONE
      INTEGER m,n,k,iRC
      REAL*8  g(m),Eig(k),Vec(m+1,k), Fact
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: Sub, Ab
      REAL*8, DIMENSION(:), ALLOCATABLE :: Eig_old, EVec, Proj, EVal
      INTEGER, DIMENSION(:), ALLOCATABLE :: Index_D
      REAL*8 Aux,Thr,Thr2,Thr3,Conv,Alpha,tmp
      real*8 ddot_
      INTEGER mk,old_mk,mink,maxk,ig,info,nTmp,iter,maxiter
      INTEGER i,j,ii,jj
      LOGICAL Last,Augmented,Reduced
      external ddot_
      PARAMETER (Thr=1.0D-6, maxiter=300, Thr2=1.0D-12, Thr3=1.0D-16)
      Real*8, Allocatable :: TmpVec(:), Diag(:), TVec(:), TAV(:),TRes(:)
      Real*8 :: Dum(1)=Zero
!
#include "print.fh"
#ifdef _DEBUGPRINT_
      INTEGER iPrint,iRout

      iRout=216
      iPrint=nPrint(iRout)
#endif
      n=m+1
!#define _DEBUGCode_
#ifdef _DEBUGCode_
           Block
             Real*8, Allocatable :: Vec(:), HM(:,:), HAug(:,:)
             Real*8, Allocatable :: EVal(:), EVec(:)
             Integer ij

             Call mma_allocate(Vec,m,Label='Vec')
             Call mma_allocate(HM,m,m,Label='HM')
             HM(:,:)=Zero
             Call mma_allocate(HAug,n,n,Label='HAug')
             HAug(:,:)=Zero

             Do i = 1, m
                Vec(:)=Zero
                Vec(i)=One
                Call SOrUpV(Vec(:),m,HM(:,i),'GRAD','BFGS')
             End Do
!            Call RecPrt('HM',' ',HM,m,m)

             Call mma_allocate(EVal,m*(m+1)/2,Label='EVal')
             Call mma_allocate(EVec,m*m,Label='EVec')
             Do i = 1, m
                Do j = 1, i
                   ij=i*(i-1)/2 + j
                   EVal(ij)=HM(i,j)
                End Do
             End Do
!
!---- Set up a unit matrix
!
             call dcopy_(m*m,[Zero],0,EVec,1)
             call dcopy_(m,[One],0,EVec,m+1)
!
!----        Compute eigenvalues and eigenvectors
!
             Call NIDiag_new(EVal,EVec,m,m)
             Call Jacord(EVal,EVec,m,m)

!            Do i = 1, m
!               ij=i*(i+1)/2
!               Write (6,*) 'Eval0(ij)=',EVal(ij)
!            End Do

             Call mma_deallocate(EVal)
             Call mma_deallocate(EVec)

             Do i = 1, m
                HAug(n,i)=g(i)
                HAug(i,n)=g(i)
                Do j = 1, m
                   HAug(i,j)=HM(i,j)
                End Do
             End Do

             Call mma_allocate(EVal,n*(n+1)/2,Label='EVal')
             Call mma_allocate(EVec,n*n,Label='EVec')
!
             Do i = 1, n
                Do j = 1, i
                   ij=i*(i-1)/2 + j
                   EVal(ij)=HAug(i,j)
                End Do
             End Do
!
!---- Set up a unit matrix
!
             call dcopy_(n*n,[Zero],0,EVec,1)
             call dcopy_(n,[One],0,EVec,n+1)
!
!----        Compute eigenvalues and eigenvectors
!
             Call NIDiag_new(EVal,EVec,n,n)
             Call Jacord(EVal,EVec,n,n)

             Do i = 1, n
                ij=i*(i+1)/2
                If (EVal(ij)<Zero) Write (6,*) 'Eval(ij)=',EVal(ij)
             End Do

             Call mma_deallocate(EVal)
             Call mma_deallocate(EVec)
             Call mma_deallocate(HAug)
             Call mma_deallocate(HM)
             Call mma_deallocate(Vec)

           End Block
#endif

#ifdef _DEBUGPRINT_
      Call NrmClc(HDiag,m,'Davidson_SCF','HDiag')
      Call NrmClc(    g,m,'Davidson_SCF','g')
!     CALL RecPrt('HDiag',' ',HDiag,1,m)
!     CALL RecPrt('g',' ',g,1,m)
#endif

!---- Initialize some parameters
!      mk   = subspace size (initially k)
!      maxk = maximum subspace size (25 if k=1)
!      mink = subspace size to reduce to when the maximum is exceeded (5 if k=1)
!
      IF (k.GT.n) CALL SysAbendMsg('Davidson_SCF','Wrong k value.','')

      mink=MIN(MAX(k+2,5),n)
      maxk=MIN(5*mink,n)
      mk=k
      iRC=0

!---- Allocate matrices
!      Sub  = Vectors (columns) defining the subspace (maximum maxk vectors)
!      Ab   = A*b vectors (A * Sub)
!      EVal  = Set of computed eigenvalues (maximum maxk elements)
!      EVec  = Set of computed eigenvectors, in the subspace (maximum maxk*maxk)
!
      CALL mma_allocate(Eig_old,k,label="Eig_old")
      Call mma_allocate(Sub,n,maxk,Label='Sub')
      Call mma_allocate(Ab ,n,maxk,Label='Ab ')
      Call mma_allocate(Proj,maxk*maxk,Label='Proj')
      Call mma_allocate(EVal,maxk     ,Label='EVal')
      Call mma_allocate(EVec,maxk*maxk,Label='EVec')
      CALL FZero(Ab,n*maxk)
      CALL FZero(EVal,maxk)
      CALL FZero(EVec,maxk*maxk)

!---- Build an index of sorted diagonal elements in A
!
      Call mma_allocate(Index_D,n,Label='Index_D')
      DO i=1,n
         Index_D(i)=i
      END DO
!
      DO i=1,n
         ii=Index_D(i)
         If (ii.eq.n) Then
            Aux=Zero
         Else
            Aux=HDiag(ii)
         End If
         ii=i
         DO j=i,n
            jj=Index_D(j)
            If (jj.eq.n) Then
               Tmp=Zero
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
#ifdef _DEBUGPRINT_
      Write (6,*) 'Index_D=',Index_D
#endif

!---- Setup the initial subspace
!      Read the non-linear-dependent columns from the initial eigenvector matrix
!      Fill up to mk with selected base vectors from the initial matrix
!       (those corresponding to the lowest diagonal elements)
!      The rest is set to zero, just in case
!
      nTmp=0
      Call mma_allocate(TmpVec,n,Label='TmpVec')
      DO i=1,k
        call dcopy_(n,Vec(1,i),1,TmpVec,1)
        CALL Add_Vector(n,nTmp,Sub,TmpVec,Thr3)
      END DO
!
      ii=0
      TmpVec(:)=Zero
!
      DO WHILE ((nTmp .LT. mk) .AND. (ii .LT. n))
         ii=ii+1
         jj=Index_D(ii)
!
!        A large value indicates a forbidden rotation.
!        A negative value indicates large rotation to another
!        global minimum. Avoid these!!!
!
         If (jj.eq.n) Then
           Aux=Zero
         Else
           Aux=HDiag(jj)
         End If
         If (Aux.lt.1.0D10.and.Aux.gt.-0.10D0) Then
            TmpVec(jj)=One
            CALL Add_Vector(n,nTmp,Sub,TmpVec,Thr3)
            TmpVec(jj)=Zero
         End If
      END DO
!
!     ig will be a global counter to loop across all n base vectors
      ig=ii
      Call mma_deallocate(TmpVec)
      CALL FZero(Sub(1,mk+1),(maxk-mk)*n)

!---- Iterative procedure starts here
!      mk     = subspace size at each iteration
!      old_mk = value of mk at the previous iteration
!
      Augmented=.FALSE.
      Reduced=.FALSE.
      Last=.FALSE.
      old_mk=0
      iter=0
      Call mma_allocate(Diag,n,Label='Diag')
      Call mma_allocate(TVec,n,Label='TVec')
      Call mma_allocate(TAV ,n,Label='TAV ')
      Call mma_allocate(TRes,n,Label='TRes')
      DO WHILE (.NOT. Last)
        iter=iter+1
        IF (iter .GT. 1) call dcopy_(k,Eig,1,Eig_old,1)
#ifdef _DEBUGPRINT_
        IF (.NOT. Reduced) THEN
          WRITE(6,'(A)') '---------------'
          WRITE(6,'(A,1X,I5)') 'Iteration',iter
        END IF
!       CALL RecPrt('Orthonormalized subspace',' ',Sub,n,mk)
#endif

!----   Compute the matrix product
!        Ab = A * Sub
!        Only the new vectors since the last iterations need to be
!        calculated
!
!       Note that the A-matrix is the augmented Hessian of a rs-rfo
!       approach. The A matrix is not explicitly stored but rather only
!       the associated gradient is. The original Hessian is implicitly
!       there and a vector corresponding to the contraction of the updated
!       Hessian and a trial vector can be computed on-the-fly.
!
        Do j=old_mk,mk-1
#ifdef _DEBUGPRINT_
           Write (6,*) 'Davidson_SCF: j,Fact=',j,Fact
           Call NrmClc(Sub(1,j+1),n,'Davidson_SCF','Sub(1,j+1)')
!          Call RecPrt('Sub',' ',Sub(1,j+1),1,n)
#endif
!
!          Pick up the contribution for the updated Hessian (BFGS update)
!
           Call SOrUpV(Sub(1,j+1),m,Ab(1,j+1),'GRAD','BFGS')
!          Call RecPrt('Ab(0)',' ',Ab(1,j+1),1,n)
           Call DScal_(m,One/Fact,Ab(1,j+1),1)
!
!          Add contribution from the gradient
!
           tmp=Sub(n,j+1)
           Call DaXpY_(m,tmp/Sqrt(Fact),g,1,Ab(1,j+1),1)
!
           Ab(n,j+1) =  DDot_(m,g,1,Sub(1,j+1),1)/Sqrt(Fact)
#ifdef _DEBUGPRINT_
           Call NrmClc(Ab(1,j+1),n,'Davidson_SCF','Ab(1,j+1)')
!          Call RecPrt('Ab',' ',Ab(1,j+1),1,n)
#endif
!
        End Do
!
!----   Compute the matrix to diagonalize (symmetric)
!        Proj = Sub^t * Ab
!        Again, only the new rows/columns are needed
!
        IF (old_mk .EQ. 0) THEN
           CALL DGeMM_('T','N',mk,mk,n,One,Sub,n,Ab,n,Zero,Proj,maxk)
        ELSE
          DO i=0,mk-1
            DO j=MAX(old_mk,i),mk-1
              Proj(1+i*maxk+j)=DDot_(n,Sub(1,j+1),1,Ab(1,i+1),1)
            END DO
          END DO
        END IF

!----   Compute the eigenvalues of the projected matrix
!        Make sure the eigenpairs are sorted
!        If the subspace has been reduced, no need to compute new eigenpairs
!
        IF (.NOT. Reduced) THEN
#ifdef _DEBUGPRINT_
          WRITE(6,'(2X,A,1X,I5)') 'Solving for subspace size:',mk
#endif
          call dcopy_(maxk*maxk,Proj,1,EVec,1)
          call dsyev_('V','L',mk,EVec,maxk,EVal,Dum,-1,info)
          If (info/=0) Then
             Write (6,*) 'info(2)/=0', info
             Call Abend()
          End If
          nTmp=Max(1,INT(Dum(1)))
          Call mma_allocate(TmpVec,nTmp,Label='TmpVec')
          call dsyev_('V','L',mk,EVec,maxk,EVal,TmpVec,nTmp,info)
          If (info/=0) Then
             Write (6,*) 'info(2)/=0', info
             Call Abend()
          End If
          Call mma_deallocate(TmpVec)
          CALL SortEig(EVal,EVec,mk,maxk,1,.false.)
          call dcopy_(k,EVal,1,Eig,1)
#ifdef _DEBUGPRINT_
!         CALL RecPrt('Current guess',' ',Eig,1,k)
#endif
        END IF
#ifdef _DEBUGPRINT_
        IF (iPrint .GE. 99) THEN
!         CALL RecPrt('Eigenvalues',' ',EVal,1,mk)
!         CALL SubRecPrt('Subspace Eigenvectors',' ',EVec,maxk,mk,mk)
          WRITE(6,*)
        END IF
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!----   Check for convergence                                          *
!
!       Converge if the change in the eigenvalues is small
!       (but if a mink size has been reached)
!       Converge if the full system has been solved
!       Stop if the number of iterations exceeds the maximum
!       Stop if no new vectors to add are found
!
!                                                                      *
!***********************************************************************
!                                                                      *
!       Compute Conv
!
        IF (iter .GT. 1) THEN
!                                                                      *
!***********************************************************************
!                                                                      *
           Conv=Zero
           DO i=1,k
              IF (ABS(Eig(i)) .GT. Thr2) THEN
                 Conv=MAX(Conv,ABS((Eig(i)-Eig_old(i))/Eig(i)))
              ELSE
                 Conv=MAX(Conv,ABS(Eig(i)-Eig_old(i)))
              END IF
           END DO
#ifdef _DEBUGPRINT_
           IF (Augmented)WRITE(6,'(2X,A,1X,G12.6)') 'Maximum relative eigenvalue change:',Conv
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
        ELSE
!                                                                      *
!***********************************************************************
!                                                                      *
           Conv=Ten*Thr
!                                                                      *
!***********************************************************************
!                                                                      *
        END IF
!                                                                      *
!***********************************************************************
!                                                                      *
!       Now check for convergence
!
        old_mk=mk
!                                                                      *
!***********************************************************************
!                                                                      *
        IF (Augmented .AND. (Conv .LE. Thr) .AND. (mk .GE. mink)) THEN
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
              WRITE(6,'(A)') 'Converged due to small change'
#endif
          Last=.TRUE.
!                                                                      *
!***********************************************************************
!                                                                      *
        ELSE IF (mk .EQ. n) THEN
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
          WRITE(6,'(A)') 'Complete system solved'
#endif
          Last=.TRUE.
!                                                                      *
!***********************************************************************
!                                                                      *
        ELSE IF (iter .GE. maxiter) THEN
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
          WRITE(6,'(A)') 'Not converged'
#endif
          Last=.TRUE.
          iRC=1

!----     Reduce the subspace size if it exceeds the maximum (maxk)
!          Sub' = Sub * Vec(1:mink)
!          Sub' should be orthonormal if Sub is orthonormal
!        (A reduction does not consume an iteration)
!       There is also a reduction if the process is stagnated
!
!                                                                      *
!***********************************************************************
!                                                                      *
        ELSE IF (mk .GT. mink .AND. (MIN(mk+k,n) .GT. maxk) .OR. (iRC .EQ. 2)) THEN
!                                                                      *
!***********************************************************************
!                                                                      *
          IF (iRC .EQ. 2) iRC=0
#ifdef _DEBUGPRINT_
          WRITE(6,'(2X,A,1X,I5)') 'Reducing search space to',mink
#endif
          Call mma_allocate(TmpVec,mink*n,Label='TmpVec')
          CALL DGeMM_('N','N',n,mink,mk,One,Sub,n,EVec,maxk,Zero,TmpVec,n)
          call dcopy_(mink*n,TmpVec,1,Sub,1)
          Call mma_deallocate(TmpVec)

!----     To make sure Sub' is orthonormal, add the vectors one by one
!
          j=0
          i=0
          DO WHILE ((j .LT. mink) .AND. (i .LT. mk))
            i=i+1
            Call Add_Vector(n,j,Sub,Sub(1,i),Thr3)
          END DO

!----     j should be mink, but who knows...
!
#ifdef _DEBUGPRINT_
          IF (j .LT. mink) THEN
            WRITE(6,'(2X,A,1X,I5)') 'Fewer vectors found:',j
          END IF
#endif
          CALL FZero(Sub(1,j+1),(maxk-j)*n)
          CALL FZero(Ab(1,j+1),(maxk-j)*n)
          CALL FZero(EVec,maxk*maxk)
          DO i=0,j-1
            EVec(1+i*(maxk+1))=One
          END DO
          mk=j
          old_mk=0
          Augmented=.FALSE.
          Reduced=.TRUE.
          iter=iter-1
!                                                                      *
!***********************************************************************
!                                                                      *
        ELSE
!                                                                      *
!***********************************************************************
!                                                                      *
!----     Expand the subspace
!          For each eigenpair i of the first k,
!          check convergence for the residuals r:
!           r = Ab * Vec(i) - Val(i) * Sub * Vec(i)
!          Add a new vector, orthonormalized with the previous vectors,
!          computed from r and the eigenpair
!          (different possible variants)
!
          Call mma_allocate(TmpVec,n,Label='TmpVec')
          Conv=Zero
!
          jj=0
          DO i=0,k-1
!            Vector in full space: Sub*Vec(i)
             Call dGeMV_('N',n,mk,One,Sub,n,EVec(1+i*maxk),1,Zero,TVec,1)
!            Product of matrix and vector: Ab*Vec(i)
             Call dGeMV_('N',n,mk,One,Ab,n,EVec(1+i*maxk),1,Zero,TAV,1)
!            Residual: (A-Val(i))*Vec(i) = Ab*Vec(i) - Val(i)*Sub*Vec(i)
             call dcopy_(n,TAV,1,TRes,1)
             call daxpy_(n,-EVal(1+i),TVec,1,TRes,1)
             Conv=MAX(Conv,DDot_(n,TRes,1,TRes,1))

!----        Scale vector, orthonormalize, and add to subspace
!
!            Diagonal matrix to scale the vectors: 1/(A(j,j)-Val(i))
             DO j=0,n-1
!               Aux=A((j+1)*(j+2)/2)-EVal(1+i)
                If (j.eq.n-1) Then
                   Aux=          -Eval(1+i)
                Else
                   Aux=HDiag(j+1)-Eval(1+i)
                End If
                If (j.eq.n-1) Then
                   Diag(1+j)=One/SIGN(MAX(ABS(Aux),Thr2),Aux)
                Else
                   If (HDiag(j+1).lt.1.0D20) Then
                      Diag(1+j)=One/SIGN(MAX(ABS(Aux),Thr2),Aux)
                   Else
                      Diag(1+j)=1.0D20
                   End If
                End If
             END DO
!
!            scale
             DO j=0,n-1
                If (Diag(1+j).lt.Ten**2) Then
                   TmpVec(1+j)=TRes(1+j)*Diag(1+j)
                Else
                   TmpVec(1+j)=Zero
                End If
             END DO
!
             Alpha=Zero
             DO j=0,n-1
                If (Diag(1+j).lt.Ten**2) Then
                   Alpha=Alpha+Diag(1+j)*TVec(1+j)**2
                End If
             END DO
             Alpha=DDot_(n,TVec,1,TmpVec,1)/Alpha
!            subtract
             DO j=0,n-1
                If (Diag(1+j).lt.Ten**2) Then
                   TVec(1+j)=TVec(1+j)*Diag(1+j)
                Else
                   TVec(1+j)=Zero
                End If
             END DO
             call daxpy_(n,-Alpha,TVec,1,TmpVec,1)
!
             IF (mk+jj .LE. n-1) THEN
                jj=mk+jj
                CALL Add_Vector(n,jj,Sub,TmpVec,Thr3)
                jj=jj-mk
             END IF
          END DO
!
#ifdef _DEBUGPRINT_
          WRITE(6,'(2X,A,1X,G12.6)') 'Maximum residual:',Conv
#endif
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
          IF ((Conv .LT. Thr3) .AND. (mk .GE. mink)) THEN
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
#ifdef _DEBUGPRINT_
             WRITE(6,'(A)') 'Converged due to small residual'
#endif
             Last=.TRUE.
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
          ELSE
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
             mk=MIN(mk+jj,n)

!----        If no new vector is found to add to the subspace, we are in
!            trouble. Try to find a non-linear-dependent base vector in
!            the original matrix
!
             IF (jj .EQ. 0) THEN
#ifdef _DEBUGPRINT_
               WRITE(6,'(A)') 'Process stagnated'
#endif
               IF (mk .LT. maxk) THEN
                  TmpVec(:n)=Zero
                  i=0
!
                  DO WHILE ((jj .LT. 1) .AND. (i .LT. n))
                     i=i+1
                     ig=MOD(ig,n)+1
                     ii=Index_D(ig)
!
!                    Avoid explicitly rotations between fermions of
!                    different types. Avoid rotations which will be
!                    large,
!
                     If (ii.eq.n) Then
                       Aux=Zero
                     Else
                       Aux=HDiag(ii)
                     End If
                     If (Aux.lt.1.0D20  .and. Aux.gt.-0.10D0) Then
                        TmpVec(ii)=One
                        jj=mk+jj
                        CALL Add_Vector(n,jj,Sub,TmpVec,Thr3)
                        TmpVec(ii)=Zero
                        jj=jj-mk
                     End If
!
                  END DO
!
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
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
          END IF
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
          Call mma_deallocate(TmpVec)
          Reduced=.FALSE.
!                                                                      *
!***********************************************************************
!                                                                      *
        END IF
!                                                                      *
!***********************************************************************
!                                                                      *
      END DO
!
      Call mma_deallocate(Diag)
      Call mma_deallocate(TVec)
      Call mma_deallocate(TAV )
      Call mma_deallocate(TRes)
      Call mma_deallocate(Index_D)

!---- Store the current lowest k eigenvectors (in the full space)
!     Vec' = Sub * Vec(1:k)

      CALL DGeMM_('N','N',n,k,mk,One,Sub,n,EVec,maxk,Zero,Vec,n)

#ifdef _DEBUGPRINT_
      Call NrmClc(Vec(1,1),m,'Davidson_SCF','Vec(1-m)')
      Call NrmClc(Vec(n,1),1,'Davidson_SCF','Vec(n)')
#endif
      Call mma_deallocate(Sub)
      Call mma_deallocate(Ab)
      Call mma_deallocate(Proj)
      Call mma_deallocate(EVal)
      Call mma_deallocate(EVec)
      CALL mma_deallocate(Eig_old)

      END SUBROUTINE Davidson_SCF
