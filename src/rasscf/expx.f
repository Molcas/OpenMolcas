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
      SUBROUTINE expx(x,x2,y,v,a,b,thmax,n)
c
c     RASSCF Program version IBM-3090: SX section
c
c     Purpose: This subroutine computes exp(X), where X is formed
c              from the super-CI coefficient matrix C for one symmetry
c              block.
c              The resulting matrix exp(X) is stored in X,
c              the input matrix X is thus destroyed.
c              exp(X) = 1 + V*a*V(+) + V*b*V(+)*X   where
c              V is the matrix of eigenvectors of the Hermitian matrix X**2
c              a is a diagonal matrix: a(ii) = cos(theta(i)) - 1
c              b is a diagonal matrix: b(ii) = sin(theta(i))/theta(i)
c              with theta(i)**2 = -eps(i)
c              eps(i) are the eigenvalues of X**2
c              X has dimension N*N
c
c          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
c
      IMPLICIT REAL*8 (a-h,o-z)
#include "output_ras.fh"
      Parameter (Routine='EXPX')
#include "warnings.fh"

      DIMENSION x(*), x2(*), v(*), a(*), b(*), y(*)

c
c     square the matrix X into the symmetric matrix X2
c
      ij=0
      i1=1
      DO i=1,n
        j1=1
        DO j=1,i
          ij=ij+1
          x2(ij)=-ddot_(n,x(i1),1,x(j1),1)
          j1=j1+n
        END DO
        i1=i1+n
      END DO
c
c     NOW DIAGONALIZE X2
c
      call dcopy_(n*n, [0.0d0], 0, v, 1)
      call dcopy_(n, [1.0d0], 0, v, n+1)
      CALL Jacob (x2, v, n, n)
c
c     CHECK FOR SMALL VALUES OF THE EIGENVALUES
c
      ii=0
      epsmin=0.0d00
      DO i=1,n
        ii=ii+i
        eps=x2(ii)
        epsmin=min(epsmin,eps)
* PAM 2008: The following test is essentially a sanity check of the code.
* But this code has always worked nicely for years, and no longer needs
* this checking. However, the test can be inadvertently triggered occasionally
* due to extreme rounding errors. Rather than using an arbitrary threshold,
* I now deactivate this test.
*        IF (eps.gt.1.0d-8) THEN
*          Write(LF,*)
*          Write(LF,'(6X,120A1)') ('=',j=1,120)
*          Write(LF,*)
*          Write(LF,'(6X,A)') 'EXPX Error: '//
*     &    'A positive eigenvalue has been found in subroutine EXPX !!!'
*          Write(LF,'(6X,A)')
*     &    'This is possible only if there is a severe malfunction     '
*          Write(LF,'(6X,A)')
*     &    'of the program. Please issue a bug report.'
*          Write(LF,'(6X,A,E20.6)') 'eps=',eps
*          Call Quit(_RC_INTERNAL_ERROR_)
*        END IF
        IF (eps.gt.0.0D0) eps=0.0D0
        IF (eps.gt.-3.0d-7) THEN
          a(i)=eps/2.0D0
          b(i)=1.0d0+eps/6.0D0
        ELSE IF (eps.gt.-3.0D-4) THEN
          a(i)=(eps/2.0D0)*(1.0D0+eps/12.0D0)
          b(i)=1.0d0+(eps/6.0D0)*(1.0D0+eps/20.0D0)
        ELSE IF (eps.gt.-1.0D-3) THEN
          a(i)=(eps/2.0D0)*(1.0D0+(eps/12.0D0)*(1.0D0+eps/30.0D0))
          b(i)=1.0d0+(eps/6.0D0)*
     &                    (1.0D0+(eps/20.0D0)*(1.0D0+eps/42.0D0))
        ELSE
          theta=sqrt(-eps)
          a(i)=cos(theta)-1.0D0
          b(i)=sin(theta)/theta
        END IF
      END DO
      thmax=sqrt(-epsmin)
c
c     NOW FORM THE MATRIX EXP(X)=1+V*A*V(+)+V*B*V(+)*X
c     FIRST FORM THE MATRIX Y = V(+)*X
c
      CALL DGEMM_('T','N',
     &            n,n,n,
     &            1.0d0,v,n,
     &            x,n,
     &            0.0d0,y,n)
c
c     MULTIPLY Y IN PLACE WITH THE DIAGONAL MATRIX B
c
      DO i=1,n
        call dscal_(n, b(i), y(i), n)
      END DO
c
c     ADD A*V(+)
c
      CALL dnaxpy (n, n, a, 1, v, 1, n, y, n, 1)
c
c     MULTIPLY V WITH Y TO OBTAIN EXP(X)-1
c     STORE THE RESULT IN THE MATRIX X
c
      CALL DGEMM_('N','N',
     &            n,n,n,
     &            1.0d0,v,n,
     &            y,n,
     &            0.0d0,x,n)
C Finally, add 1 to obtain exp(x)
      do i=1,n**2,(n+1)
        x(i)=x(i)+1.0D0
      end do
c
      RETURN
      END
