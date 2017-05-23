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
* Copyright (C) 2013,2014, Ignacio Fdez. Galvan                        *
************************************************************************
*-------------------------------------------------------------------------------
*    SUBROUTINES FOR SUPERPOSITION OF MOLECULES
*
*    I. Fdez. Galvan, Uppsala 2013
*
* Get_RMSD_w  - Find minimum possible RMSD, no change in structures
* Get_RMSD    - Ditto, with equal weights
* Superpose_w - Superpose structures, first structure modified
* Superpose   - Ditto, with equal weights
*
*--- Based on:
*    Acta Chrystallogr. Sec. A 61 (2005), 478
*    J. Comput. Chem. 31 (2010), 1561
*-------------------------------------------------------------------------------
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>Get\_RMSD\_w</NAME>
*   <SYNTAX>Call Get\_RMSD\_w(x,y,w,nAt,RMSD)</Syntax>
*   <ARGUMENTS>
*     \Argument{x}{Cartesian coordinates of the first structure}{Real*8 (3,nAt)}{in}
*     \Argument{y}{Cartesian coordinates of the second structure}{Real*8 (3,nAt)}{in}
*     \Argument{w}{Weights for each atom}{Real*8 (nAt)}{in}
*     \Argument{nAt}{Number of atoms in the structures}{Integer}{in}
*     \Argument{RMSD}{Minimum weighted RMSD between the two structures}{Real*8}{out}
*   </ARGUMENTS>
*   <PURPOSE>Compute the optimal weighted RMSD between two structures</PURPOSE>
*   <DEPENDENCIES>get\_rotation</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Compute the optimal RMSD between two structures, with weights.
*     Coordinates are not changed, but the RMSD corresponds to the aligned structures.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE Get_RMSD_w(x,y,w,nAt,RMSD)
      IMPLICIT NONE
      INTEGER nAt
      REAL*8 x(3,nAt), y(3,nAt), w(nAt), RMSD
      CALL get_rotation(x,y,w,nAt,RMSD,.FALSE.)
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>Get\_RMSD</NAME>
*   <SYNTAX>Call Get\_RMSD(x,y,nAt,RMSD)</Syntax>
*   <ARGUMENTS>
*     \Argument{x}{Cartesian coordinates of the first structure}{Real*8 (3,nAt)}{in}
*     \Argument{y}{Cartesian coordinates of the second structure}{Real*8 (3,nAt)}{in}
*     \Argument{nAt}{Number of atoms in the structures}{Integer}{in}
*     \Argument{RMSD}{Minimum RMSD between the two structures}{Real*8}{out}
*   </ARGUMENTS>
*   <PURPOSE>Compute the optimal RMSD between two structures</PURPOSE>
*   <DEPENDENCIES>Get\_RMSD\_w</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Compute the optimal RMSD between two structures, with equal weights.
*     Coordinates are not changed, but the RMSD corresponds to the aligned structures.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE Get_RMSD(x,y,nAt,RMSD)
      IMPLICIT NONE
#include "real.fh"
#include "stdalloc.fh"
      INTEGER nAt, iAt
      REAL*8 x(3,nAt), y(3,nAt), RMSD
      REAL*8, DIMENSION(:), ALLOCATABLE :: w
      CALL mma_allocate(w,nAt)
      DO iAt=1,nAt
        w(iAt)=One
      END DO
      CALL Get_RMSD_w(x,y,w,nAt,RMSD)
      CALL mma_deallocate(w)
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>Superpose\_w</NAME>
*   <SYNTAX>Call Superpose\_w(x,y,w,nAt,RMSD,RMax)</Syntax>
*   <ARGUMENTS>
*     \Argument{x}{Cartesian coordinates of the first structure}{Real*8 (3,nAt)}{inout}
*     \Argument{y}{Cartesian coordinates of the second structure}{Real*8 (3,nAt)}{in}
*     \Argument{w}{Weights for each atom}{Real*8 (nAt)}{in}
*     \Argument{nAt}{Number of atoms in the structures}{Integer}{in}
*     \Argument{RMSD}{Weighted RMSD between the two final structures}{Real*8}{out}
*     \Argument{RMax}{Maximum atom-wise (weighted) distance between the two structures}{Real*8}{out}
*   </ARGUMENTS>
*   <PURPOSE>Superpose two structures</PURPOSE>
*   <DEPENDENCIES>get\_rotation</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Superpose two structures, with weights.
*     The coordinates of the first structure are modified.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE Superpose_w(x,y,w,nAt,RMSD,RMax)
      IMPLICIT NONE
#include "real.fh"
      INTEGER nAt, iAt
      REAL*8 x(3,nAt), y(3,nAt), w(nAt), RMSD, RMax, r2
      CALL get_rotation(x,y,w,nAt,RMSD,.TRUE.)
*---- Calculate the maximum weighted distance between atoms
*     (not sure what meaning it has with w!=1)
      RMax=Zero
      Do iAt=1,nAt
        r2=(x(1,iAt)-y(1,iAt))**2+
     &     (x(2,iAt)-y(2,iAt))**2+
     &     (x(3,iAt)-y(3,iAt))**2
        RMax=MAX(RMax,w(iAt)*r2)
      END DO
      RMax=SQRT(RMax)
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>Superpose</NAME>
*   <SYNTAX>Call Superpose(x,y,nAt,RMSD,RMax)</Syntax>
*   <ARGUMENTS>
*     \Argument{x}{Cartesian coordinates of the first structure}{Real*8 (3,nAt)}{inout}
*     \Argument{y}{Cartesian coordinates of the second structure}{Real*8 (3,nAt)}{in}
*     \Argument{nAt}{Number of atoms in the structures}{Integer}{in}
*     \Argument{RMSD}{RMSD between the two final structures}{Real*8}{out}
*     \Argument{RMax}{Maximum atom-wise distance between the two structures}{Real*8}{out}
*   </ARGUMENTS>
*   <PURPOSE>Superpose two structures</PURPOSE>
*   <DEPENDENCIES>Superpose\_w</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Superpose two structures, with equal weights.
*     The coordinates of the first structure are modified.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE Superpose(x,y,nAt,RMSD,RMax)
      IMPLICIT NONE
#include "real.fh"
#include "stdalloc.fh"
      INTEGER nAt, iAt
      REAL*8 x(3,nAt), y(3,nAt), RMSD, RMax
      REAL*8, DIMENSION(:), ALLOCATABLE :: w
      CALL mma_allocate(w,nAt)
      DO iAt=1,nAt
        w(iAt)=One
      END DO
      CALL Superpose_w(x,y,w,nAt,RMSD,RMax)
      CALL mma_deallocate(w)
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>get\_rotation</NAME>
*   <SYNTAX>Call get\_rotation(x,y,w,nAt,RMSD,rotate)</Syntax>
*   <ARGUMENTS>
*     \Argument{x}{Cartesian coordinates of the first structure}{Real*8 (3,nAt)}{inout}
*     \Argument{y}{Cartesian coordinates of the second structure}{Real*8 (3,nAt)}{in}
*     \Argument{w}{Weights for each atom}{Real*8 (nAt)}{in}
*     \Argument{nAt}{Number of atoms in the structures}{Integer}{in}
*     \Argument{RMSD}{Minimum weighted RMSD between the two final structures}{Real*8}{out}
*     \Argument{rotate}{Flag for changing the coordinates or not}{Logical}{in}
*   </ARGUMENTS>
*   <PURPOSE>Compute the RMSD and superpose two structures</PURPOSE>
*   <DEPENDENCIES>center\_mol, inner\_prod, inner\_mat, build\_polynomial, find\_lambda, build\_K\_matrix, get\_eigenvalue, apply\_rotation, DScal, DCopy</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Calculate the minimum weighted RMSD between two structures and, optionally,
*     rotate and translate the first one to be superposed onto the other
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE get_rotation(x,y,w,nAt,RMSD,rotate)
      IMPLICIT NONE
#include "real.fh"
#include "stdalloc.fh"
      INTEGER nAt, iAt, i
      REAL*8 x(3,nAt), y(3,nAt), w(nAt), RMSD,
     &       wTot, Gx, Gy, c_x(3), c_y(3), Mxy(3,3), Kxy(4,4), lambda,
     &       inner_prod, q(4), DDot_
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: xCen, yCen
      LOGICAL rotate
*---- if USE_QCD is defined, use the real QCD method, otherwise use conventional
*     method to locate the largest eigenvalue (may be more robust in some cases)
*define USE_QCD
#ifdef USE_QCD
      REAL*8 Poly(5)
#else
      REAL*8 wTmp(100)
#endif
      CALL mma_allocate(xCen,3,nAt)
      CALL mma_allocate(yCen,3,nAt)
*---- Center the two structures and calculate their inner products
      CALL center_mol(x,w,nAt,c_x,xCen)
      CALL center_mol(y,w,nAt,c_y,yCen)
      Gx=inner_prod(xCen,w,nAt)
      Gy=inner_prod(yCen,w,nAt)
      CALL inner_mat(xCen,yCen,w,nAt,Mxy)
*---- Find the largest eigenvalue (always lower or equal to (Gx+Gy)/2)
#ifdef USE_QCD
      CALL build_polynomial(Mxy,Poly)
      lambda=Half*(Gx+Gy)
      CALL find_lambda(Poly,lambda)
#else
      CALL build_K_matrix(Mxy,Kxy)
      i=0
      call dsyev_('N','U',4,Kxy,4,q,wTmp,100,i)
      lambda=q(4)
#endif
*---- Calculate the optimized RMSD value
      wTot=Zero
      DO iAt=1,nAt
        wTot=wTot+w(iAt)
      END DO
      RMSD=SQRT(ABS(Gx+Gy-Two*lambda)/wTot)
*---- The rotation is only performed if an actual superposition is requested
      IF (rotate) THEN
*----   The rotation quaternion is obtained as the eigenvector of K corresponding
*       to the eigenvalue lambda.
        CALL build_K_matrix(Mxy,Kxy)
        CALL get_eigenvector(Kxy,lambda,q)
*----   Apply the rotation to the first structure. Before rotation, the quaternion
*       must be normalized and the rotation inverted (change of sign in the 1st component)
        call dscal_(4,One/SQRT(DDot_(4,q,1,q,1)),q,1)
        q(1)=-q(1)
        CALL apply_rotation(xCen,nAt,q)
*----   Translate the structure to match the second (reference) one
        DO iAt=1,nAt
          DO i=1,3
            xCen(i,iAt)=xCen(i,iAt)+c_y(i)
          END DO
        END DO
*----   Copy the aligned structures back in the original matrices
        call dcopy_(3*nAt,xCen,1,x,1)
      END IF
      CALL mma_deallocate(xCen)
      CALL mma_deallocate(yCen)
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>center\_mol</NAME>
*   <SYNTAX>Call center\_mol(x,w,nAt,c,xCen)</Syntax>
*   <ARGUMENTS>
*     \Argument{x}{Cartesian coordinates of the structure}{Real*8 (3,nAt)}{in}
*     \Argument{w}{Weights for each atom}{Real*8 (nAt)}{in}
*     \Argument{nAt}{Number of atoms in the structures}{Integer}{in}
*     \Argument{c}{Weighted center}{Real*8 (3)}{out}
*     \Argument{xCen}{Centered Cartesian coordinates}{Real*8 (3,nAt)}{out}
*   </ARGUMENTS>
*   <PURPOSE>Center a structure on its weighted geometric center</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Center a structure on its weighted geometric center.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE center_mol(x,w,nAt,c,xCen)
      IMPLICIT NONE
#include "real.fh"
      INTEGER nAt, iAt, i
      REAL*8 x(3,nAt), w(nAt), c(3), xCen(3,nAt), wTot
*---- Compute the total weight
      wTot=Zero
      DO iAt=1,nAt
        wTot=wTot+w(iAt)
      END DO
      DO i=1,3
*----   Calculate the center in each dimension (x, y, z)
        c(i)=Zero
        DO iAt=1,nAt
          c(i)=c(i)+w(iAt)*x(i,iAt)
        END DO
        c(i)=c(i)/wTot
*----   Center the structure in each dimension
        DO iAt=1,nAt
          xCen(i,iAt)=x(i,iAt)-c(i)
        END DO
      END DO
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>inner\_prod</NAME>
*   <SYNTAX>inner\_prod(x,w,nAt)</Syntax>
*   <ARGUMENTS>
*     \Argument{x}{Cartesian coordinates of the structure}{Real*8 (3,nAt)}{in}
*     \Argument{w}{Weights for each atom}{Real*8 (nAt)}{in}
*     \Argument{nAt}{Number of atoms in the structures}{Integer}{in}
*   </ARGUMENTS>
*   <PURPOSE>Calculate the inner product of a single structure</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Real*8 function that returns the inner product of the weighted coordinates of a single structure.
*     $\mbox{inner\_prod} = \sum_{i=1,N} w_i ( x_i^2 + y_i^2 + z_i^2 )$.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      FUNCTION inner_prod(x,w,nAt)
      IMPLICIT NONE
#include "real.fh"
      INTEGER nAt, iAt
      REAL*8 inner_prod, x(3,nAt), w(nAt), r
      r=Zero
      DO iAt=1,nAt
        r=r+w(iAt)*(x(1,iAt)**2+x(2,iAt)**2+x(3,iAt)**2)
      END DO
      inner_prod=r
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>inner\_mat</NAME>
*   <SYNTAX>Call inner\_mat(x,y,w,nAt,M)</Syntax>
*   <ARGUMENTS>
*     \Argument{x}{Cartesian coordinates of the first structure}{Real*8 (3,nAt)}{in}
*     \Argument{x}{Cartesian coordinates of the second structure}{Real*8 (3,nAt)}{in}
*     \Argument{w}{Weights for each atom}{Real*8 (nAt)}{in}
*     \Argument{nAt}{Number of atoms in the structures}{Integer}{in}
*     \Argument{M}{Inner product matrix}{Real*8 (3,3)}{out}
*   </ARGUMENTS>
*   <PURPOSE>Calculate the inner product matrix of two structures structure</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Calculate the inner product matrix of the weighted coordinates of two structures.
*     $M_{ij} = \sum_{k=1,N} w_i (A_i)_k (B_j)_k$, where $a,b = \{x, y, z\}$.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE inner_mat(x,y,w,nAt,M)
      IMPLICIT NONE
#include "real.fh"
      INTEGER nAt, iAt, i, j
      REAL*8 x(3,nAt), y(3,nAt), w(nAt), M(3,3)
      DO i=1,3
        DO j=1,3
          M(j,i)=Zero
          DO iAt=1,nAt
            M(j,i)=M(j,i)+w(iAt)*x(j,iAt)*y(i,iAt)
          END DO
        END DO
      END DO
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>build\_K\_matrix</NAME>
*   <SYNTAX>Call build\_K\_matrix(M,K)</Syntax>
*   <ARGUMENTS>
*     \Argument{M}{Inner product matrix}{Real*8 (3,3)}{in}
*     \Argument{K}{Key matrix K}{Real*8 (4,4)}{out}
*   </ARGUMENTS>
*   <PURPOSE>Build the key matrix K from the M matrix</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Build the key matrix K from the inner product matrix M.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE build_K_matrix(M,K)
      IMPLICIT NONE
      INTEGER i, j
      REAL*8 M(3,3), K(4,4)
*---- Compute the unique elements
      K(1,1)=M(1,1)+M(2,2)+M(3,3)
      K(1,2)=M(2,3)-M(3,2)
      K(1,3)=M(3,1)-M(1,3)
      K(1,4)=M(1,2)-M(2,1)
      K(2,2)=M(1,1)-M(2,2)-M(3,3)
      K(2,3)=M(1,2)+M(2,1)
      K(2,4)=M(1,3)+M(3,1)
      K(3,3)=M(2,2)-M(1,1)-M(3,3)
      K(3,4)=M(2,3)+M(3,2)
      K(4,4)=M(3,3)-M(1,1)-M(2,2)
*---- Make the matrix symmetric
      DO i=2,4
        DO j=1,i-1
          K(i,j)=K(j,i)
        END DO
      END DO
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>build\_polynomial</NAME>
*   <SYNTAX>Call build\_polynomial(M,P)</Syntax>
*   <ARGUMENTS>
*     \Argument{M}{Inner product matrix}{Real*8 (3,3)}{in}
*     \Argument{P}{The polynomial coefficients}{Real*8 (5)}{out}
*   </ARGUMENTS>
*   <PURPOSE>Calculate the coefficients for the characteristic polynomial</PURPOSE>
*   <DEPENDENCIES>determinant3, DDot</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Calculate the coefficients for the characteristic polynomial of the K matrix
*     directly from the elements of the M matrix.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE build_polynomial(M,P)
      IMPLICIT NONE
#include "real.fh"
      REAL*8 M(3,3), P(5), a1, a2, a3, a4, a5, a6, DDot_, determinant3
*---- The x^4 and x^3 coefficients are constant
      P(5)=One
      P(4)=Zero
*---- The x^2 coefficient is the sum of all squared elements, times -2
      P(3)=-Two*DDot_(9,M,1,M,1)
*---- The x^1 coefficient is the determinant of M, times -8
      P(2)=-Eight*determinant3(M)
*---- The x^0 coefficient has a rather cumbersome expansion
      a1=(M(1,2)**2+M(1,3)**2-M(2,1)**2-M(3,1)**2)**2
      a2=(-M(1,1)**2+M(2,2)**2+M(3,3)**2+M(2,3)**2+M(3,2)**2-
     &    Two*(M(2,2)*M(3,3)-M(2,3)*M(3,2)))*
     &   (-M(1,1)**2+M(2,2)**2+M(3,3)**2+M(2,3)**2+M(3,2)**2+
     &    Two*(M(2,2)*M(3,3)-M(2,3)*M(3,2)))
      a3=(-(M(1,3)+M(3,1))*(M(2,3)-M(3,2))+
     &     (M(1,2)-M(2,1))*(M(1,1)-M(2,2)-M(3,3)))*
     &   (-(M(1,3)-M(3,1))*(M(2,3)+M(3,2))+
     &     (M(1,2)-M(2,1))*(M(1,1)-M(2,2)+M(3,3)))
      a4=(-(M(1,3)+M(3,1))*(M(2,3)+M(3,2))-
     &     (M(1,2)+M(2,1))*(M(1,1)+M(2,2)-M(3,3)))*
     &   (-(M(1,3)-M(3,1))*(M(2,3)-M(3,2))-
     &     (M(1,2)+M(2,1))*(M(1,1)+M(2,2)+M(3,3)))
      a5=( (M(1,2)+M(2,1))*(M(2,3)+M(3,2))+
     &     (M(1,3)+M(3,1))*(M(1,1)-M(2,2)+M(3,3)))*
     &   (-(M(1,2)-M(2,1))*(M(2,3)-M(3,2))+
     &     (M(1,3)+M(3,1))*(M(1,1)+M(2,2)+M(3,3)))
      a6=( (M(1,2)+M(2,1))*(M(2,3)-M(3,2))+
     &     (M(1,3)-M(3,1))*(M(1,1)-M(2,2)-M(3,3)))*
     &   (-(M(1,2)-M(2,1))*(M(2,3)+M(3,2))+
     &     (M(1,3)-M(3,1))*(M(1,1)+M(2,2)-M(3,3)))
      P(1)=a1+a2+a3+a4+a5+a6
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>find\_lambda</NAME>
*   <SYNTAX>Call find\_lambda(P,r)</Syntax>
*   <ARGUMENTS>
*     \Argument{P}{The polynomial coefficients}{Real*8 (5)}{in}
*     \Argument{r}{The found root}{Real*8)}{inout}
*   </ARGUMENTS>
*   <PURPOSE>Find one root of a 4th degree polynomial</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Find one real root of a 4th degree polynomial using the Newton-Raphson method.
*     On input, the variable r contains the initial guess for the root.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE find_lambda(P,r)
      IMPLICIT NONE
#include "real.fh"
      REAL*8 P(5), r, thr, thrr, r_old, f, df
      INTEGER loop, mxloop, j
      PARAMETER ( thr=1.0D-11, mxloop=100 )
      r_old=MAX(Two*r,Ten)
*---- Find the root with the Newton-Raphson method, starting from the initial value
      loop=0
      thrr=thr*r
      DO WHILE ((ABS(r-r_old).GT.thrr).AND.(loop.LT.mxloop))
        r_old=r
*----   Compute the values of the polynomial (f) and its derivative (df)
*       at the current value of r
        f=P(5)
        df=Zero
        DO j=4,1,-1
          df=df*r+f
          f=f*r+P(j)
        END DO
*----   Update the value of r from the values of f and df
*       (special case if the derivative is zero)
        IF (ABS(df).LT.thrr) THEN
          IF (ABS(f).LT.thr) THEN
            r_old=r
          ELSE
*---        If the df is zero and r is not a root, we are in trouble,
*           just move the point a tad
            r=r-SIGN(Two*thrr,f)
          END IF
        ELSE
          r=r-f/df
        END IF
        thrr=thr*r
*----   Increase count to avoid neverending loop
        loop=loop+1
      END DO
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>get\_eigenvector</NAME>
*   <SYNTAX>Call get\_eigenvector(K,r,V)</Syntax>
*   <ARGUMENTS>
*     \Argument{K}{The K matrix}{Real*8 (4,4)}{in}
*     \Argument{r}{The eigenvalue}{Real*8)}{in}
*     \Argument{V}{The eigenvector}{Real*8 (4)}{out}
*   </ARGUMENTS>
*   <PURPOSE>Find an eigenvector of the K matrix, given its eigenvalue</PURPOSE>
*   <DEPENDENCIES>cofactor, DDot</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS>The matrix K is destroyed on output</SIDE_EFFECTS>
*   <DESCRIPTION>
*     Compute the eigenvector of the K matrix corresponding to one eigenvalue.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE get_eigenvector(K,r,V)
      IMPLICIT NONE
#include "real.fh"
      REAL*8 K(4,4), r, V(4), n, thr, DDot_, cofactor
      INTEGER i, j
      PARAMETER ( thr=1.0D-12 )
*---- Get the matrix K-rI
      K(1,1)=K(1,1)-r
      K(2,2)=K(2,2)-r
      K(3,3)=K(3,3)-r
      K(4,4)=K(4,4)-r
*---- Find the eigenvector as any non-zero column of adj(K-rI)
      n=Zero
      DO j=1,4
        IF (n.LT.thr) THEN
*----     Since adj(X) is the transpose of the cofactors matrix, we take a row from K
          DO i=1,4
            V(i)=cofactor(K,j,i)
          END DO
*----     Get the norm of the vector, if not too small, no further vector will be computed
          n=DDot_(4,V,1,V,1)
        END IF
      END DO
*---- If the norm is still too small, return the identity quaternion (no rotation)
      IF (n.LT.thr) THEN
        V(1)=One
        V(2)=Zero
        V(3)=Zero
        V(4)=Zero
      END IF
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>cofactor</NAME>
*   <SYNTAX>cofactor(M,i,j)</Syntax>
*   <ARGUMENTS>
*     \Argument{M}{The matrix}{Real*8 (4,4)}{in}
*     \Argument{i}{The row of the cofactor}{Integer}{in}
*     \Argument{j}{The column of the cofactor}{Integer}{in}
*   </ARGUMENTS>
*   <PURPOSE>Return a cofactor from a 4x4 matrix</PURPOSE>
*   <DEPENDENCIES>determinant3</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Real*8 function that returns a cofactor (determinant of a submatrix) from
*     a 4x4 matrix.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      FUNCTION cofactor(M,i,j)
      IMPLICIT NONE
#include "real.fh"
      REAL*8 cofactor, M(4,4), A(3,3), f, determinant3
      INTEGER i, j, ii, jj
*---- Get the submatrix as four blocks, depending on whether the indices are
*     greater or smaller than the element position
      DO ii=1,i-1
        DO jj=1,j-1
          A(ii,jj)=M(ii,jj)
        END DO
      END DO
      DO ii=1,i-1
        DO jj=j+1,4
          A(ii,jj-1)=M(ii,jj)
        END DO
      END DO
      DO ii=i+1,4
        DO jj=1,j-1
          A(ii-1,jj)=M(ii,jj)
        END DO
      END DO
      DO ii=i+1,4
        DO jj=j+1,4
          A(ii-1,jj-1)=M(ii,jj)
        END DO
      END DO
*---- The cofactor is +1/-1 times the determinant of the submatrix
      f=(-One)**(i+j)
      cofactor=f*determinant3(A)
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>determinant3</NAME>
*   <SYNTAX>determinant3(M)</Syntax>
*   <ARGUMENTS>
*     \Argument{M}{The matrix}{Real*8 (3,3)}{in}
*   </ARGUMENTS>
*   <PURPOSE>Return the determinant of a 3x3 matrix</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Real*8 function that returns the determinant of a 3x3 matrix.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      FUNCTION determinant3(M)
      IMPLICIT NONE
      REAL*8 determinant3, M(3,3)
      determinant3=(M(1,1)*M(2,2)*M(3,3)+
     &              M(1,2)*M(2,3)*M(3,1)+
     &              M(1,3)*M(2,1)*M(3,2))-
     &             (M(1,1)*M(2,3)*M(3,2)+
     &              M(1,2)*M(2,1)*M(3,3)+
     &              M(1,3)*M(2,2)*M(3,1))
      END
*
*-------------------------------------------------------------------------------
* <DOC>
*   <NAME>apply\_rotation</NAME>
*   <SYNTAX>Call apply\_rotation(x,nAt,q)</Syntax>
*   <ARGUMENTS>
*     \Argument{x}{Cartesian coordinates of the structure}{Real*8 (3,nAt)}{inout}
*     \Argument{nAt}{Number of atoms in the structure}{Integer}{in}
*     \Argument{q}{A unit quaternion representing a rotation}{Real*8 (4)}{in}
*   </ARGUMENTS>
*   <PURPOSE>Rotate a structure with a quaternion</PURPOSE>
*   <DEPENDENCIES>DDot, DCopy</DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Apply a rotation, given by a unit quaternion, to all atoms in a structure.
*   </DESCRIPTION>
* </DOC>
*-------------------------------------------------------------------------------
      SUBROUTINE apply_rotation(x,nAt,q)
      IMPLICIT NONE
#include "real.fh"
      INTEGER nAt, iAt, i
      REAL*8 x(3,nAt), q(4), MRot(3,3), v(3), DDot_
*---- Compute the rotation matrix from the quaternion
*     (transposed or not, that depends on the rotation direction
*      and storage order... this version seems to work)
      MRot(1,1)=q(1)**2+q(2)**2-q(3)**2-q(4)**2
      MRot(2,2)=q(1)**2-q(2)**2+q(3)**2-q(4)**2
      MRot(3,3)=q(1)**2-q(2)**2-q(3)**2+q(4)**2
      MRot(2,1)=Two*(q(2)*q(3)+q(1)*q(4))
      MRot(1,2)=Two*(q(2)*q(3)-q(1)*q(4))
      MRot(3,2)=Two*(q(3)*q(4)+q(1)*q(2))
      MRot(2,3)=Two*(q(3)*q(4)-q(1)*q(2))
      MRot(1,3)=Two*(q(2)*q(4)+q(1)*q(3))
      MRot(3,1)=Two*(q(2)*q(4)-q(1)*q(3))
*---- Apply the matrix to every atom in the structure
      DO iAt=1,nAt
        DO i=1,3
          v(i)=DDot_(3,MRot(1,i),1,x(1,iAt),1)
        END DO
        call dcopy_(3,v,1,x(1,iAt),1)
      END DO
      END
