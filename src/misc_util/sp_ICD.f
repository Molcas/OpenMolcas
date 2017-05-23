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
*   <NAME>Sp\_ICD</NAME>
*   <SYNTAX>Call Sp\_ICD(n,A,ija,B,ijb)</Syntax>
*   <ARGUMENTS>
*     \Argument{n}{Size of the square matrices}{Integer}{in}
*     \Argument{A}{Input matrix in sparse format}{Real*8 (*)}{in}
*     \Argument{ija}{Index vector of matrix A}{Integer (*)}{in}
*     \Argument{B}{Output matrix in sparse format}{Real*8 (*)}{out}
*     \Argument{ijb}{Index vector of matrix B}{Integer (*)}{out}
*   </ARGUMENTS>
*   <PURPOSE>Perform an incomplete Cholesky decomposition of sparse matrix A</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Perform an incomplete Cholesky decomposition of a symmetric matrix A, in sparse format.
*     On output, matrix B contains a lower triangular matrix such that:
*       $ A \simeq B B^T$
*     ``Incomplete'' means that only non-zero elements in A will be computed in B.
*     This is useful as a preconditioner.
*     The decomposition can be done in place, if B=A, but only if A is stored as symmetric.
*   </DESCRIPTION>
* </DOC>
*-----------------------------------------------------------------------
      SUBROUTINE Sp_ICD(n,A,ija,B,ijb)
      IMPLICIT NONE
      INTEGER n, ija(*), ijb(*), nijb, i, j, k, kk, kkb, l
      REAL*8 A(*), B(*), Ljk, Thr
      INTEGER idLoc
      EXTERNAL idLoc
      PARAMETER (Thr=1.0D-12)
      LOGICAL Sym, GoOn
#include "real.fh"

      Sym=(A(n+1).GT.0.0D0)
      IF (idLoc(A(1)).EQ.idLoc(B(1))) THEN
        IF (.NOT.Sym) THEN
          CALL SysAbendMsg('Sp_ICD',
     &                     'In-place decomposition only allowed with '
     &                   //'symmetric-stored matrix.','')
        END IF
      END IF
      nijb=n+1
      ijb(1)=n+2
      DO i=1,n
        B(i)=A(i)
c
c       Loop all elements in row i
        DO k=ija(i),ija(i+1)-1
          j=ija(k)
          IF (j.LT.i) THEN
            nijb=nijb+1
            B(nijb)=A(k)
            ijb(nijb)=ija(k)
c
c           Loop all previous elements of row i
            DO kk=ija(i),k-1
              Ljk=Zero
              GoOn=.TRUE.
              l=ijb(j)
c
c             This loop to find an element in row j that belongs
c             to the same column as each of the parent loop
              DO WHILE (GoOn)
                IF (ijb(l).GE.j) GoOn=.FALSE.
                IF (ijb(l).EQ.ija(kk)) THEN
                  Ljk=B(l)
                  GoOn=.FALSE.
                END IF
                l=l+1
                IF (l.GE.ijb(j+1)) GoOn=.FALSE.
              END DO
c
c             As the number of elements per row in A and B can be
c             different, the proper offset must be calculated
              kkb=kk-ija(i)+ijb(i)
              B(nijb)=B(nijb)-B(kkb)*Ljk
            END DO
            IF (B(j).GT.Thr) THEN
              B(nijb)=B(nijb)/B(j)
            ELSE
              B(nijb)=Zero
            END IF
            B(i)=B(i)-B(nijb)**2
          END IF
        END DO
        B(i)=SQRT(ABS(B(i)))
        ijb(i+1)=nijb+1
      END DO
c
c     The lower triangular matrix is not symmetric, obviously
      B(n+1)=Zero

      END
