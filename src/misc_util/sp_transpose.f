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
*   <NAME>Sp\_Transpose</NAME>
*   <SYNTAX>Call Sp\_Transpose(n,A,ija,B,ijb,nij)</Syntax>
*   <ARGUMENTS>
*     \Argument{n}{Size of the matrices}{Integer}{in}
*     \Argument{A}{Input matrix, in sparse format}{Real*8 (*)}{in}
*     \Argument{ija}{Index vector of matrix A}{Integer (*)}{in}
*     \Argument{B}{Output matrix, in sparse format}{Real*8 (*)}{out}
*     \Argument{ijb}{Index vector of matrix B}{Integer (*)}{out}
*     \Argument{nij}{Length of the vectors}{Integer}{in}
*   </ARGUMENTS>
*   <PURPOSE>Transposes a matrix in sparse format</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Saves in B the transpose of the input matrix A, both in sparse format.
*   </DESCRIPTION>
* </DOC>
*-----------------------------------------------------------------------
      SUBROUTINE Sp_Transpose(n,A,ija,B,ijb,nij)
      IMPLICIT NONE
      INTEGER n, ija(*), ijb(*), nij, i, j, k, kk
      INTEGER, DIMENSION(:), ALLOCATABLE :: ia
      REAL*8 A(*), B(*)
#include "real.fh"
#include "stdalloc.fh"

      IF (A(n+1).GT.Zero) THEN
        call dcopy_(nij,A,1,B,1)
        CALL ICopy(nij,ija,1,ijb,1)
      ELSE
        CALL mma_allocate(ia,nij)
c
c       Create an index of the rows in A
        DO i=1,n
          B(i)=A(i)
          DO k=ija(i),ija(i+1)-1
            ia(k)=i
          END DO
        END DO
c
c       Lookup each column in A, save in B as a row
        ijb(1)=n+2
        kk=ijb(1)
        DO j=1,n
          DO k=ija(1),nij
            IF (ija(k).EQ.j) THEN
              ijb(kk)=ia(k)
              B(kk)=A(k)
              kk=kk+1
            END IF
          END DO
          ijb(j+1)=kk
        END DO
        B(n+1)=Zero
        CALL mma_deallocate(ia)
      END IF

      END
