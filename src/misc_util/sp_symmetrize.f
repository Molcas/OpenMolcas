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
*   <NAME>Sp\_Symmetrize</NAME>
*   <SYNTAX>Call Sp\_Symmetrize(n,A,ija,B,ijb)</Syntax>
*   <ARGUMENTS>
*     \Argument{n}{Size of the matrix}{Integer}{in}
*     \Argument{A}{Input matrix, in sparse format}{Real*8 (*)}{in}
*     \Argument{ija}{Index vector of matrix A}{Integer (*)}{in}
*     \Argument{B}{Output matrix, in sparse format}{Real*8 (*)}{out}
*     \Argument{ijb}{Index vector of matrix B}{Integer (*)}{out}
*   </ARGUMENTS>
*   <PURPOSE>Converts a sparse matrix to symmetric format</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Converts a matrix A to a matrix B, stored in symmetric mode.
*     Only the lower triangle of A is stored.
*   </DESCRIPTION>
* </DOC>
*-----------------------------------------------------------------------
      SUBROUTINE Sp_Symmetrize(n,A,ija,B,ijb)
      IMPLICIT NONE
      INTEGER n, nijb, ija(*), ijb(*), i, j, k
      REAL*8 A(*), B(*)
#include "real.fh"

      ijb(1)=n+2
      nijb=n+1
      DO i=1,n
        B(i)=A(i)
        DO k=ija(i),ija(i+1)-1
          j=ija(k)
          IF (j.LT.i) THEN
            nijb=nijb+1
            B(nijb)=A(k)
            ijb(nijb)=ija(k)
          END IF
        END DO
        ijb(i+1)=nijb+1
      END DO
      B(n+1)=One

      END
