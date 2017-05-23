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
*   <NAME>Sp\_Unpack</NAME>
*   <SYNTAX>Call Sp\_Unpack(n,Sp,ij\_Sp,A)</Syntax>
*   <ARGUMENTS>
*     \Argument{n}{Size of the matrix}{Integer}{in}
*     \Argument{Sp}{Input matrix, in sparse format}{Real*8 (*)}{in}
*     \Argument{ij\_Sp}{Index vector of matrix Sp}{Integer (*)}{in}
*     \Argument{A}{Output matrix, in dense (conventional) format}{Real*8 (n,n)}{out}
*   </ARGUMENTS>
*   <PURPOSE>Converts a sparse matrix to a dense (conventional) one</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     Unpacks the sparse matrix Sp into a conventional dense matrix A.
*   </DESCRIPTION>
* </DOC>
*-----------------------------------------------------------------------
      SUBROUTINE Sp_Unpack(n,Sp,ij_Sp,A)
      IMPLICIT NONE
      INTEGER n, ij_Sp(*), i, j, k
      REAL*8 Sp(*), A(n,n)
#include "real.fh"

      CALL DZero(n*n,A)
      IF (Sp(n+1).GT.Zero) THEN
        DO i=1,n
          A(i,i)=Sp(i)
          DO k=ij_Sp(i),ij_Sp(i+1)-1
            j=ij_Sp(k)
            A(i,j)=Sp(k)
            A(j,i)=Sp(k)
          END DO
        END DO
      ELSE
        DO i=1,n
          A(i,i)=Sp(i)
          DO k=ij_Sp(i),ij_Sp(i+1)-1
            j=ij_Sp(k)
            A(i,j)=Sp(k)
          END DO
        END DO
      END IF

      END
