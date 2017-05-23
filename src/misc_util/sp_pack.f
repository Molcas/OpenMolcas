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
*   <NAME>Sp\_Pack</NAME>
*   <SYNTAX>Call Sp\_Pack(n,A,nmax,Sp,ij\_Sp,Sym,Thr)</Syntax>
*   <ARGUMENTS>
*     \Argument{n}{Size of the matrix}{Integer}{in}
*     \Argument{A}{Input matrix, in dense (conventional) format}{Real*8 (n,n)}{in}
*     \Argument{nmax}{Maximum size of the stored vectors}{Integer}{in}
*     \Argument{Sp}{Output matrix in sparse format}{Real*8 (*)}{out}
*     \Argument{ij\_Sp}{Index vector for the sparse output matrix}{Integer (*)}{out}
*     \Argument{Sym}{Flag specifying whether to use the symmetric format}{Logical}{in}
*     \Argument{Thr}{Threshold to discard small elements}{Real*8}{in}
*   </ARGUMENTS>
*   <PURPOSE>Store a matrix in sparse format</PURPOSE>
*   <DEPENDENCIES></DEPENDENCIES>
*   <AUTHOR>I. Fdez. Galvan</AUTHOR>
*   <MODIFIED_BY></MODIFIED_BY>
*   <SIDE_EFFECTS></SIDE_EFFECTS>
*   <DESCRIPTION>
*     An input matrix A will be stored in a sparse format, as two vectors Sp and ij\_Sp.
*     Only elements with an absolute value larger than Thr will be kept.
*     An error is raised if the number of stored elements is larger than nmax.
*     The sizes of Sp and ij\_Sp must be at least n+k+1, where k is the number of
*     non_zero off-diagonal elements of A.
*   </DESCRIPTION>
* </DOC>
*-----------------------------------------------------------------------
      SUBROUTINE Sp_Pack(n,A,nmax,Sp,ij_Sp,Sym,Thr)
      IMPLICIT NONE
      INTEGER n, nmax, ij_Sp(*), i, j, nij
      REAL*8 A(n,n), Sp(*), Thr
      LOGICAL Sym
#include "real.fh"

      ij_Sp(1)=n+2
      nij=n+1
      IF (Sym) THEN
        DO i=1,n
          Sp(i)=A(i,i)
          DO j=1,i-1
            IF (ABS(A(i,j)).LE.Thr) THEN
              nij=nij+1
              IF (nij.GT.nmax) THEN
                CALL SysAbendMsg('Sp_Pack',
     &          'Too many non-zero elements.','')
              END IF
              Sp(nij)=A(i,j)
              ij_Sp(nij)=j
            END IF
          END DO
          ij_Sp(i+1)=nij+1
        END DO
        Sp(n+1)=One
      ELSE
        DO i=1,n
          Sp(i)=A(i,i)
          DO j=1,n
            IF ((j.NE.i).AND.(ABS(A(i,j)).LE.Thr)) THEN
              nij=nij+1
              IF (nij.GT.nmax) THEN
                CALL SysAbendMsg('Sp_Pack',
     &          'Too many non-zero elements.','')
              END IF
              Sp(nij)=A(i,j)
              ij_Sp(nij)=j
            END IF
          END DO
          ij_Sp(i+1)=nij+1
        END DO
        Sp(n+1)=Zero
      END IF

      END
