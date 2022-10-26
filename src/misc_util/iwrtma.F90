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
! Copyright (C) 1987, Jeppe Olsen                                      *
!               1989, Markus P. Fuelscher                              *
!***********************************************************************
      SUBROUTINE IWRTMA(IMAT,NROW,NCOL,MAXROW,MAXCOL)
!
!     AUTHOR: J. OLSEN, UNIV. OF LUND, SWEDEN, APRIL 1987
!     MODIFICATIONS: INCLUSION INTO THE RASSCF METHOD
!                    M.P. FUELSCHER, UNIV. OF LUND, SWEDEN, MAY 1989
!
!     WRITE A MATRIX OF INTEGER VALUES
!
      DIMENSION IMAT(MAXROW,MAXCOL)
!
      DO 100 I = 1, NROW
        Write(6,1110) (IMAT(I,J),J= 1,NCOL)
100   CONTINUE
1110  FORMAT(/,1X,8I10,/,(1X,8I10))
!
      RETURN
      END
