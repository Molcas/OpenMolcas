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
      FUNCTION GETH2A(I,J,K,L,TUVX)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TUVX(*)
      NI=MAX(I,J)
      NJ=MIN(I,J)
      IJ=NJ+NI*(NI-1)/2
      NK=MAX(K,L)
      NL=MIN(K,L)
      KL=NL+NK*(NK-1)/2
      NIJ=MAX(IJ,KL)
      NKL=MIN(IJ,KL)
      NIJKL=NKL+NIJ*(NIJ-1)/2
      GETH2A=TUVX(NIJKL)
      RETURN
      END
