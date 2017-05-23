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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MLTUNF(LST,X,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),Y(*)
      DIMENSION LST(4,NLST1)
#include "sigma.fh"

C Given a list with entries LST(4,ITEM), ITEM=1,NLST1,
C with entries called L1,L2,L3,L4 for given ITEM, and
C an array of the form Y(p,q), compute the matrix
C    X(p,L1,L2) := Add V*Y(p,L3), p=1..LEN1
C where V=VAL1(L4), looped over ITEM=1,NLST1.
C Note: Arrays are addressed by strides given in common.
*     write(*,*)' In MLTUNF. List:'
*     write(*,'(1x,4i5)') ((LST(I,J),I=1,4),J=1,NLST1)
*     write(*,'(1x,a,3i5)')'INCX:',INCX1,INCX2,INCX3
*     write(*,'(1x,a,3i5)')'INCY:',INCY1,INCY2
*     write(*,'(1x,a,3i5)')'LEN1:',LEN1
      DO ILST=1,NLST1
        L1=LST(1,ILST)
        L2=LST(2,ILST)
        L3=LST(3,ILST)
        L4=LST(4,ILST)
        V=VAL1(L4)
        IX=1+INCX2*(L1-1)+INCX3*(L2-1)
        IY=1+INCY2*(L3-1)
        CALL DAXPY_(LEN1,V,Y(IY),INCY1,X(IX),INCX1)
      END DO
      RETURN
      END
