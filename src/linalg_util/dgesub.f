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
      Subroutine dgeSub(A,LDA,FORMA,B,LDB,FORMB,C,LDC,M,N)
C
C     MATRIX subtraction FOR GENERAL MATRICES
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 FORMA,FORMB
      REAL*8 A(*), B(*), C(*)

      IF (FORMA.EQ.'N' .AND. FORMB.EQ.'N') THEN
         Do 10 iRow=0,m-1
         Do 10 iCol=0,n-1
            c(iRow+iCol*ldc+1)=a(iRow+iCol*lda+1)-b(iRow+iCol*ldb+1)
   10    Continue
      ELSE IF (FORMA.EQ.'T' .AND. FORMB.EQ.'N') THEN
         Do 20 iRow=0,m-1
         Do 20 iCol=0,n-1
            c(iRow+iCol*ldc+1)=a(iCol+iRow*lda+1)-b(iRow+iCol*ldb+1)
   20    Continue
      ELSE IF (FORMA.EQ.'N' .AND. FORMB.EQ.'T') THEN
         Do 30 iRow=0,m-1
         Do 30 iCol=0,n-1
            c(iRow+iCol*ldc+1)=a(iRow+iCol*lda+1)-b(iCol+iRow*ldb+1)
   30    Continue
      ELSE IF (FORMA.EQ.'T' .AND. FORMB.EQ.'T') THEN
         Do 40 iRow=0,m-1
         Do 40 iCol=0,n-1
            c(iRow+iCol*ldc+1)=a(iCol+iRow*lda+1)-b(iCol+iRow*ldb+1)
   40    Continue
      ELSE
         Write(6,*) 'Error when calling DGESUB, forma=',
     &     FormA,'   formb=',FormB
         Call Abend()
      END IF
      RETURN
      END
