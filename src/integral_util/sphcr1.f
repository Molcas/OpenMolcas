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
* Copyright (C) 1990,1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine SphCr1(Win,ijkla,
     &                  Scrt,nScrt,
     &                  Coeff3,kCar,kSph,Tr3,Pr3,
     &                  Coeff4,lCar,lSph,Tr4,Pr4,Wout,mcd)
************************************************************************
*                                                                      *
* Object : to transform the two-electron integrals from cartesian      *
*          gaussians to real spherical harmonic gaussians.             *
*                                                                      *
*          Observe that most of the time Win and Wout will overlap.    *
*                                                                      *
* Called from: TwoEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              DGEMM_   (ESSL)                                         *
*              RecPrt                                                  *
*              DGeTMO   (ESSL)                                         *
*              DCopy    (ESSL)                                         *
*              GetMem                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN.                                         *
*             Modified to back projection to cartesian gaussians,      *
*             January '92.                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
      Real*8 Win(ijkla*kSph*lSph), Scrt(nScrt),
     &       Coeff3(kCar,kCar), Coeff4(lCar,lCar),
     &       Wout(mcd*ijkla)
      Logical Tr3, Pr3, Tr4, Pr4
*
      iRout = 59
      iPrint = nPrint(iRout)
*     iQ = 0
*     Call qEnter('SphCr1')
*     Call RecPrt(' In SphCr1: P(AB|CD) ',' ',Win,ijkla,kSph*lSph)
      If (Tr3.and.Tr4) Then
*        Call RecPrt(' Right contraction',' ',Coeff4,lCar,lSph)
*--------Starting with IJKL,AB,CD transforming to d,IJKL,AB,C
*        Call xxDGeMul(Coeff4,lCar,'N',
*    &               Win,ijkla*kSph,'T',
*    &               Scrt,lCar,
*    &               lCar,lSph,ijkla*kSph)
         Call NTMul(Coeff4,Win,Scrt,lCar,lSph,ijkla*kSph)
*
*        Call RecPrt(' In SphCr: P(AB|Cd) ',' ',Scrt,lCar*ijkla,kSph)
*        Call RecPrt(' Left contraction',' ',Coeff3,kCar,kSph)
*--------Transform d,IJKL,AB,C to cd,IJKL,AB
*        Call xxDGeMul(Coeff3,kCar,'N',
*    &               Scrt,lCar*ijkla,'T',
*    &               Wout,kCar,
*    &               kCar,kSph,lCar*ijkla)
         Call NTMul(Coeff3,Scrt,Wout,kCar,kSph,lCar*ijkla)
      Else If (Tr4) Then
*        Call RecPrt(' Right contraction',' ',Coeff4,lCar,lSph)
*--------Starting with IJKL,AB,cD transforming to d,IJKL,AB,c
*        Call xxDGeMul(Coeff4,lCar,'N',
*    &               Win,ijkla*kCar,'T',
*    &               Scrt,lCar,
*    &               lCar,lSph,ijkla*kCar)
         Call NTMul(Coeff4,Win,Scrt,lCar,lSph,ijkla*kCar)
*--------Transpose d,IJKL,AB,c to cd,IJKL,AB
         Call DGeTMO(Scrt,lCar*ijkla,lCar*ijkla,kCar,Wout,kCar)
      Else If (Tr3) Then
*--------Transpose IJKL,AB,C,d to d,IJKL,AB,C
         Call DGeTMO(Win,ijkla*kSph,ijkla*kSph,lCar,Scrt,lCar)
*
*        Call RecPrt(' Left contraction',' ',Coeff3,kCar,kSph)
*        Transform d,IJKL,AB,c to cd,IJKL,AB
*        Call xxDGeMul(Coeff3,kCar,'N',
*    &               Scrt,lCar*ijkla,'T',
*    &               Wout,kCar,
*    &               kCar,kSph,lCar*ijkla)
         Call NTMul(Coeff3,Scrt,Wout,kCar,kSph,lCar*ijkla)
      Else
*---------Transpose IJKL,AB,cd to cd,IJKL,AB
          If (kCar*lCar.ne.1) Then
             call dcopy_(ijkla*kCar*lCar,Win,1,Scrt,1)
             Call DGeTMO(Scrt,ijkla,ijkla,kCar*lCar,Wout,kCar*lCar)
          Else
             call dcopy_(ijkla*kCar*lCar,Win,1,Scrt,1)
             call dcopy_(ijkla*kCar*lCar,Scrt,1,Wout,1)
          End If
      End If
*
*     Call RecPrt(' In SphCr1: P(AB|cd)  ',' ',Wout,mcd,ijkla)
*     Call GetMem(' Exit SphCr1','CHECK','REAL',iDum,iDum)
*     Call qExit('SphCr1')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_logical(Pr3)
         Call Unused_logical(Pr4)
      End If
      End
#ifdef _OLD_CODE_
      Subroutine NTMul(A,B,C,nRowA,nColA,nRowB)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 A(nRowA,nColA), B(nRowB,nColA),
     &       C(nRowA,nRowB)
*
      nCache=(64/8)*1024
      mCache=(nCache*3)/4 - nRowA*nColA
      Incj=mCache/(nRowA+nColA)
*
*-----Sectioning of long index
*
      Do jj = 1, nRowB, Incj
         njVec=Min(Incj,nRowB-jj+1)
*
         Do i = 1, nRowA
*-----------Set target to zero
            Do j = jj, jj+njVec-1
               C(i,j) = Zero
            End Do
            Do k = 1, nColA
               If (A(i,k).ne.Zero) Then
                  Do j = jj, jj+njVec-1
                     C(i,j) = C(i,j) + A(i,k)*B(j,k)
                  End Do
               End If
            End Do
         End Do
*
      End Do    ! End of sectioning
*
      Return
      End
#else
c--------------------------------------------------------------------
      subroutine ntmul(a,b,r,ncol,nlink,nrow)
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (mxind=2000)
      dimension r(ncol,*),a(ncol,*),b(nrow,*),ind(mxind)
c
      do 100 i=1,ncol
c
      nnot=0
      do k=1,min(nlink,mxind)
        if ( a(i,k) .ne. 0.0d0 ) then
          nnot = nnot + 1
          ind(nnot) = k
        end if
      end do
c
      do 20 j=1,nrow-15,16
      s1=0.0D0
      s2=0.0D0
      s3=0.0D0
      s4=0.0D0
      s5=0.0D0
      s6=0.0D0
      s7=0.0D0
      s8=0.0D0
      t1=0.0D0
      t2=0.0D0
      t3=0.0D0
      t4=0.0D0
      t5=0.0D0
      t6=0.0D0
      t7=0.0D0
      t8=0.0D0
      do 15 k=1,nnot
      s1=s1+a(i,ind(k))*b(j,  ind(k))
      s2=s2+a(i,ind(k))*b(j+1,ind(k))
      s3=s3+a(i,ind(k))*b(j+2,ind(k))
      s4=s4+a(i,ind(k))*b(j+3,ind(k))
      s5=s5+a(i,ind(k))*b(j+4,ind(k))
      s6=s6+a(i,ind(k))*b(j+5,ind(k))
      s7=s7+a(i,ind(k))*b(j+6,ind(k))
      s8=s8+a(i,ind(k))*b(j+7,ind(k))
      t1=t1+a(i,ind(k))*b(j+8,ind(k))
      t2=t2+a(i,ind(k))*b(j+9,ind(k))
      t3=t3+a(i,ind(k))*b(j+10,ind(k))
      t4=t4+a(i,ind(k))*b(j+11,ind(k))
      t5=t5+a(i,ind(k))*b(j+12,ind(k))
      t6=t6+a(i,ind(k))*b(j+13,ind(k))
      t7=t7+a(i,ind(k))*b(j+14,ind(k))
      t8=t8+a(i,ind(k))*b(j+15,ind(k))
15    continue
      r(i,j  )=s1
      r(i,j+1)=s2
      r(i,j+2)=s3
      r(i,j+3)=s4
      r(i,j+4)=s5
      r(i,j+5)=s6
      r(i,j+6)=s7
      r(i,j+7)=s8
      r(i,j+8)=t1
      r(i,j+9)=t2
      r(i,j+10)=t3
      r(i,j+11)=t4
      r(i,j+12)=t5
      r(i,j+13)=t6
      r(i,j+14)=t7
      r(i,j+15)=t8
20    continue
c
      nr1=mod(nrow,16)
      if(nr1.eq.0) goto 100
      j=nrow-nr1+1
c
      if(nr1.ge.8) then
      s1=0.0D0
      s2=0.0D0
      s3=0.0D0
      s4=0.0D0
      s5=0.0D0
      s6=0.0D0
      s7=0.0D0
      s8=0.0D0
      do 25 k=1,nnot
      s1=s1+a(i,ind(k))*b(j,  ind(k))
      s2=s2+a(i,ind(k))*b(j+1,ind(k))
      s3=s3+a(i,ind(k))*b(j+2,ind(k))
      s4=s4+a(i,ind(k))*b(j+3,ind(k))
      s5=s5+a(i,ind(k))*b(j+4,ind(k))
      s6=s6+a(i,ind(k))*b(j+5,ind(k))
      s7=s7+a(i,ind(k))*b(j+6,ind(k))
      s8=s8+a(i,ind(k))*b(j+7,ind(k))
25    continue
      r(i,j  )=s1
      r(i,j+1)=s2
      r(i,j+2)=s3
      r(i,j+3)=s4
      r(i,j+4)=s5
      r(i,j+5)=s6
      r(i,j+6)=s7
      r(i,j+7)=s8
      nr1=nr1-8
      j=j+8
      end if
c
      if(nr1.ge.4) then
      s1=0.0D0
      s2=0.0D0
      s3=0.0D0
      s4=0.0D0
      do 35 k=1,nnot
      s1=s1+a(i,ind(k))*b(j,  ind(k))
      s2=s2+a(i,ind(k))*b(j+1,ind(k))
      s3=s3+a(i,ind(k))*b(j+2,ind(k))
      s4=s4+a(i,ind(k))*b(j+3,ind(k))
35    continue
      r(i,j  )=s1
      r(i,j+1)=s2
      r(i,j+2)=s3
      r(i,j+3)=s4
      nr1=nr1-4
      j=j+4
      end if
c
      If (nr1.eq.1) Then
         s1=0.0D0
         do 45 k=1,nnot
         s1=s1+a(i,ind(k))*b(j,ind(k))
45       continue
         r(i,j)=s1
      Else If (nr1.eq.2) Then
         s1=0.0D0
         s2=0.0D0
         do 50 k=1,nnot
         s1=s1+a(i,ind(k))*b(j,  ind(k))
         s2=s2+a(i,ind(k))*b(j+1,ind(k))
50       continue
         r(i,j  )=s1
         r(i,j+1)=s2
      Else If (nr1.eq.3) Then
         s1=0.0D0
         s2=0.0D0
         s3=0.0D0
         do 55 k=1,nnot
         s1=s1+a(i,ind(k))*b(j,  ind(k))
         s2=s2+a(i,ind(k))*b(j+1,ind(k))
         s3=s3+a(i,ind(k))*b(j+2,ind(k))
55       continue
         r(i,j  )=s1
         r(i,j+1)=s2
         r(i,j+2)=s3
      Else If (nr1.gt.3) Then
         Call WarningMessage(2,'nr1.gt.3')
         Call Abend()
      End If
 100  Continue
*
      if(mxind.ge.nlink) return
      Call WarningMessage(2,'MxInd.lt.nLink')
      Write (6,*) 'mxind,nlink=',mxind,nlink
      Call Abend()
      Return
      End
#endif
