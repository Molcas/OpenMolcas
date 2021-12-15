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
      Subroutine RRINT(K,ALFA,A,BETA,R0,GRINT,lmax)
      Implicit Real*8(A-H,O-Z)
#include "real.fh"
#include "welcom.fh"
      Real*8 grint(0:lmax,lmax),rri(0:kmax+2)
*
      M=K+1
      EXPA=-A*A*ALFA-BETA*R0*R0
      AEXP=(ALFA+BETA)
      BEXP1 = -Two*(BETA*R0+A*ALFA)/AEXP
      BEXP2 = -Two*(BETA*R0-A*ALFA)/AEXP
      If (A.eq.Zero) Then
         Test =Zero
      Else
         If (R0.eq.Zero) Then
            Test=One
         Else
            TEST=A*ALFA
         End If
      End If
CFUE  IF(TEST.LT..02D+00)GO TO 900
      IF(TEST.LT..005D+00)GO TO 900
C     Write (*,*) ' Large A'
C.....K=0 ONE CONTRIBUTION SS-INTEGRAL
      Do 40 i=0,k
         rri(i)=qrint(i+1,aexp,bexp1,expa)*DBLE((-1)**i)
     &          -qrint(i+1,aexp,bexp2,expa)
 40   Continue
C     Call RecPrt(' In RRInt: rri',' ',rri,k+1,1)
      AL=One/(Two*ALFA*A)
      Do 41 i=0,k
         mmax=i/2
         Do 42 m=1,mmax+1
c.....calculate integral ri(i,m)
            fiintm=fiint(m-1,0)
            grint(i,m)=Zero
            Do 43 n=1,m
               Bi=binom(m-1,n-1)*DBLE((-1)**(n+1))
               ind=i-(m-n)*2
               Do 44 kk=0,ind
                  ggg=fac(ind)/fac(ind-kk)*al**(kk+1)
                  grint(i,m)=grint(i,m)+bi*ggg*rri(i-kk)*fiintm
 44            Continue
43          Continue
42       Continue
41    Continue
      go to 100
900   Continue
C     Write (*,*) ' SERIES EXPANSION FOR SMALL A'
c
C.....SERIES EXPANSION FOR SMALL A
c
C.... K=0 FIRST
      l = (k+1)/2
      EXP1=-(ALFA*A*A+BETA*R0*R0)
      BEXP=-Two*BETA*R0/(ALFA+BETA)
      Do 45 i=0,l+2
         rri(i)=qrint(2*(i+1),AExp,BExp,Exp1)
 45   Continue
C     Call RecPrt(' rri',' ',rri,l+3,1)
      pi4=pi*Four
      AA  = Two *(A*Alfa)
      AA2 = Two *(A*Alfa)**2
      AA3 = Four*(A*Alfa)**3
      AA4 = Two *(A*Alfa)**4
      AA5 = Four*(A*Alfa)**5
      GRINT(0,1) = pi4*(               rri(0)
     &           +      AA2/Three    * rri(1)
     &           +      AA4/15.0D+00 * rri(2)  )
      IF(K.EQ.0)  go to 100
      Do 20 ll=1,l
         Do kk = 1, ll+1
*
            tmp1= fiint(kk-1,0)/fiint(0,0)
            tmp2= -Two*DBLE(kk-1)
            tmp = Zero
            Do mm = 0, kk-1
               tmp3 = tmp1*binom(kk-1,mm)*(-One)**mm
               tmp4 = tmp2 + DBLE(mm)*Two + One
               tmp = tmp + tmp3 * (
     &              + One /(      DBLE(2*ll  )+tmp4)      * rri(ll)
     &              + AA2/(       DBLE(2*ll+2)+tmp4)      * rri(ll+1)
     &              + AA4/(Three*(DBLE(2*ll+4)+tmp4))     * rri(ll+2))
            End Do
            Grint(ll*2,  kk) = pi4 * tmp
*
            if (kk.eq.1) cycle
            tmp1= fiint(kk-2,0)/fiint(0,0)
            tmp=Zero
            Do mm = 1, kk-1
               tmp3 = tmp1*binom(kk-2,mm-1)*(-One)**(mm+1)
               tmp4 = tmp2 + DBLE(mm)*Two + One
               tmp = tmp - tmp3 * (
     &              + AA /       (DBLE(2*ll)  +tmp4)      * rri(ll)
     &              + AA3/(Three*(DBLE(2*ll+2)+tmp4))     * rri(ll+1)
     &              + AA5/(15.0D0*(DBLE(2*ll+4)+tmp4))    * rri(ll+2))
            End Do
            Grint(ll*2-1,kk-1) = pi4 * tmp
            Grint(ll*2-1,kk) = Zero
*
         End Do
 20   Continue
 100  Continue
*
*     Call RecPrt(' In RRint:grint',' ',grint,1+lmax,lmax)
      Return
      End
