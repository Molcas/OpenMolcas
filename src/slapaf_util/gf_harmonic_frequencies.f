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
      Subroutine GF_Harmonic_Frequencies(G,GInv,Tmp1,Tmp2,
     &                                   EVec,EVal,RedM,iNeg,nX,mInter)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "constants2.fh"
      Real*8 G(nX**2),GInv(nX**2), Tmp1(mInter,mInter), Tmp2(nX**2),
     &       EVec(2*mInter,mInter),  EVal(2*mInter), RedM(mInter)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute eigenvectors and eigenfunctions
*
      Call FZero(Tmp1,mInter**2)
      call dcopy_(nX,[One],0,Tmp1,mInter+1)
      Call NIDiag_new(Tmp2,Tmp1,mInter,mInter,0)
*
      Call FZero(EVal,2*mInter)
      Call FZero(EVec,2*mInter**2)
*
*     Move over eigenvalue and eigenvectors, note that the eigenvectors
*     are transformed back to Cartesian coordinates from mass-weighted
*     Cartesians.
*
      Do iX = 1, mInter
         iiT= iX*(iX+1)/2
         EVal((iX-1)*2+1)=Tmp2(iiT)
         call dcopy_(mInter,Tmp1(1,iX),1,EVec(1,iX),2)
*
         r2 = Zero
         Do jX = 1, mInter
            jj = (jX-1)*mInter + jX
            tmp=EVec((jX-1)*2+1,iX)
            tmp=tmp*Sqrt(G(jj))
            EVec((jX-1)*2+1,iX)=tmp
            r2 = r2 + tmp**2
         End Do
         Call DScal_(mInter,One/Sqrt(r2),EVec(1,iX),2)
      End Do
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('EVal',' ',EVal,2,mInter)
      Call RecPrt('EVec',' ',EVec,mInter*2,mInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the harmonic frequencies
*                                                                      *
************************************************************************
*                                                                      *
      iNeg=0
      Do iHarm = 1, 2*mInter, 2
         jHarm = (iHarm+1)/2
         temp = EVal(iHarm)
*
*        Fix imaginary frequencies
*
         If (temp.ge.Zero) Then
            EVal(jHarm) = Sqrt(temp)*autocm
         Else
            iNeg=iNeg+1
            EVal(jHarm) = -Sqrt(Abs(temp))*autocm
         End If
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('Converted EVal',' ',EVal,1,mInter)
#endif
*
*-----Normalize over the metric of the masses
*
      Do iHarm = 1, mInter
         call dcopy_(mInter,EVec(1,iHarm),2,Tmp1,1)
         Call DGEMM_('N','N',
     &               mInter,1,mInter,
     &               1.0d0,GInv,mInter,
     &               Tmp1,mInter,
     &               0.0d0,Tmp2,mInter)
         r2=DDot_(mInter,Tmp1,1,Tmp2,1)
         RedM(iHarm)=r2
         r2=One/Sqrt(r2)
         Call DScal_(mInter,r2,EVec(1,iHarm),2)
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('Normal coordinates (Q)',' ', EVec,mInter*2,mInter)
#endif
*
*-----Order, from low to high. Put translations and rotations last.
*
      Do iHarm = 1, mInter-1
         Test_i=EVal(iHarm)
         If (Abs(Test_i).lt.1.0D-3) Test_i = 1.0D+5
         Do jHarm = iHarm+1, mInter
            Test_j=EVal(jHarm)
            If (Abs(Test_j).lt.1.0D-3) Test_j = 1.0D+5
            If (Test_j.lt.Test_i) Then
               rlow=Test_i
               Test_i=Test_j
               Test_j=rlow
               rlow=EVal(iHarm)
               EVal(iHarm)=EVal(jHarm)
               EVal(jHarm)=rLow
               rlow=RedM(iHarm)
               RedM(iHarm)=RedM(jHarm)
               RedM(jHarm)=rLow
               Call DSwap_(mInter,EVec(1,iHarm),2,EVec(1,jHarm),2)
            End If
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('Frequencies (cm-1)',' ',EVal,1,mInter)
      Call RecPrt('Reduced masses (u)',' ',RedM,1,mInter)
      Call RecPrt('Normal Coordinates',' ',EVec,mInter*2,mInter)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
