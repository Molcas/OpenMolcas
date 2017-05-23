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
      Subroutine FREQ(nX,H,nDeg,nrvec,TmpAux,Tmp3,EVec,EVal,iNeg)
      Implicit Real*8 (a-h,o-z)
      Real*8 H(*), TmpAux(*),Tmp3(nX,nX),
     &       EVec(2*nX,nX),
     &       EVal(2*nX)

#include "constants2.fh"
      Integer nrvec(*),ndeg(*)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

*
*----- form the Mass Weighted Cartesian force constant matrix.
*
      iprint=0
      Do i = 1, nX
       Do j=1,nX
            rm=rmass(nrvec(i))
            If (rm.eq.0.0d0) rm=1.0d7
            Tmp3(i,j) = sqrt(DBLE(nDeg(i)*nDeg(j)))*
     &                  H(itri(i,j))/rm
       End Do
      End Do
*
*-----Compute eigenvectors and eigenfunctions
*
      nAux = 2 * nX
      iOpt=1
      islct=0
      If ( nX.gt.0 ) then
        Call Not_DGeEv(iOpt,Tmp3,nX,
     &             EVal,EVec,nX,
     &             iSlct,nX,TmpAux,nAux)
      End If
*
*-----Compute the harmonic frequencies
*
      iNeg=0
      Do 649 iHarm = 1, 2*nX, 2
         jHarm = (iHarm+1)/2
         temp = EVal(iHarm)
         If (temp.ge.0.0d0) Then
            EVal(jHarm) = Sqrt(temp)*autocm
         Else
            iNeg=iNeg+1
            EVal(jHarm) = -Sqrt(Abs(temp))*autocm
         End If
 649  Continue
      If (iPrint.ge.99) Call RecPrt('Converted EVal',' ',EVal,1,nX)
      If (iPrint.ge.99) Call RecPrt('Normalized EVec',' ',
     &                               EVec,nX*2,nX)
*
*-----Normalize
*
      Do iHarm=  1, nX
        r2=0.0d0
        Do i=1,nx
            rm=rmass(nrvec(i))/UTOAU
            r2=r2+Evec(2*(i-1)+1,iHarm)*Evec(2*(i-1)+1,iHarm)*rm
        End Do
        r2=1.0d0/Sqrt(r2)
        Call DScal_(nX,r2,EVec(1,iHarm),2)
      End Do
*
*-----Order, from low to high.
*
      Do iHarm = 1, nX-1
         Do jHarm = iHarm+1, nX
            If (EVal(jHarm).lt.EVal(iHarm)) Then
               rlow=EVal(iHarm)
               EVal(iHarm)=EVal(jHarm)
               EVal(jHarm)=rLow
               Call DSwap_(nX,EVec(1,iHarm),2,EVec(1,jHarm),2)
            End If
         End Do
      End Do
*
      Return
      End
