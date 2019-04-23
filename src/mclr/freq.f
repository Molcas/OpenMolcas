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
      Subroutine FREQ(nX,H,nDeg,nrvec,TmpAux,Tmp3,EVec,EVal,RedM,iNeg)
      Implicit Real*8 (a-h,o-z)
      Real*8 H(*), TmpAux(*),Tmp3(nX,nX),
     &       EVec(2*nX,nX),
     &       EVal(2*nX),
     &       RedM(nX)

#include "real.fh"
#include "constants2.fh"
#include "stdalloc.fh"
      Integer nrvec(*),ndeg(*)
      Logical Found
      Real*8, Allocatable :: Mass(:)
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

*
*----- read masses from runfile
*
      Call Qpg_dArray('Isotopes',Found,nIsot)
      If (.Not.Found) Then
        Write(6,*) 'No masses found on RunFile'
        Call AbEnd()
      End If
      Call mma_allocate(Mass,nIsot)
      Call Get_dArray('Isotopes',Mass,nIsot)
*
*----- form the Mass Weighted Cartesian force constant matrix.
*
      iprint=0
      Do i = 1, nX
         rm=Mass(nrvec(i))
         If (rm.eq.Zero) rm=1.0d7
         Do j=1,nX
            Tmp3(i,j) = sqrt(DBLE(nDeg(i)*nDeg(j)))*
     &                  H(itri(i,j))/rm
       End Do
      End Do
*
*-----Compute eigenvectors and eigenfunctions
*
      nAux = 2 * nX
      iOpt=1
      If ( nX.gt.0 ) then
        Call Not_DGeEv(iOpt,Tmp3,nX,
     &             EVal,EVec,nX,
     &             nX,TmpAux,nAux)
      End If
*
*-----Compute the harmonic frequencies
*
      iNeg=0
      Do 649 iHarm = 1, 2*nX, 2
         jHarm = (iHarm+1)/2
         temp = EVal(iHarm)
         If (temp.ge.Zero) Then
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
      Do iHarm = 1, nX
        r2=Zero
        Do i=1,nx
            rm=Mass(nrvec(i))/UTOAU
            r2=r2+Evec(2*(i-1)+1,iHarm)*Evec(2*(i-1)+1,iHarm)*rm
        End Do
        RedM(iHarm)=r2
        r2=One/Sqrt(r2)
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
               rlow=RedM(iHarm)
               RedM(iHarm)=RedM(jHarm)
               RedM(jHarm)=rLow
               Call DSwap_(nX,EVec(1,iHarm),2,EVec(1,jHarm),2)
            End If
         End Do
      End Do
*
      Call mma_deallocate(Mass)
*
      Return
      End
