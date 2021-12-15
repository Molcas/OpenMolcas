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
      Subroutine Mk_tVt(TInt,nTheta_All,tVt,nTheta,List2,mData,
     &                  iPrm,nPrm,iAng,jAng,nk,nl,Indkl,nkl,
     &                  iAL,nA,nB)
      Implicit Real*8 (a-h,o-z)
      Real*8 TInt(nTheta_All,nTheta_All), tVt(nTheta,nTheta)
      Integer List2(mData,nTheta_All), iPrm(nPrm), Indkl(nkl),
     &        iAL(nA,nB)
*
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('Mk_tVt: TInt',' ',TInt,nTheta_All,nTheta_All)
      Call iVcPrt('iPrm',' ',iPrm,nPrm)
      Call iVcPrt('Indkl',' ',Indkl,nkl)
      Call iVcPrt('iAL',' ',iAL,nA*nB)
#endif
      Call FZero(tVt,nTheta**2)
      iA=iAng+1
      jA=jAng+1
      Do iTheta_All = 1, nTheta_All
         kComp=List2(3,iTheta_All)
         lComp=List2(4,iTheta_All)
         ik=   List2(5,iTheta_All)
         il=   List2(6,iTheta_All)
         If (iAng.eq.jAng) Then
            iTheta_Full = ik*(ik-1)/2+il
         Else
            iTheta_Full = (il-1)*nk + ik
         End If
         If (
C    &       iAL(kComp,lComp).eq.1 .and.
     &       kComp.eq.iA .and. lComp.eq.jA .and.
     &       iPrm(iTheta_Full).eq.1) Then
            iTheta=Indkl(iTheta_Full)
*
            Do jTheta_All = 1, nTheta_All
               mComp=List2(3,jTheta_All)
               nComp=List2(4,jTheta_All)
               im   =List2(5,jTheta_All)
               in   =List2(6,jTheta_All)
               If (iAng.eq.jAng) Then
                  jTheta_Full = im*(im-1)/2+in
               Else
                  jTheta_Full = (in-1)*nk + im
               End If
               If (
C    &             iAL(mComp,nComp).eq.1 .and.
     &             mComp.eq.iA .and. nComp.eq.jA .and.
     &             iPrm(jTheta_Full).eq.1) Then
                  jTheta=Indkl(jTheta_Full)
*
                  tVt(iTheta,jTheta) = tVt(iTheta,jTheta)
     &                               +     TInt(iTheta_All,jTheta_All)
C    &                               + Abs(TInt(iTheta_All,jTheta_All))
*
               End If
            End Do
*
         End If
      End Do
*
#ifdef _DEBUGPRINT_
      Call RecPrt('tVt',' ',tVt,nTheta,nTheta)
#else
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(iAL)
#endif
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nl)
      End
