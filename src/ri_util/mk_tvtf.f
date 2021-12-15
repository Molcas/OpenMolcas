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
      Subroutine Mk_tVtF(TInt,nTheta_All,tVtF,nTheta,List2,mData,
     &                   iPrm,nPrm,
     &                  iAng,jAng,nk,nl,Indkl,nkl,
     &                  nTheta_Full,
     &                   iAL,nA,nB)
      Implicit Real*8 (a-h,o-z)
      Real*8 TInt(nTheta_All,nTheta_All), tVtF(nTheta,nTheta_Full)
      Integer List2(mData,nTheta_All), Indkl(nkl), iPrm(nPrm),
     &        iAL(nA,nB)
*
#ifdef _DEBUGPRINT_
      Call RecPrt('TInt',' ',TInt,nTheta_all,nTheta_all)
      Call iVcPrt('Indkl',' ',Indkl,nkl)
#endif
      Call FZero(tVtF,nTheta*nTheta_Full)
      iA=iAng+1
      jA=jAng+1
      Do iTheta_All = 1, nTheta_All
         kComp= List2(3,iTheta_All)
         lComp= List2(4,iTheta_All)
         ik=    List2(5,iTheta_All)
         il=    List2(6,iTheta_All)
         If (iAng.eq.jAng) Then
            iTheta_Full = ik*(ik-1)/2+il
         Else
            iTheta_Full = (il-1)*nk + ik
         End If
         If (
     &       iPrm(iTheta_Full).eq.1
C    &       .and. iAL(kComp,lComp).eq.1
     &       .and. kComp.eq.iA .and. lComp.eq.jA
     &      ) Then
            iTheta=Indkl(iTheta_Full)
*
            Do jTheta_All = 1, nTheta_All
               mComp= List2(3,jTheta_All)
               nComp= List2(4,jTheta_All)
C              If (iAL(mComp,nComp).eq.1) Then
               If (mComp.eq.iA.and.nComp.eq.jA) Then
                  jk=List2(5,jTheta_All)
                  jl=List2(6,jTheta_All)
                  If (iAng.eq.jAng) Then
                     jTheta_Full = jk*(jk-1)/2+jl
                  Else
                     jTheta_Full = (jl-1)*nk + jk
                  End If
*
                  tVtF(iTheta,jTheta_Full)=tVtF(iTheta,jTheta_Full)
C    &                                 +Abs(TInt(iTheta_All,jTheta_All))
     &                                 +    TInt(iTheta_All,jTheta_All)
*
               End If
            End Do
*
         End If
      End Do
#ifdef _DEBUGPRINT_
      Call RecPrt('tVtF',' ',tVtF,nTheta,nTheta_Full)
#endif
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nl)
         Call Unused_integer_array(iAL)
      End If
      End
