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
      Subroutine Mk_TInt_P(TInt_p,nTInt_p,
     &                     TP,nTP,
     &                     iAL,nCompA,nCompB,
     &                     List2_p,nList2_p,
     &                     mData,iAng,jAng,npk,npl,
     &                     List_TP)
      Implicit Real*8 (a-h,o-z)
      Real*8 TInt_p(nTInt_p,nTInt_p), TP(nTP,nTP)
      Integer iAL(nCompA,nCompB), List2_p(mData,nList2_p),
     &        List_TP(2,nTP)
*
      iA=iAng+1
      jA=jAng+1
      Call FZero(TP,nTP**2)
      Do iList2_p = 1, nList2_p
         kAng= List2_p(1,iList2_p)
         lAng= List2_p(2,iList2_p)
         kComp=List2_p(3,iList2_p)
         lComp=List2_p(4,iList2_p)
C        Write (6,*) 'kComp,lComp=',kComp,lComp
         If (
     &       kAng.eq.iAng .and. lAng.eq.jAng
C    &       .and. iAL(kComp,lComp).eq.1
     &       .and. kComp.eq.iA .and. lComp.eq.jA
     &      ) Then
*
            k=List2_p(5,iList2_p)
            l=List2_p(6,iList2_p)
            If (iAng.eq.jAng) Then
               iTP = k*(k-1)/2 + l
            Else
               iTP = (l-1)*npk + k
            End If
            List_TP(1,iTP)=k
            List_TP(2,iTP)=l
*
            Do jList2_p = 1, nList2_p
               mAng= List2_p(1,jList2_p)
               nAng= List2_p(2,jList2_p)
               mComp=List2_p(3,jList2_p)
               nComp=List2_p(4,jList2_p)
C              Write (6,*) 'mComp,nComp=',mComp,nComp
               If (
     &             mAng.eq.iAng .and. nAng.eq.jAng
C    &             .and. iAL(mComp,nComp).eq.1
     &             .and. mComp.eq.iA .and. nComp.eq.jA
     &            ) Then
*
                  m=List2_p(5,jList2_p)
                  n=List2_p(6,jList2_p)
                  If (iAng.eq.jAng) Then
                     jTP = m*(m-1)/2 + n
                  Else
                     jTP = (n-1)*npk + m
                  End If
*
                  TP(iTP,jTP) = TP(iTP,jTP)
     &                        +     TInt_P(iList2_p,jList2_p)
C    &                        + Abs(TInt_P(iList2_p,jList2_p))
*
               End If
*
            End Do
*
         End If
      End Do
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iAL)
         Call Unused_integer(npl)
      End If
      End
