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
       Subroutine Fix_Exponents(nP,mP,nC,Exp,CoeffC,CoeffP)
       Implicit Real*8 (a-h,o-z)
       Real*8 Exp(nP), CoeffC(nP*nC*2), CoeffP(nP*nP*2)
*
       mP = nP
*
       Call Fix_Exp(nP,mP,nC,Exp,CoeffC,CoeffP)
*
       If (mP.ne.nP) Then
          nCP=nC*nP
          Do iC = 2, nC
             Do iP = 1, mP
                iFrom = (iC-1)*nP + iP
                iTo   = (iC-1)*mP + iP
                CoeffC(iTo    )=CoeffC(iFrom    )
                CoeffC(iTo+nCP)=CoeffC(iFrom+nCP)
             End Do
          End Do
          mCP=nC*mP
          Do iCP = 1, mCP
             iFrom = iCP + nCP
             iTo   = iCP + mCP
             CoeffC(iTo) = CoeffC(iFrom)
          End Do
*
          nPP=nP*nP
          Do iC = 2, mP
             Do iP = 1, mP
                iFrom = (iC-1)*nP + iP
                iTo   = (iC-1)*mP + iP
                CoeffP(iTo    )=CoeffP(iFrom    )
                CoeffP(iTo+nCP)=CoeffP(iFrom+nCP)
             End Do
          End Do
          mCP=mP*mP
          Do iCP = 1, mCP
             iFrom = iCP + nCP
             iTo   = iCP + mCP
             CoeffP(iTo) = CoeffP(iFrom)
          End Do
       End If
*
       Return
       End
       Subroutine Fix_Exp(nP,mP,nC,Exp,CoeffC,CoeffP)
       Implicit Real*8 (a-h,o-z)
       Real*8 Exp(nP), CoeffC(nP,nC,2), CoeffP(nP,nP,2)
*
*      First, put the exponents with all zero coefficients
*      at the end.
*
       Do iP = nP, 1, -1
*
          iSkip=1
          Do iC = 1, nC
             If (Abs(CoeffC(iP,iC,1)).ne.0.0D0) iSkip=0
          End Do
*
          If (iSkip.eq.1) Then
             If (iP.lt.mP) Then
                Temp   =Exp(iP)
                Exp(iP)=Exp(mP)
                Exp(mP)=Temp
                Temp           =CoeffP(iP,iP,1)
                CoeffP(iP,iP,1)=CoeffP(mP,mP,1)
                CoeffP(mP,mP,1)=Temp
                Temp           =CoeffP(iP,iP,2)
                CoeffP(iP,iP,1)=CoeffP(mP,mP,2)
                CoeffP(mP,mP,1)=Temp
                Do iC = 1, nC
                   Temp            = CoeffC(iP,iC,1)
                   CoeffC(iP,iC,1) = CoeffC(mP,iC,1)
                   CoeffC(mP,iC,1) = Temp
                   Temp            = CoeffC(iP,iC,2)
                   CoeffC(iP,iC,2) = CoeffC(mP,iC,2)
                   CoeffC(mP,iC,2) = Temp
                End Do
             End If
             mP = mP -1
          End If
*
       End Do
*
*      Second, order from largest to smallest
*
       Do iP = 1, mP-1
          Do jP = iP+1, mP
             If (Exp(jP).gt.Exp(ip)) Then
                Temp   =Exp(iP)
                Exp(iP)=Exp(jP)
                Exp(jP)=Temp
                Temp           =CoeffP(iP,iP,1)
                CoeffP(iP,iP,1)=CoeffP(jP,jP,1)
                CoeffP(jP,jP,1)=Temp
                Temp           =CoeffP(iP,iP,2)
                CoeffP(iP,iP,1)=CoeffP(jP,jP,2)
                CoeffP(jP,jP,1)=Temp
                Do iC = 1, nC
                   Temp            = CoeffC(iP,iC,1)
                   CoeffC(iP,iC,1) = CoeffC(jP,iC,1)
                   CoeffC(jP,iC,1) = Temp
                   Temp            = CoeffC(iP,iC,2)
                   CoeffC(iP,iC,2) = CoeffC(jP,iC,2)
                   CoeffC(jP,iC,2) = Temp
                End Do
             End If
          End Do
       End Do
*
       Return
       End
