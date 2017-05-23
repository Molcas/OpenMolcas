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
      Subroutine BestMatch(nConstr,nOrb,Occ,Match,ldim)

      Implicit Real*8 (a-h,o-z)
      Integer nConstr, nOrb, ldim, Match(2,ldim)
      Real*8  Occ(nOrb)
*
      jConstr=1

10    Continue

      delta0=2.0d0
      Do i=1,nOrb
         Do j=1,i-1
            xOcc=Occ(i)+Occ(j)
            delta=abs(2.0d0-xOcc)
            If (delta.lt.delta0) Then
               delta0=delta
               If (Occ(i).gt.Occ(j)) Then
                  Match(1,jConstr)=i
                  Match(2,jConstr)=j
               Else
                  Match(1,jConstr)=j
                  Match(2,jConstr)=i
               EndIf
            EndIf
         End Do
      End Do
*
      If (jConstr.lt.nConstr) Then ! note: Occ array destroyed here
         k=Match(1,jConstr)
         Occ(k)=-42.0d0
         l=Match(2,jConstr)
         Occ(l)=-42.0d0
         jConstr=jConstr+1
         Go To 10
      EndIf
*
      Return
      End
