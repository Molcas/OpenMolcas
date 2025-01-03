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
      SubRoutine ABXpY(Array1,Array2,idsym)
      use MCLR_Data, only: ipMO, NA
      use input_mclr, only: nSym,nAsh,nIsh,nOrb
      Implicit None
      Integer idsym
      Real*8 Array1(*),Array2(*)

      Integer i, j, itri
      Integer iS, jS, kS, lS, ijS
      Integer iA, jA, kA, lA
      Integer iAsh, jAsh, kAsh, lAsh
      Integer iiA, ij, kl, ijkl, ip1
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      Do iS=1,nSym
       Do iA=1,Nash(is)
        iAsh=nA(is)+iA
        iiA=nIsh(is)+iA
        Do jS=1,nSym
         ijs=iEOR(is-1,js-1)+1
         Do jA=1,Nash(js)
          jAsh=nA(js)+jA
          ij=itri(iash,jash)
          If (iAsh.ge.jash) Then
          Do kS=1,nSym
           Do kA=1,Nash(ks)
            kAsh=nA(ks)+kA
             ls=ieor(iEOR(kS-1,ijs-1),idsym-1)+1
             Do lA=1,Nash(ls)
              lAsh=nA(ls)+lA
              If (kAsh.ge.lash) Then
              kl=itri(kAsh,lash)
              If (ij.ge.kl) Then
              ijkl=itri(ij,kl)
              ip1=ipMO(js,ks,ls)+
     &             nOrb(is)*nAsh(js)*((lA-1)*nAsh(kS)+kA-1)+
     &             nOrb(is)*(ja-1)+iia-1
              Array2(ijkl)=array1(ip1)+array2(ijkl)
              End If
              End If
            End Do
           End Do
          End Do
          End If
         End Do
        End Do
       End Do
      End Do
      End SubRoutine ABXpY
