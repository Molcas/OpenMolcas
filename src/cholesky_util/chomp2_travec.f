************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_TraVec(VecAO,VecMO,COcc,CVir,Scr,lScr,
     &                         iSyCho,iSyCO,iSyCV,iLoc)
C
C     Thomas Bondo Pedersen, Dec. 2004.
C
C     Purpose: compute ai-vector from reduced set AO vector.
C
#include "implicit.fh"
      Real*8 VecAO(*), VecMO(*), COcc(*), CVir(*)
      Real*8 Scr(lScr)
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Character*13 SecNam
      Parameter (SecNam = 'ChoMP2_TraVec')

      Real*8 Fac(0:1)
      Data Fac /0.5D0,1.0D0/

      iRS2F(i,j)=iWork(ip_iRS2F-1+2*(j-1)+i)
      IndRed(i,j)=iWork(ip_IndRed-1+nnBstrT(1)*(j-1)+i)
      MulD2h(i,j)=iEor(i-1,j-1)+1

      If (iLoc.lt.2 .or. iLoc.gt.3) Then
         Write(6,*) SecNam,': illegal iLoc = ',iLoc
         Call ChoMP2_Quit(SecNam,'iLoc out of bounds!',' ')
      End If
      iSyScr = MulD2h(iSyCho,iSyCO)
      If (lScr .lt. nT1AOT(iSyScr)) Then
         Write(6,*) SecNam,': insufficient scratch space lScr = ',lScr
         Write(6,*) SecNam,': needed                          = ',
     &              nT1AOT(iSyScr)
         Call ChoMP2_Quit(SecNam,'Insufficient scratch space',' ')
      Else
         Call Cho_dZero(Scr,nT1AOT(iSyScr))
      End If

C     First half-transformation step:
C     Scr(i,alpha) = sum_beta VecAO(alpha,beta)*COcc(i,beta)
C     ------------------------------------------------------

      If (iSyCho .eq. 1) Then

         Do iAlBe = 1,nnBstR(iSyCho,iLoc)

            jAlBe  = IndRed(iiBstR(iSyCho,iLoc)+iAlBe,iLoc)
            iAlpha = iRS2F(1,jAlBe)
            iBeta  = iRS2F(2,jAlBe)

            iSymAl = 1
            Do iSym = nSym,2,-1
               If (iAlpha .gt. iBas(iSym)) Then
                  iSymAl = iSym
                  Go To 998
               End If
            End Do
  998       iSymBe = iSymAl
            iSymi  = MulD2h(iSymBe,iSyCO)

            jAlpha = iAlpha - iBas(iSymAl)
            jBeta  = iBeta  - iBas(iSymBe)

            AOVal  = Fac(min(abs(iAlpha-iBeta),1))*VecAO(iAlBe)
            kOffAl = iT1AOT(iSymi,iSymAl) + nOcc(iSymi)*(jAlpha - 1)
            kOffBe = iT1AOT(iSymi,iSymBe) + nOcc(iSymi)*(jBeta  - 1)
            Do i = 1,nOcc(iSymi)
               Scr(kOffAl+i) = Scr(kOffAl+i) + AOVal*COcc(kOffBe+i)
               Scr(kOffBe+i) = Scr(kOffBe+i) + AOVal*COcc(kOffAl+i)
            End Do

         End Do

      Else

         Do iAlBe = 1,nnBstR(iSyCho,iLoc)

            jAlBe  = IndRed(iiBstR(iSyCho,iLoc)+iAlBe,iLoc)
            iAlpha = iRS2F(1,jAlBe)
            iBeta  = iRS2F(2,jAlBe)

            iSymAl = 1
            Do iSym = nSym,2,-1
               If (iAlpha .gt. iBas(iSym)) Then
                  iSymAl = iSym
                  Go To 999
               End If
            End Do
  999       iSymBe = MulD2h(iSymAl,iSyCho)

            jAlpha = iAlpha - iBas(iSymAl)
            jBeta  = iBeta  - iBas(iSymBe)

            AOVal  = VecAO(iAlBe)

            iSymi  = MulD2h(iSymBe,iSyCO)
            kOffAl = iT1AOT(iSymi,iSymAl) + nOcc(iSymi)*(jAlpha - 1)
            kOffBe = iT1AOT(iSymi,iSymBe) + nOcc(iSymi)*(jBeta  - 1)
            Do i = 1,nOcc(iSymi)
               Scr(kOffAl+i) = Scr(kOffAl+i) + AOVal*COcc(kOffBe+i)
            End Do

            iSymi  = MulD2h(iSymAl,iSyCO)
            kOffAl = iT1AOT(iSymi,iSymAl) + nOcc(iSymi)*(jAlpha - 1)
            kOffBe = iT1AOT(iSymi,iSymBe) + nOcc(iSymi)*(jBeta  - 1)
            Do i = 1,nOcc(iSymi)
               Scr(kOffBe+i) = Scr(kOffBe+i) + AOVal*COcc(kOffAl+i)
            End Do

         End Do

      End If

C     Second half-transformation step:
C     VecMO(a,i) = sum_alpha CVir(alpha,a)*Scr(i,alpha)
C     -------------------------------------------------

      Do iSymi = 1,nSym

         iSyma  = MulD2h(iSymi,iSyCho)
         iSymAl = MulD2h(iSyma,iSyCV)

         nTotAl = nBas(iSymAl)
         nTota  = nVir(iSyma)
         nToti  = nOcc(iSymi)

         If (nToti.gt.0 .and. nTota.gt.0 .and. nTotAl.gt.0) Then
            kOff1 = iAOVir(iSymAl,iSyma) + 1
            kOff2 = iT1AOT(iSymi,iSymAl) + 1
            kOff3 = iT1am(iSyma,iSymi)   + 1
            Call DGEMM_('T','T',nVir(iSyma),nOcc(iSymi),nBas(iSymAl),
     &                 1.0D0,CVir(kOff1),nTotAl,Scr(kOff2),nToti,
     &                 0.0D0,VecMO(kOff3),nTota)
         End If

      End Do

      End
