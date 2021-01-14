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
* Copyright (C) 2010, Jonas Bostrom                                    *
************************************************************************
      SubRoutine ChoMP2g_TraVec(VecAO,VecMO,COrb1,COrb2,Scr,lScr,
     &                         iSyCho,iSyCO,iSyCV,iLoc,
     &                         iMoType1,iMoType2)
C
C     Jonas Bostrom, Feb 2010
C
C     Purpose: compute pq-vector from reduced set AO vector.
C
      use ChoArr, only: iRS2F
#include "implicit.fh"
      Real*8 VecAO(*), VecMO(*), COrb1(*), COrb2(*)
      Real*8 Scr(lScr)
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2.fh"
#include "chomp2g.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Character*13 SecNam
      Parameter (SecNam = 'ChoMP2_TraVec')

      Real*8 Fac(0:1)
      Data Fac /0.5D0,1.0D0/

      IndRed(i,j)=iWork(ip_IndRed-1+nnBstrT(1)*(j-1)+i)
      MulD2h(i,j)=iEor(i-1,j-1)+1

*     Check what type of Cholesky vector to make (fro-occ, occ-occ.....)
      iVecType = iMoType2 + (iMoType1-1)*nMoType

      If (iLoc.lt.2 .or. iLoc.gt.3) Then
         Write(6,*) SecNam,': illegal iLoc = ',iLoc
         Call ChoMP2_Quit(SecNam,'iLoc out of bounds!',' ')
      End If
      iSyScr = MulD2h(iSyCho,iSyCO)
      If (lScr .lt. nMoAo(iSyScr,iMoType1)) Then
         Write(6,*) SecNam,': insufficient scratch space lScr = ',lScr
         Write(6,*) SecNam,': needed                          = ',
     &              nMoAo(iSyScr,iMoType1)
         Call ChoMP2_Quit(SecNam,'Insufficient scratch space',' ')
      Else
         Call Cho_dZero(Scr,nMoAo(iSyScr,iMoType1))
      End If

C     First half-transformation step:
C     Scr(i,alpha) = sum_beta VecAO(alpha,beta)*COrb1(i,beta)
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
            iSymP  = MulD2h(iSymBe,iSyCO)

            jAlpha = iAlpha - iBas(iSymAl)
            jBeta  = iBeta  - iBas(iSymBe)

            AOVal  = Fac(min(abs(iAlpha-iBeta),1))*VecAO(iAlBe)
            kOffAl = iMoAo(iSymP,iSymAl,iMoType1) +  nMo(iSymP,iMoType1)
     &                                       * (jAlpha - 1)
            kOffBe = iMoAo(iSymP,iSymBe,iMoType1) + nMo(iSymP,iMoType1)
     &                                       * (jBeta  - 1)
            Do i = 1,nMo(iSymP,iMoType1)
               Scr(kOffAl+i) = Scr(kOffAl+i) + AOVal*COrb1(kOffBe+i)
               Scr(kOffBe+i) = Scr(kOffBe+i) + AOVal*COrb1(kOffAl+i)
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

            iSymP  = MulD2h(iSymBe,iSyCO)
            kOffAl = iMoAo(iSymP,iSymAl,iMoType1)
     &             + nMo(iSymP,iMoType1)*(jAlpha - 1)
            kOffBe = iMoAo(iSymP,iSymBe,iMoType1)
     &             + nMo(iSymP,iMoType1)*(jBeta  - 1)
            Do i = 1,nMo(iSymP,iMoType1)
               Scr(kOffAl+i) = Scr(kOffAl+i) + AOVal*COrb1(kOffBe+i)
            End Do

            iSymP  = MulD2h(iSymAl,iSyCO)
            kOffAl = iMoAo(iSymP,iSymAl,iMoType1)
     &             + nMo(iSymP,iMoType1)*(jAlpha - 1)
            kOffBe = iMoAo(iSymP,iSymBe,iMoType1)
     &             + nMo(iSymP,iMoType1)*(jBeta  - 1)
            Do i = 1,nMo(iSymP,iMoType1)
               Scr(kOffBe+i) = Scr(kOffBe+i) + AOVal*COrb1(kOffAl+i)
            End Do

         End Do

      End If

C     Second half-transformation step:
C     VecMO(q,p) = sum_alpha COrb2(alpha,q)*Scr(p,alpha)
C     -------------------------------------------------

      Do iSymp = 1,nSym

         iSymq  = MulD2h(iSymp,iSyCho)
         iSymAl = MulD2h(iSymq,iSyCV)

         nTotAl = nBas(iSymAl)
         nTotq  = nMo(iSymq,iMoType2)
         nTotp  = nMo(iSymp,iMoType1)

         If (nTotp.gt.0 .and. nTotq.gt.0 .and. nTotAl.gt.0) Then
            kOff1 = iAoMo(iSymAl,iSymq,iMoType2) + 1
            kOff2 = iMoAo(iSymp,iSymAl,iMoType1) + 1
            kOff3 = iMoMo(iSymq,iSymp,iVecType)   + 1
            Call DGEMM_('T','T',nMo(iSymq,iMoType2),nMo(iSymp,iMoType1),
     &                 nBas(iSymAl),1.0D0,COrb2(kOff1),nTotAl,
     &                 Scr(kOff2),nTotp,0.0D0,VecMO(kOff3),nTotq)
         End If

      End Do

      End
