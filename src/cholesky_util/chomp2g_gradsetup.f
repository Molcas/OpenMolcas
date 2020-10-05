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
* Copyright (C) 2010-2012, Jonas Bostrom                               *
*               2013, Roland Lindh                                     *
************************************************************************
      Subroutine ChoMP2g_GradSetup(irc,CMO)
************************************************************************
*                                                                      *
*     Author: Jonas Bostrom (2010-2012)                                *
*             Roland Lindh, Implementation of symmetry (some while in  *
*                           New Orleans during the 245th National ACS  *
*                           meeting, 7-11 April 2013).                 *
************************************************************************
#include "implicit.fh"
#include "WrkSpc.fh"
#include "chomp2g.fh"
#include "chomp2.fh"
#include "cholesky.fh"
#include "choorb.fh"
*
      Integer nOccAll(8), iOffCInv(8), iOffLRo(8,8), iOffLRb(8,8),
     &        nLRo(8), nLRb(8), iOffD(8), iAdrR(8), iOff_WJL(8),
     &        nB3(8), iOffB3(8,8), iOffCMOo(8), iAdrB(8)
      Real*8 CMO(*)
*
      Character*5 fname
      Character*6 fname2
*
      Character*9  ThisNm
      Character*17 SecNam
      Integer LuB(2), LuA(2)
      Parameter (SecNam = 'ChoMP2g_GradSetup', ThisNm = 'GradSetup')
*                                                                      *
************************************************************************
*                                                                      *
      iWr=1
      iRd=2
*                                                                      *
************************************************************************
*                                                                      *
*     Compute some sizes of some arrays
*
      nBas2    = 0
      lTriDens = 0
      lRecDens = 0
      nOrbBas = 0
      nOccBas = 0
      nVirBas = 0
      iAdr=1
      jAdr=1
      iOff_JL = 0
      iCMOo = 0
      Do iSym = 1, nSym
         iOffCMOo(iSym)=iCMOo
         iCMOo = iCMOo + nBas(iSym)*nOcc(iSym)
         iOff_WJL(iSym)=iOff_JL
         iOff_JL = iOff_JL + nMP2Vec(iSym)*NumCho(iSym)
         iAdrR(iSym)=iAdr
         iAdrB(iSym)=jAdr
         iOffD(iSym) = nBas2
         nBas2 = nBas2 + nBas(iSym)*nBas(iSym)
         lTriDens = lTriDens + nBas(iSym)*(nBas(iSym)+1)/2
         lRecDens = lRecDens + nBas(iSym)*nBas(iSym)
         nOccAll(iSym) = nFro(iSym) + nOcc(iSym)
         iOffCInv(iSym) = nOrbBas
         nOrbBas = nOrbBas + nOrb(iSym)*nBas(iSym)
         nOccBas = nOccBas + nOcc(iSym)*nBas(iSym)
         nVirBas = nVirBas + nVir(iSym)*nBas(iSym)
         nOO=0
         nBB=0
         iB3 = 0
         Do jSym = 1, nSym
            kSym = iEor(iSym-1,jSym-1)+1
*
            iOffLRo(iSym,jSym) = nOO
            iOffLRb(iSym,jSym) = nBB
            nOO = nOO + nOrb(jSym)*nOrb(kSym)
            nBB = nBB + nBas(jSym)*nBas(kSym)
            iOffB3(iSym,jSym)=iB3
            iB3 = iB3 + nBas(jSym)*nOcc(kSym)
         End Do
         nB3(iSym)=iB3
         nLRo(iSym)=nOO
         nLRb(iSym)=nBB
         iAdr = iAdr + nBB*nMP2Vec(iSym)
         jAdr = jAdr + nBB*NumCho(iSym)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('CMO_inv','Allo','Real',ip_CMO_inv,nOrbBas)
      Call FZero(Work(ip_CMO_inv),nOrbBas)
      Call GetMem('CMO_o  ','Allo','Real',ip_CMO_o  ,nOccBas)
      Call GetMem('CMO_v  ','Allo','Real',ip_CMO_v  ,nVirBas)
*                                                                      *
************************************************************************
*                                                                      *
*     Get memory for density matricies and overlap matrix
*
      Call GetMem('MP2Density','Allo','Real',ip_MP2Density,lRecDens)
      Call GetMem('SCFDensity','Allo','Real',ip_SCFDensity,lRecDens)
      Call GetMem('Overlap', 'Allo','Real', ip_Smat, lRecDens)
*
      Call GetMem('MP2TTotDensity','Allo','Real',ip_MP2TTotDensity,
     &            lTriDens)
      Call GetMem('MP2TDensity','Allo','Real',ip_MP2TDensity,lTriDens)
      Call GetMem('SCFTDensity','Allo','Real',ip_SCFTDensity,lTriDens)
      Call Getmem('OverlapT','Allo','Real', ip_STmat,lTriDens)
*                                                                      *
************************************************************************
*                                                                      *
*     Get the Variational MP2 Density and the HF density
*
      Call Get_D1ao_Var(Work(ip_MP2TTotDensity),lTriDens)
      Call Get_D1ao(Work(ip_SCFTDensity), lTriDens)
*                                                                      *
************************************************************************
*                                                                      *
*     Get the overlap matrix
*
      iSymlbl=1
      Call RdOne(irc,6,'Mltpl  0',1,Work(ip_STmat),iSymlbl)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the differential MP2 density
      Do i = 1, lTriDens
         Work(ip_MP2TDensity + i-1)=2.0d0*(Work(ip_MP2TTotDensity+i-1)-
     &                                Work(ip_SCFTDensity + i-1))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Modified for triangular storage to quadratic storage. Some of the
*     triangular storage is folded, i.e off diagonal element are
*     multiplied with 2.
*
      index = 0
      iOff = 0
      Do iSym = 1, nSym
         Do i = 1, nBas(iSym)
            Do j = 1, i-1
               Work(ip_MP2Density+ j-1 + nBas(iSym)*(i-1)+iOff) =
     &                            Work(ip_MP2TDensity+index)/2.0d0
               Work(ip_MP2Density+ i-1 + nBas(iSym)*(j-1)+iOff) =
     &                            Work(ip_MP2TDensity+index)/2.0d0
               Work(ip_SCFDensity+ j-1 + nBas(iSym)*(i-1)+iOff) =
     &                            Work(ip_SCFTDensity+index)/2.0d0
               Work(ip_SCFDensity+ i-1 + nBas(iSym)*(j-1)+iOff) =
     &                            Work(ip_SCFTDensity+index)/2.0d0
               Work(ip_Smat + j-1 + nBas(iSym)*(i-1)+iOff) =
     &                         Work(ip_STmat+index)
               Work(ip_Smat + i-1 + nBas(iSym)*(j-1)+iOff) =
     &                         Work(ip_STmat+index)
*
               index = index + 1
            End Do
            Work(ip_MP2Density+ i-1 + nBas(iSym)*(i-1)+iOff) =
     &                      Work(ip_MP2TDensity+index)
            Work(ip_SCFDensity+ i-1 + nBas(iSym)*(i-1)+iOff) =
     &                      Work(ip_SCFTDensity+index)
            Work(ip_Smat + i-1 + nBas(iSym)*(i-1)+iOff) =
     &                      Work(ip_STmat+index)
*
            index = index + 1
         End Do
         iOff = iOff + nBas(iSym)**2
      End Do
*
*     Deallocate temporary memory
*
      Call GetMem('MP2TTotDensity','Free','Real',ip_MP2TTotDensity,
     &            lTriDens)
      Call GetMem('MP2TDensity','Free','Real',ip_MP2TDensity,lTriDens)
      Call GetMem('SCFTDensity','Free','Real',ip_SCFTDensity,lTriDens)
      Call Getmem('OverlapT','Free','Real', ip_STmat,lTriDens)
*                                                                      *
************************************************************************
*                                                                      *
*     Calculate inverse CMO-matrix as C^-1 = C^T*S
*     --------------------------------------------
*
      iOff1 = 0
      iOff2 = 0
      Do iSym = 1, nSym
         Call dGemm_('T','N', nOrb(iSym),nBas(iSym),nBas(iSym),
     &              1.0d0,CMO(1+iOff2), nBas(iSym),
     &                    Work(ip_Smat+iOff1), nBas(iSym),
     &              0.0d0,Work(ip_CMO_inv+iOff2), nOrb(iSym))
*
         iOff1 =  iOff1 + nBas(iSym)**2
         iOff2 =  iOff2 + nOrb(iSym)*nBas(iSym)
      End Do
*
      Call GetMem('Overlap', 'Free','Real', ip_Smat, lRecDens)
*                                                                      *
************************************************************************
*                                                                      *
*     Calculate occupied and virtual inverse CMO-matrices
*     ---------------------------------------------------
*
      k_v = 1
      k_o = 1
      iOff2 = 0
      Do iSym = 1, nSym
         Do i = 1, nOrb(iSym)
            Do j = 1, nBas(iSym)
               If(i.gt. nFro(iSym) .and. (i .le. nOccAll(iSym))) Then
                  Work(ip_CMO_o+k_o-1) = CMO((i-1)*nBas(iSym)+j+iOff2)
                  k_o = k_o+1
               Else If(i.gt. nOccAll(iSym)) Then
                  Work(ip_CMO_v+k_v-1) = CMO((i-1)*nBas(iSym)+j+iOff2)
                  k_v = k_v+1
               End If
            End Do
         End Do
         iOff2 =  iOff2 + nOrb(iSym)*nBas(iSym)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Open Cholesky vector files.
*     ---------------------------
*
      iTypL = 1
      iTypR = 2
*
      Do iSym = 1, nSym
         Call ChoMP2_OpenF(1,iTypL,iSym)
         Call ChoMP2_OpenF(1,iTypR,iSym)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Open some temporaty files
*
      iSeed = 7
      LuXVec = IsFreeUnit(iSeed)
      FName='TMPV1'
      Call DaName_MF_WA(LuXVec,Fname)

      iSeed = 7
      LuYVec = isFreeUnit(iSeed)
      FName='TMPV2'
      Call DaName_MF_WA(LuYVec,Fname)

      iSeed = 7
      LuRVec = isFreeUnit(iSeed)
      FName='TMPV3'
      Call DaName_MF_WA(LuRVec,Fname)

      LuLVec = isFreeUnit(iSeed)
      FName='TMPV4'
      Call DaName_MF_WA(LuLVec,Fname)

      iSeed = 7
      LuBTmp = isFreeUnit(iSeed)
      FName='TMPV6'
      Call DaName_MF_WA(LuBTmp,Fname)
*                                                                      *
************************************************************************
*                                                                      *
*     Open files for the terms in Eqs. 35, 36, 40, and 41.
*     i=1 Coulombic and i=2 Exchange contributions.
*     This can be merged into singel 2-center and 3-center
*     contributions
*
      Do i = 1,2
*
         iSeed = 7
         LuA(i) = isFreeUnit(iSeed)
         Write(Fname2,'(A5,I1)') 'AMP2V',i
         Call DaName_MF_WA(LuA(i),Fname2)
*
         iSeed = 7
         LuB(i) = isFreeUnit(iSeed)
         Write(Fname2,'(A5,I1)') 'BMP2V',i
         Call DaName_MF_WA(LuB(i),Fname2)
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Max number of vectors in one batch. This is not foolproof, if you
*     run the code with too little ram you might get problems
*     ------------------------------------------------------------------
*
*     This has to be recoded to be dynamic. To be done after the
*     symmetry has been implemented in the code.
*
      maxvalue = 10
*                                                                      *
************************************************************************
*                                                                      *
      nVec = 0
      Do iSym = 1, nSym
         nVec = min(max(NumCho(iSym),nMP2Vec(iSym),nVec),maxvalue)
         If(nVec .lt. 1) Then
            Call ChoMP2_Quit(SecNam,'Insufficient memory','[1]')
         End If
      End Do
      iSym=1  ! Here we are!
*
*     Allocate memory for L-vectors
*     -----------------------------
*
      iVecFF = 1
      iVecOF = 2
      iVecVF = 3
      iVecFO = 4
      iVecOO = 5
      iVecVO = 6
      iVecFV = 7
      iVecOV = 8
      iVecVV = 9
*
      lCJK = nMoMo(iSym,iVecFF)*nVec
      Call GetMem('CJK','Allo','Real',ip_CJK, lCJK)
*
      lCKi = nMoMo(iSym,iVecOF)*nVec
      Call GetMem('CKi','Allo','Real',ip_CKi, lCKi)
*
      lCKa = nMoMo(iSym,iVecVF)*nVec
      Call GetMem('CKa','Allo','Real',ip_CKa, lCKa)
*
      lCiK = nMoMo(iSYm,iVecFO)*nVec
      Call GetMem('CiK','Allo','Real',ip_CiK, lCiK)
*
      lCij = nMoMo(iSym,iVecOO)*nVec
      Call GetMem('Cij','Allo','Real',ip_Cij, lCij)
*
      lCia = nMoMo(iSym,iVecVO)*nVec
      Call GetMem('Cia','Allo','Real',ip_Cia, lCia)
*
      lCaK = nMoMo(iSym,iVecFV)*nVec
      Call GetMem('CaK','Allo','Real',ip_CaK, lCaK)
*
      lCai = nMoMo(iSym,iVecOV)*nVec
      Call GetMem('Cai','Allo','Real',ip_Cai, lCai)
*
      lCab = nMoMo(iSym,iVecVV)*nVec
      Call GetMem('Cab','Allo','Real',ip_Cab, lCab)
*
      lCpq = nOrb(iSym)*nOrb(iSym)*nVec
      Call GetMem('Cpq','Allo','Real',ip_Cpq, lCpq)
*
*     The Cholesky vectors of the MP2 amplitudes,
*     see Eq. 26.
*
      lRia = nOcc(iSym)*nVir(iSym)*nVec
      Call GetMem('Ria','Allo','Real',ip_Ria, lRia)
*
      lCpn = nOrb(iSym)*nBas(iSym)*nVec
      Call GetMem('Cpn','Allo','Real',ip_Cpn, lCpn)
*
      lRin = nOcc(iSym)*nBas(iSym)*nVec
      Call GetMem('Rin','Allo','Real',ip_Rin, lRin)
*
      lCmn = nBas(iSym)*nBas(iSym)*nVec
      Call GetMem('Cmn','Allo','Real',ip_Cmn, lCmn)
*
      lRmn = nBas(iSym)*nBas(iSym)*nVec
      Call GetMem('Rmn','Allo','Real',ip_Rmn, lRmn)
*
      lUkn = nBas(iSym)*nBas(iSym)*nVec
      Call GetMem('Ukn','Allo','Real',ip_Ukn, lUkn)
*
      lVkn = nBas(iSym)*nBas(iSym)*nVec
      Call GetMem('Vkn','Allo','Real',ip_Vkn, lVkn)
*
      lWJL = NumCho(iSym)*nMP2Vec(iSym)
      Call GetMem('WJL','Allo','Real',ip_WJL, lWJL)
      Call FZero(Work(ip_WJL),lWJL)
*
      lWmjKJ = nVec*nVec*nOcc(iSym)*nBas(iSym)
      Call GetMem('WmjKJ','Allo','Real',ip_WmjKJ, lWmjKJ)
*
      lB3jl = nBas(iSym)*nOcc(iSym)*nVec
      Call GetMem('B3jl','Allo','Real',ip_B3jl, lB3jl)
*
      lB3kl = nBas(iSym)*nBas(iSym)*nVec
      Call GetMem('B3kl','Allo','Real',ip_B3kl, lB3kl)
*                                                                      *
************************************************************************
*                                                                      *
*     Start loop over batches of C_pq^J vectors
*
      Do iSym = 1, nSym
      Do iiJ = 1, NumCho(iSym), nVec
         nJ  = Min(nVec, NumCho(iSym) - (iiJ -1) )
*
*        Read C_pq^J-vectors from disk
*        -----------------------------

         lTot = nMoMo(iSym,iVecFF)*nJ
         iAdr = 1 + nMoMo(iSym,iVecFF)*(iiJ-1) + iAdrOff(iSym,iVecFF)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Work(ip_CJK),lTot,iAdr)

         lTot = nMoMo(iSym,iVecOF)*nJ
         iAdr = 1 + nMoMo(iSym,iVecOF)*(iiJ-1) + iAdrOff(iSym,iVecOF)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Work(ip_CKi),lTot,iAdr)

         lTot = nMoMo(iSym,iVecVF)*nJ
         iAdr = 1 + nMoMo(iSym,iVecVF)*(iiJ-1) + iAdrOff(iSym,iVecVF)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Work(ip_CKa),lTot,iAdr)

         lTot = nMoMo(iSym,iVecFO)*nJ
         iAdr = 1 + nMoMo(iSym,iVecFO)*(iiJ-1) + iAdrOff(iSym,iVecFO)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Work(ip_CiK),lTot,iAdr)

         lTot = nMoMo(iSym,iVecOO)*nJ
         iAdr = 1 + nMoMo(iSym,iVecOO)*(iiJ-1) + iAdrOff(iSym,iVecOO)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Work(ip_Cij),lTot,iAdr)

         lTot = nMoMo(iSym,iVecVO)*nJ
         iAdr = 1 + nMoMo(iSym,iVecVO)*(iiJ-1) + iAdrOff(iSym,iVecVO)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Work(ip_Cia),lTot,iAdr)

         lTot = nMoMo(iSym,iVecFV)*nJ
         iAdr = 1 + nMoMo(iSym,iVecFV)*(iiJ-1) + iAdrOff(iSym,iVecFV)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Work(ip_CaK),lTot,iAdr)

         lTot = nMoMo(iSym,iVecOV)*nJ
         iAdr = 1 + nMoMo(iSym,iVecOV)*(iiJ-1) + iAdrOff(iSym,iVecOV)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Work(ip_Cai),lTot,iAdr)

         lTot = nMoMo(iSym,iVecVV)*nJ
         iAdr = 1 + nMoMo(iSym,iVecVV)*(iiJ-1) + iAdrOff(iSym,iVecVV)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Work(ip_Cab),lTot,iAdr)

*        Put the C_pq^J-vectors together to one large vector
*        ------------------------------------------------

         iOff3 = 0
         Do iJ = 1, nJ
*
            iOffFF=0
            iOffOF=0
            iOffFO=0
            iOffOO=0
            iOffVF=0
            iOffFV=0
            iOffVO=0
            iOffOV=0
            iOffVV=0
            iOffBB=0
            Do jSym = 1, nSym
               kSym=iEor(iSym-1,jSym-1)+1
            Do i = 1, nFro(kSym)
               iOff1 = nMoMo(iSym,iVecFF)*(iJ-1) + (i-1)*nFro(jSym)
               Call dCopy_(nFro(jSym),Work(ip_CJK+iOff1+iOffFF),1,
     &                               Work(ip_Cpq+iOff3+iOffBB),1)
               iOff3 = iOff3+nFro(jSym)
               iOff1 = nMoMo(iSym,iVecOF)*(iJ-1) + (i-1)*nOcc(jSym)
               Call dCopy_(nOcc(jSym),Work(ip_CKi+iOff1+iOffOF),1,
     &                               Work(ip_Cpq+iOff3+iOffBB),1)
               iOff3 = iOff3+nOcc(jSym)
               iOff1 = nMoMo(iSym,iVecVF)*(iJ-1) + (i-1)*nVir(jSym)
               Call dCopy_(nVir(jSym),Work(ip_CKa+iOff1+iOffVF),1,
     &                               Work(ip_Cpq+iOff3+iOffBB),1)
               iOff3 = iOff3+nVir(jSym)
            End Do
*
            Do i = 1, nOcc(kSym)
               iOff1 = nMoMo(iSym,iVecFO)*(iJ-1) + (i-1)*nFro(jSym)
               Call dCopy_(nFro(jSym),Work(ip_CiK+iOff1+iOffFO),1,
     &                               Work(ip_Cpq+iOff3+iOffBB),1)
               iOff3 = iOff3 + nFro(jSym)
               iOff1 = nMoMo(iSym,iVecOO)*(iJ-1) + (i-1)*nOcc(jSym)
               Call dCopy_(nOcc(jSym),Work(ip_Cij+iOff1+iOffOO),1,
     &                               Work(ip_Cpq+iOff3+iOffBB),1)
               iOff3 = iOff3 + nOcc(jSym)
               iOff1 = nMoMo(iSym,iVecVO)*(iJ-1) + (i-1)*nVir(jSym)
               Call dCopy_(nVir(jSym),Work(ip_Cia+iOff1+iOffVO),1,
     &                               Work(ip_Cpq+iOff3+iOffBB),1)
               iOff3 = iOff3 + nVir(jSym)
            End Do
*
            Do i = 1, nVir(kSym)
               iOff1 = nMoMo(iSym,iVecFV)*(iJ-1) + (i-1)*nFro(jSym)
               Call dCopy_(nFro(jSym),Work(ip_CaK+iOff1+iOffFV),1,
     &                               Work(ip_Cpq+iOff3+iOffBB),1)
               iOff3 = iOff3 + nFro(jSym)
               iOff1 = nMoMo(iSym,iVecOV)*(iJ-1) + (i-1)*nOcc(jSym)
               Call dCopy_(nOcc(jSym),Work(ip_Cai+iOff1+iOffOV),1,
     &                               Work(ip_Cpq+iOff3+iOffBB),1)
               iOff3 = iOff3 + nOcc(jSym)
               iOff1 = nMoMo(iSym,iVecVV)*(iJ-1) + (i-1)*nVir(jSym)
               Call dCopy_(nVir(jSym),Work(ip_Cab+iOff1+iOffVV),1,
     &                               Work(ip_Cpq+iOff3+iOffBB),1)
               iOff3 = iOff3 + nVir(jSym)
            End Do
*
            iOffFF = iOffFF + nFro(jSym)*nFro(kSym)
            iOffOF = iOffOF + nOcc(jSym)*nFro(kSym)
            iOffFO = iOffFO + nFro(jSym)*nOcc(kSym)
            iOffOO = iOffFO + nOcc(jSym)*nOcc(kSym)
            iOffVF = iOffVF + nVir(jSym)*nFro(kSym)
            iOffFV = iOffFV + nFro(jSym)*nVir(kSym)
            iOffVO = iOffVO + nVir(jSym)*nOcc(kSym)
            iOffOV = iOffOV + nOcc(jSym)*nVir(kSym)
            iOffVV = iOffVV + nVir(jSym)*nVir(kSym)
            iOffBB = iOffBB + nBas(jSym)*nBas(kSym)
            End Do
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
            Do jSym = 1, nSym
               kSym = iEor(iSym-1,jSym-1)+1
*
*           Do backtransformations of the Cholesky-vectors to get C_mn.
*           ----------------------------------------------------------
*
*           (C^-1)^T x C_pq^J = C_nq^J
*
            iOff1 = nLRo(iSym)*(iJ-1) + iOffLRo(iSym,jSym)
            Call dGemm_('T', 'N',nBas(jSym),nOrb(kSym), nOrb(jSym),
     &                 1.0d0,Work(ip_CMO_inv+iOffCInv(jSym)),nOrb(jSym),
     &                       Work(ip_Cpq+iOff1), nOrb(jSym),
     &                 0.0d0,Work(ip_Cpn), nBas(jSym))
*
*           C_nq^J x  (C^-1) = C_nm^J
*
            iOff2 = nLRb(iSym)*(iJ-1) + iOffLRb(iSym,jSym)
            Call dGemm_('N','N',nBas(jSym),nBas(kSym), nOrb(kSym),
     &                 1.0d0,Work(ip_Cpn),nBas(jSym),
     &                       Work(ip_CMO_inv+iOffCInv(kSym)),nOrb(kSym),
     &                 0.0d0,Work(ip_Cmn+iOff2), nBas(jSym))
*
            End Do
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
            Do jSym = 1, nSym
               kSym = iEor(iSym-1,jSym-1)+1
*
*           C_nm^J x D(MP2)_mk = U_nk^J, Eq. (43)
*
            iOff = nLRb(iSym)*(iJ-1) + iOffLRb(iSym,jSym)
            iOff0= iOffD(kSym)
            Call dGemm_('N','N',nBas(jSym),nBas(kSym), nBas(kSym),
     &                 1.0d0,Work(ip_Cmn+iOff),nBas(jSym),
     &                       Work(ip_MP2Density+iOff0), nBas(kSym),
     &                 0.0d0,Work(ip_Ukn+iOff), nBas(jSym))
*
*           D(HF)_kn x C_nm^J = V_km^J, Eq. (42)
*
            iOff0= iOffD(jSym)
            Call dGemm_('N','N',nBas(jSym),nBas(kSym), nBas(jSym),
     &                 1.0d0,Work(ip_SCFDensity+iOff0), nBas(jSym),
     &                       Work(ip_Cmn+iOff), nBas(jSym),
     &                 0.0d0,Work(ip_Vkn+iOff), nBas(jSym))
*
            End Do
         End Do ! iJ
*                                                                      *
************************************************************************
*                                                                      *
*        Write C_nm^J, U_nk^J, and V_km^J to disk.
*
         lTot = nLRb(iSym)*nJ
*
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuLVec,iWr,Work(ip_Cmn),lTot,iAdr)
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuXVec,iWr,Work(ip_Ukn),lTot,iAdr)
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuYVec,iWr,Work(ip_Vkn),lTot,iAdr)
*                                                                      *
************************************************************************
*                                                                      *
         Do lSym = 1, nSym
         Do iiK = 1, nMP2Vec(lSym), nVec
            nK = Min(nVec, nMP2Vec(lSym) - (iiK-1))
*
*           Read the R-Vectors
*
            lTot = nMoMo(lSym,iVecOV)*nK
            iAdr = nMoMo(lSym,iVecOV)*(iiK-1) + 1
            Call dDaFile(lUnit_F(lSym,iTypR),iRd,Work(ip_Ria),lTot,iAdr)

*           Do backtransformations of the Amplitude-vectors to get Rmnx.
*           ------------------------------------------------------------
*           See Eq. (39).
*
            iOff_Ria = 0
            iOff_Rin = 0
            iOff_Rmn = 0
            Do iK = 1, nK

               Do jSym = 1, nSym
                  kSym = iEor(lSym-1,jSym-1)+1
*
*                 C_na x R_ai^K = R_ni^K  (Temporary storage)
*
               iOff1 = iOff_Ria
               iOff_Ria = iOff_Ria + nVir(jSym)*nOcc(kSym)
               iOff2 = iOff_Rin
               iOff_Rin = iOff_Rin + nBas(jSym)*nOcc(kSym)
               Call dGemm_('N','N',nBas(jSym),nOcc(kSym), nVir(jSym),
     &              1.0d0,Work(ip_CMO_v),nBas(jSym),
     &                    Work(ip_Ria+iOff1), nVir(jSym),
     &              0.0d0,Work(ip_Rin+iOff2), nBas(jSym))
*
*              R_ni^K x C_mi^T = R_nm^K
*
               iOff3 = iOff_Rmn
               iOff_Rmn = iOff_Rmn + nBas(jSym)*nBas(kSym)
               Call dGemm_('N','T',nBas(jSym),nBas(kSym), nOcc(kSym),
     &                 1.0d0,Work(ip_Rin+iOff2),nBas(jSym),
     &                       Work(ip_CMO_o), nBas(kSym),
     &                 0.0d0,Work(ip_Rmn+iOff3), nBas(jSym))
*
               End Do ! jSym
            End Do    ! iK
*
*           Write  R_mn^K to disk
*
            lTot = nLRb(iSym)*nK
            iAdr = iAdrR(iSym) +  nLRb(iSym)*(iiK-1)
            Call dDaFile(LuRVec,iWr,Work(ip_Rmn),lTot,iAdr)
*                                                                      *
************************************************************************
*                                                                      *
*           Sum_mn (R_mn^K)^T x C_mn^J = W_KJ  (Eq. 38)
*
            If (iSym.eq.lSym) Then
            iOffZ = iOff_WJL(iSym) + nMP2Vec(iSym)*(iiJ-1) + iiK-1
            Call dGemm_('T','N',nK,nJ,nLRb(iSym),
     &                 1.0d0,Work(ip_Rmn),nLRb(iSym),
     &                       Work(ip_Cmn),nLRb(iSym),
     &                 0.0d0,Work(ip_WJL+iOffZ), nMP2Vec(iSym))
            End If
*                                                                      *
************************************************************************
*                                                                      *
            Do iJ = 1, nJ
               iOff_Rin=0
               Do iK = 1, nK
*
                  Do jSym = 1, nSym
                     kSym = iEor(iSym-1,jSym-1)+1
                     mSym = iEor(lSym-1,kSym-1)+1
*
                     iOffL = nLRb(iSym)*(iJ-1) + iOffLRb(iSym,jSym)
                     iOffR = iOff_Rin
                     iOff_Rin = iOff_Rin + nBas(kSym)*nOcc(mSym)
*
*                 C_mn^J x R_ni^K = W_mi^JK, see Eq. (44)
*
                  Call dGemm_('N','N',nBas(jSym),nOcc(mSym),nBas(kSym),
     &                        1.0d0,Work(ip_Cmn+iOffL),nBas(jSym),
     &                              Work(ip_Rin+iOffR),nBas(kSym),
     &                        0.0d0,Work(ip_WmjKJ),nBas(jSym))
*
                  Fac = 1.0d0
                  If((iK .eq. 1) .and. (iiK .eq. 1)) Fac = 0.0d0
*
*                 R_nm^K x W_mi^JK, last two summation in Eqs. 40 and 41
*
                  mSym2 = iEor(lSym-1,jSym-1)+1
                  iOffB = nB3(iSym)*(iJ-1) + iOffB3(iSym,mSym2)
                  iOffR = nLRb(lSym)*(iK-1) + iOffLRb(lSym,mSym2)
                  Call dGemm_('N','N',nBas(mSym2),nOcc(mSym),nBas(jSym),
     &                       1.0d0,Work(ip_Rmn+iOffR),nBas(mSym2),
     &                             Work(ip_WmjKJ), nBas(jSym),
     &                       Fac,  Work(ip_B3jl+iOffB),nBas(mSym2))
*
                  End Do ! jSym
*
               End Do ! iK
            End Do    ! iJ
*
         End Do ! iiK
         End Do ! lSym

         Do iJ = 1, nJ
*
            Do jSym = 1, nSym
               kSym = iEor(iSym-1,jSym-1)+1
*
            iOffB1 = nB3(iSym)*(iJ-1) + iOffB3(iSym,jSym)
            iOffB2 = nLRb(iSym)*(iJ-1) + iOffLRb(iSym,jSym)
            iOff   = iOffCMOo(kSym)
*
*           Complete the 3rd RHS term in Eq. 40
*
            Call dGemm_('N','T',nBas(jSym),nBas(kSym),nOcc(kSym),
     &                -8.0d0,Work(ip_B3jl+iOffB1) , nBas(jSym),
     &                       Work(ip_CMO_o+iOff),nBas(kSym),
     &                 0.0d0,Work(ip_B3kl+iOffB2),nBas(jSym))
*
            End Do ! jSym
         End Do    ! iJ
*
*        Write matrix to file. Note that it has the same structure as
*        the R and L matrices.
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuBTmp,iWr,Work(ip_B3kl),lTot,iAdr)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End Do ! iiJ
      End Do ! iSym
*
      Call GetMem('B3jl','Free','Real',ip_B3jl, lB3jl)
*
      Call GetMem('Rin','Free','Real',ip_Rin, lRin)
      Call GetMem('Ria','Free','Real',ip_Ria, lRia)
*
      Call GetMem('Cpq','Free','Real',ip_Cpq, lCpq)
      Call GetMem('Cpn','Free','Real',ip_Cpn, lCpn)
*
      Call GetMem('CJK','Free','Real',ip_CJK, lCJK)
      Call GetMem('CKi','Free','Real',ip_CKi, lCKi)
      Call GetMem('CKa','Free','Real',ip_CKa, lCKa)
      Call GetMem('CiK','Free','Real',ip_CiK, lCiK)
      Call GetMem('Cij','Free','Real',ip_Cij, lCij)
      Call GetMem('Cia','Free','Real',ip_Cia, lCia)
      Call GetMem('CaK','Free','Real',ip_CaK, lCaK)
      Call GetMem('Cai','Free','Real',ip_Cai, lCai)
      Call GetMem('Cab','Free','Real',ip_Cab, lCab)
      Call GetMem('SCFDensity','Free','Real',ip_SCFDensity,lRecDens)
*
      iSym = 1
*
      lB1kl = nBas(iSym)*nBas(iSym)*nVec
      Call GetMem('B1kl','Allo','Real',ip_B1kl, lB1kl)
*
      lB3kl_s = nBas(iSym)*nBas(iSym)*nVec
      Call GetMem('B3kl_s','Allo','Real',ip_B3kl_s, lB3kl_s)
*
      lA1 = NumCho(iSym)*NumCho(iSym)
      Call GetMem('A1','Allo','Real',ip_A1, lA1)
      Call FZero(Work(ip_A1),lA1)
*
      lA2 = NumCho(iSym)*NumCho(iSym)
      Call GetMem('A2','Allo','Real',ip_A2, lA2)
      Call FZero(Work(ip_A2),lA2)
*
      lB2kl = nBas(iSym)*nBas(iSym)*nVec
      Call GetMem('B2kl','Allo','Real',ip_B2kl, lB2kl)
*
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Start loop over batches of J
*
      iOff_A = 0
      Do iSym = 1, nSym
      Do iiJ = 1, NumCho(iSym), nVec
         nJ = Min(nVec, NumCho(iSym) - (iiJ-1))
*
*        Read V_km^J
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuYVec,iRd,Work(ip_Vkn),lTot,iAdr)
*
         Do iiK = 1, NumCho(iSym), nVec
            nK = Min(nVec, NumCho(iSym) - (iiK-1))
*
*           Read U_km_K
*
            lTot = nLRb(iSym)*nK
            iAdr = iAdrB(iSym) + nLRb(iSym)*(iiK-1)
            Call dDaFile(LuXVec,iRd,Work(ip_Ukn),lTot,iAdr)
*
*           U_km^K x V_km^J, 1st term RHS eq. 41

            iOffA = iiK-1 + (iiJ-1)*NumCho(iSym) + iOff_A
            Fact = 1.0D0
            Call dGemm_('T','N', nK, nJ,nLRb(iSym),
     &                 Fact,Work(ip_Ukn),nLRb(iSym),
     &                       Work(ip_Vkn),nLRb(iSym),
     &                 0.0d0,Work(ip_A1+iOffA),NumCho(iSym))
*
         End Do
*                                                                      *
************************************************************************
*                                                                      *
         Do iiK = 1, nMP2Vec(iSym), nVec
            nK = Min(nVec, nMP2Vec(iSym) - (iiK-1))
*
*           Read R_mn^K
*
            lTot = nLRb(iSym)*nK
            iAdr = iAdrR(iSym) + nLRb(iSym)*(iiK-1)
            Call dDaFile(LuRVec,iRd,Work(ip_Rmn),lTot,iAdr)
*
*           R_mn_K x W_JK, 2nd RHS term Eq. 35
*
            Fac = 1.0D0
            If (iiK.eq.1) Fac = 0.0D0
            iOffZ = iOff_WJL(iSym) + iiK-1 + (iiJ-1)*nMp2Vec(iSym)
            Call dGemm_('N','N',nLRb(iSym),nJ,nK,
     &                 -8.0d0,Work(ip_Rmn), nLRb(iSym),
     &                        Work(ip_WJL+iOffZ),nMP2Vec(iSym),
     &                  Fac,  Work(ip_B2kl),nLRb(iSym))
         End Do
*
*        Write to disk
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuB(2),iWr,Work(ip_B2kl),lTot,iAdr)
*                                                                      *
************************************************************************
*                                                                      *
*        Read 3rd RHS term in Eq. 40
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuBTmp,iRd,Work(ip_B3kl),lTot,iAdr)
*
*        Symmetrize
*
         Do iJ = 1, nJ
            iOffL = nLRb(iSym)*(iJ-1)
            Do jSym = 1, nSym
               kSym= iEor(iSym-1,jSym-1)+1
*
            Do k = 1, nBas(jSym)
               Do l = 1, nBas(kSym)
                  kl = k + nBas(jSym)*(l-1) + iOffLRb(iSym,jSym)
                  kl_s = kl
                  lk = l + nBas(kSym)*(k-1) + iOffLRb(iSym,kSym)
                  Work(ip_B3kl_s+iOffL+kl_s-1) =
     &               ( Work(ip_B3kl+iOffL+kl-1) +
     &                 Work(ip_B3kl+iOffL+lk-1) )/2
               End Do
            End Do
*
            End Do ! jSym
         End Do
*
*        Write the symmetrized 3rd term to disk
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuBTmp,iWr,Work(ip_B3kl_s),lTot,iAdr)
*                                                                      *
************************************************************************
*                                                                      *
*
         Do iJ = 1, nJ
            iOff = nLRb(iSym)*(iJ-1)
*
*           V_kn^J x D(MP2)_kn, 2nd RHS term Eq. 40
*
            Do jSym = 1, nSym
               kSym = iEor(iSym-1,jSym-1) + 1
*
               iOff1 = iOff + iOffLRb(iSym,jSym)
               iOff2 = iOffD(kSym)
            Call dGemm_('N','N',nBas(jSym), nBas(kSym), nBas(kSym),
     &                 1.0d0,Work(ip_Vkn+iOff1), nBas(jSym),
     &                       Work(ip_MP2Density+iOff2), nBas(kSym),
     &                 0.0d0,Work(ip_B1kl+iOff1), nBas(jSym))
            End Do
         End Do
*
*        Compound 2nd and 3rd RHS term in Eq. 40.

         Call DaXpY_(nLRb(iSym)*nJ,1.0d0,Work(ip_B3kl),1,
     &                                     Work(ip_B1kl),1)
*
*        Write compounded terms to disk
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuB(1),iWr,Work(ip_B1kl),lTot,iAdr)
*                                                                      *
************************************************************************
*                                                                      *
      End Do ! iiJ
         iOff_A = iOff_A + NumCho(iSym)**2
      End Do ! iSym
*
      Call GetMem('MP2Density','Free','Real',ip_MP2Density,lRecDens)
      iSym = 1
*                                                                      *
************************************************************************
*                                                                      *
*     Add nonseparable coulombic terms to A matrix
*     --------------------------------------------
*
*     -4 W_MJ x W_MK    (don't ask me why the code say 8!)
*
*     Second term RHS Eq. 36, note the different order of the indecies.
*
      iOff_A2 = 0
      Do iSym = 1, nSym
         iOff1 = iOff_WJL(iSym)
         Fact=-8.0D0
      Call dGemm_('T','N', NumCho(iSym),NumCho(iSym),nMP2Vec(iSym),
     &            Fact ,Work(ip_WJL+iOff1), nMP2Vec(iSym),!-2
     &                  Work(ip_WJL+iOff1), nMP2Vec(iSym),
     &            0.0d0,Work(ip_A2+iOff_A2),NumCho(iSym))
*
         iOff_A2 = iOff_A2 + NumCho(iSym)**2
      End Do
*
      iSym = 1
*                                                                      *
************************************************************************
*                                                                      *
*     Start loop over batches of L
*
***** Nonseparable exchange, third loop, transformation and
*
      iOff_A = 0
      Do iSym = 1, nSym
      Do iiJ = 1, NumCho(iSym), nVec
         nJ = Min(nVec,NumCho(iSym)-(iiJ-1))
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nlRb(iSym)*(iiJ-1)
         Call dDaFile(LuBTmp,iRd,Work(ip_B3kl_s),lTot,iAdr)
*
         Do iiK = 1, NumCho(iSym), nVec
            nK = Min(nVec,NumCho(iSym)-(iiK-1))
*
            lTot = nLRb(iSym)*nK
            iAdr = iAdrB(iSym) + nLRb(iSym)*(iiK-1)
            Call dDaFile(LuLVec,iRd,Work(ip_Cmn),lTot,iAdr)
*
*           B_nm^L x C_nm^K
*
*           This is the second term of the RHS of Eq. (41)
*
            iOffA = iiJ-1 + (iiK-1)*NumCho(iSym) + iOff_A
            Fact = 1.0D0
            Call dGemm_('T', 'N', nJ, nK, nLRb(iSym),
     &                 Fact ,Work(ip_B3kl_s),nLRb(iSym),
     &                       Work(ip_Cmn), nLRb(iSym),
     &                 1.0d0,Work(ip_A1+iOffA),NumCho(iSym))
*
         End Do
*
      End Do
         iOff_A = iOff_A + NumCho(iSym)**2
      End Do ! iSym
*                                                                      *
************************************************************************
*                                                                      *
*     Write the Coulomic and Exchange contributions to the 2-center
*     integrals to disk. These terms could be merged to just one.
*
      iAdr1 = 1
      iAdr2 = 1
      iOff  = 0
      Do iSym = 1, nSym
         lTot = NumCho(iSym)*NumCho(iSym)
         Call dDaFile(LuA(1),iWr,Work(ip_A1+iOff),lTot,iAdr1)
         Call dDaFile(LuA(2),iWr,Work(ip_A2+iOff),lTot,iAdr2)
         iOff = iOff + lTot
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Close some files
*
      Call DaClos(LuXVec)
      Call DaClos(LuYVec)
      Call DaClos(LuRVec)
      Call DaClos(LuLVec)
      Call DaClos(LuBTmp)
      Do i = 1,2
         Call DaClos(LuA(i))
         Call DaClos(LuB(i))
      End Do
*
      iClos = 2
      Do iSym = 1, nSym
         Call ChoMP2_OpenF(iClos,iTypR,iSym)
         Call ChoMP2_OpenF(iClos,iTypL,iSym)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocate memory
*
      Call GetMem('Cmn','Free','Real',ip_Cmn, lCmn)
*
      Call GetMem('Rmn','Free','Real',ip_Rmn, lRmn)
*
      Call GetMem('Ukn','Free','Real',ip_Ukn, lUkn)
      Call GetMem('Vkn','Free','Real',ip_Vkn, lVkn)
*
      Call GetMem('B1kl','Free','Real',ip_B1kl, lB2kl)
      Call GetMem('A1','Free','Real',ip_A1, lA1)
      Call GetMem('B2kl','Free','Real',ip_B2kl, lB2kl)
      Call GetMem('A2','Free','Real',ip_A2, lA2)
      Call GetMem('B3kl','Free','Real',ip_B3kl, lB3kl)
      Call GetMem('B3kl_s','Free','Real',ip_B3kl_s, lB3kl_s)
      Call GetMem('WmjKJ','Free','Real',ip_WmjKJ, lWmjKJ)
      Call GetMem('WJL','Free','Real',ip_WJL, lWJL)
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('CMO_v  ','Free','Real',ip_CMO_v  ,nVirBas)
      Call GetMem('CMO_o  ','Free','Real',ip_CMO_o  ,nOccBas)
      Call GetMem('CMO_inv','Free','Real',ip_CMO_inv,nOrbBas)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
