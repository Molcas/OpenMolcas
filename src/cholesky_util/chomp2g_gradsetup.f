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
#include "real.fh"
#include "chomp2g.fh"
#include "chomp2.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"
*
      Integer nOccAll(8), iOffCInv(8), iOffLRo(8,8), iOffLRb(8,8),
     &        nLRo(8), nLRb(8), iOffD(8), iAdrR(8), iOff_WJL(8),
     &        nB3(8), iOffB3(8,8), iOffCMOo(8), iAdrB(8)
      Real*8 CMO(*)
*
      Character(LEN=5) fname
      Character(LEN=6) fname2
*
      Character(LEN=9), Parameter:: ThisNm = 'GradSetup'
      Character(LEN=17), Parameter:: SecNam = 'ChoMP2g_GradSetup'

      Integer LuB(2), LuA(2)

      Real*8, Allocatable:: CMO_Inv(:), CMO_o(:), CMO_v(:)
      Real*8, Allocatable:: MP2Density(:), SCFDensity(:), SMat(:)
      Real*8, Allocatable:: MP2TTotDensity(:)
      Real*8, Allocatable:: MP2TDensity(:)
      Real*8, Allocatable:: SCFTDensity(:)
      Real*8, Allocatable:: STMat(:)

      Real*8, Allocatable:: CJK(:), CKi(:), CKa(:), CiK(:), Cij(:)
      Real*8, Allocatable:: Cia(:), CaK(:), Cai(:), Cab(:), Cpq(:)

      Real*8, Allocatable:: Ria(:), Cpn(:), Rin(:), Cmn(:), Rmn(:)
      Real*8, Allocatable:: Ukn(:), Vkn(:), WJL(:), WmjKJ(:)
      Real*8, Allocatable:: B3jl(:), B3kl(:), B3kl_s(:)
      Real*8, Allocatable:: B1kl(:), B2kl(:)
      Real*8, Allocatable:: A1(:), A2(:)
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
      Call mma_allocate(CMO_inv,nOrbBas,Label='CMO_Inv')
      CMO_inv(:)=0.0D0
      Call mma_allocate(CMO_o,nOccBas,Label='CMO_o')
      Call mma_allocate(CMO_v,nVirBas,Label='CMO_v')
*                                                                      *
************************************************************************
*                                                                      *
*     Get memory for density matricies and overlap matrix
*
      Call mma_allocate(MP2Density,lRecDens,Label='MP2Density')
      Call mma_allocate(SCFDensity,lRecDens,Label='SCFDensity')
      Call mma_allocate(SMat,lRecDens,Label='SMat')
*
      Call mma_allocate(MP2TTotDensity,lTriDens,Label='MP2TTotDensity')
      Call mma_allocate(MP2TDensity,lTriDens,Label='MP2TDensity')
      Call mma_allocate(SCFTDensity,lTriDens,Label='SCFTDensity')
      Call mma_allocate(STMat,lTriDens,Label='STMat')
*                                                                      *
************************************************************************
*                                                                      *
*     Get the Variational MP2 Density and the HF density
*
      Call Get_D1ao_Var(MP2TTotDensity,lTriDens)
      Call Get_dArray_chk('D1ao',SCFTDensity,lTriDens)
*                                                                      *
************************************************************************
*                                                                      *
*     Get the overlap matrix
*
      iSymlbl=1
      Call RdOne(irc,6,'Mltpl  0',1,STmat,iSymlbl)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the differential MP2 density
      MP2TDensity(:)=2.0d0*(MP2TTotDensity(:)- SCFTDensity(:))
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
               index = index + 1
               MP2Density(j + nBas(iSym)*(i-1)+iOff) =
     &                            MP2TDensity(index)/2.0d0
               MP2Density(i + nBas(iSym)*(j-1)+iOff) =
     &                            MP2TDensity(index)/2.0d0
               SCFDensity(j + nBas(iSym)*(i-1)+iOff) =
     &                            SCFTDensity(index)/2.0d0
               SCFDensity(i + nBas(iSym)*(j-1)+iOff) =
     &                            SCFTDensity(index)/2.0d0
               Smat(j + nBas(iSym)*(i-1)+iOff) = STmat(index)
               Smat(i + nBas(iSym)*(j-1)+iOff) = STmat(index)
*
            End Do
            index = index + 1
            MP2Density(i + nBas(iSym)*(i-1)+iOff) = MP2TDensity(index)
            SCFDensity(i + nBas(iSym)*(i-1)+iOff) = SCFTDensity(index)
            Smat(i + nBas(iSym)*(i-1)+iOff) = STmat(index)
*
         End Do
         iOff = iOff + nBas(iSym)**2
      End Do
*
*     Deallocate temporary memory
*
      Call mma_deallocate(MP2TTotDensity)
      Call mma_deallocate(MP2TDensity)
      Call mma_deallocate(SCFTDensity)
      Call mma_deallocate(STMat)
*                                                                      *
************************************************************************
*                                                                      *
*     Calculate inverse CMO-matrix as C^-1 = C^T*S
*     --------------------------------------------
*
      iOff1 = 1
      iOff2 = 1
      Do iSym = 1, nSym
         Call dGemm_('T','N', nOrb(iSym),nBas(iSym),nBas(iSym),
     &              1.0d0,CMO(iOff2), nBas(iSym),
     &                    Smat(iOff1), nBas(iSym),
     &              0.0d0,CMO_inv(iOff2), nOrb(iSym))
*
         iOff1 =  iOff1 + nBas(iSym)**2
         iOff2 =  iOff2 + nOrb(iSym)*nBas(iSym)
      End Do
*
      Call mma_deallocate(SMat)
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
                  CMO_o(k_o) = CMO((i-1)*nBas(iSym)+j+iOff2)
                  k_o = k_o+1
               Else If(i.gt. nOccAll(iSym)) Then
                  CMO_v(k_v) = CMO((i-1)*nBas(iSym)+j+iOff2)
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
      Call mma_allocate(CJK,lCJK,Label='CJK')
*
      lCKi = nMoMo(iSym,iVecOF)*nVec
      Call mma_allocate(CKi,lCKi,Label='CKi')
*
      lCKa = nMoMo(iSym,iVecVF)*nVec
      Call mma_allocate(CKa,lCKa,Label='CKa')
*
      lCiK = nMoMo(iSYm,iVecFO)*nVec
      Call mma_allocate(CiK,lCiK,Label='CiK')
*
      lCij = nMoMo(iSym,iVecOO)*nVec
      Call mma_allocate(Cij,lCij,Label='Cij')
*
      lCia = nMoMo(iSym,iVecVO)*nVec
      Call mma_allocate(Cia,lCia,Label='Cia')
*
      lCaK = nMoMo(iSym,iVecFV)*nVec
      Call mma_allocate(CaK,lCaK,Label='CaK')
*
      lCai = nMoMo(iSym,iVecOV)*nVec
      Call mma_allocate(Cai,lCai,Label='Cai')
*
      lCab = nMoMo(iSym,iVecVV)*nVec
      Call mma_allocate(Cab,lCab,Label='Cab')
*
      lCpq = nOrb(iSym)*nOrb(iSym)*nVec
      Call mma_allocate(Cpq,lCpq,Label='Cpq')
*
*     The Cholesky vectors of the MP2 amplitudes,
*     see Eq. 26.
*
      lRia = nOcc(iSym)*nVir(iSym)*nVec
      Call mma_allocate(Ria,lRia,Label='Ria')
*
      lCpn = nOrb(iSym)*nBas(iSym)*nVec
      Call mma_allocate(Cpn,lCpn,Label='Cpn')
*
      lRin = nOcc(iSym)*nBas(iSym)*nVec
      Call mma_allocate(Rin,lRin,Label='Rin')
*
      lCmn = nBas(iSym)*nBas(iSym)*nVec
      Call mma_allocate(Cmn,lCmn,Label='Cmn')
*
      lRmn = nBas(iSym)*nBas(iSym)*nVec
      Call mma_allocate(Rmn,lRmn,Label='Rmn')
*
      lUkn = nBas(iSym)*nBas(iSym)*nVec
      Call mma_allocate(Ukn,lUkn,Label='Ukn')
*
      lVkn = nBas(iSym)*nBas(iSym)*nVec
      Call mma_allocate(Vkn,lVkn,Label='Vkn')
*
      lWJL = NumCho(iSym)*nMP2Vec(iSym)
      Call mma_allocate(WJL,lWJL,Label='WJL')
      WJL(:)=Zero
*
      lWmjKJ = nVec*nVec*nOcc(iSym)*nBas(iSym)
      Call mma_allocate(WmjKJ,lWmjKJ,Label='WmjKJ')
*
      lB3jl = nBas(iSym)*nOcc(iSym)*nVec
      Call mma_allocate(B3jl,lB3jl,Label='B3jl')
*
      lB3kl = nBas(iSym)*nBas(iSym)*nVec
      Call mma_allocate(B3kl,lB3kl,Label='B3kl')
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
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,CJK,lTot,iAdr)

         lTot = nMoMo(iSym,iVecOF)*nJ
         iAdr = 1 + nMoMo(iSym,iVecOF)*(iiJ-1) + iAdrOff(iSym,iVecOF)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,CKi,lTot,iAdr)

         lTot = nMoMo(iSym,iVecVF)*nJ
         iAdr = 1 + nMoMo(iSym,iVecVF)*(iiJ-1) + iAdrOff(iSym,iVecVF)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,CKa,lTot,iAdr)

         lTot = nMoMo(iSym,iVecFO)*nJ
         iAdr = 1 + nMoMo(iSym,iVecFO)*(iiJ-1) + iAdrOff(iSym,iVecFO)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,CiK,lTot,iAdr)

         lTot = nMoMo(iSym,iVecOO)*nJ
         iAdr = 1 + nMoMo(iSym,iVecOO)*(iiJ-1) + iAdrOff(iSym,iVecOO)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Cij,lTot,iAdr)

         lTot = nMoMo(iSym,iVecVO)*nJ
         iAdr = 1 + nMoMo(iSym,iVecVO)*(iiJ-1) + iAdrOff(iSym,iVecVO)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Cia,lTot,iAdr)

         lTot = nMoMo(iSym,iVecFV)*nJ
         iAdr = 1 + nMoMo(iSym,iVecFV)*(iiJ-1) + iAdrOff(iSym,iVecFV)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,CaK,lTot,iAdr)

         lTot = nMoMo(iSym,iVecOV)*nJ
         iAdr = 1 + nMoMo(iSym,iVecOV)*(iiJ-1) + iAdrOff(iSym,iVecOV)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Cai,lTot,iAdr)

         lTot = nMoMo(iSym,iVecVV)*nJ
         iAdr = 1 + nMoMo(iSym,iVecVV)*(iiJ-1) + iAdrOff(iSym,iVecVV)
         Call dDaFile(lUnit_F(iSym,iTypL),iRd,Cab,lTot,iAdr)

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
               Call dCopy_(nFro(jSym),CJK(1+iOff1+iOffFF),1,
     &                               Cpq(1+iOff3+iOffBB),1)
               iOff3 = iOff3+nFro(jSym)
               iOff1 = nMoMo(iSym,iVecOF)*(iJ-1) + (i-1)*nOcc(jSym)
               Call dCopy_(nOcc(jSym),Cki(1+iOff1+iOffOF),1,
     &                               Cpq(1+iOff3+iOffBB),1)
               iOff3 = iOff3+nOcc(jSym)
               iOff1 = nMoMo(iSym,iVecVF)*(iJ-1) + (i-1)*nVir(jSym)
               Call dCopy_(nVir(jSym),Cka(1+iOff1+iOffVF),1,
     &                               Cpq(1+iOff3+iOffBB),1)
               iOff3 = iOff3+nVir(jSym)
            End Do
*
            Do i = 1, nOcc(kSym)
               iOff1 = nMoMo(iSym,iVecFO)*(iJ-1) + (i-1)*nFro(jSym)
               Call dCopy_(nFro(jSym),CiK(1+iOff1+iOffFO),1,
     &                               Cpq(1+iOff3+iOffBB),1)
               iOff3 = iOff3 + nFro(jSym)
               iOff1 = nMoMo(iSym,iVecOO)*(iJ-1) + (i-1)*nOcc(jSym)
               Call dCopy_(nOcc(jSym),Cij(1+iOff1+iOffOO),1,
     &                               Cpq(1+iOff3+iOffBB),1)
               iOff3 = iOff3 + nOcc(jSym)
               iOff1 = nMoMo(iSym,iVecVO)*(iJ-1) + (i-1)*nVir(jSym)
               Call dCopy_(nVir(jSym),Cia(1+iOff1+iOffVO),1,
     &                               Cpq(1+iOff3+iOffBB),1)
               iOff3 = iOff3 + nVir(jSym)
            End Do
*
            Do i = 1, nVir(kSym)
               iOff1 = nMoMo(iSym,iVecFV)*(iJ-1) + (i-1)*nFro(jSym)
               Call dCopy_(nFro(jSym),CaK(1+iOff1+iOffFV),1,
     &                               Cpq(1+iOff3+iOffBB),1)
               iOff3 = iOff3 + nFro(jSym)
               iOff1 = nMoMo(iSym,iVecOV)*(iJ-1) + (i-1)*nOcc(jSym)
               Call dCopy_(nOcc(jSym),Cai(1+iOff1+iOffOV),1,
     &                               Cpq(1+iOff3+iOffBB),1)
               iOff3 = iOff3 + nOcc(jSym)
               iOff1 = nMoMo(iSym,iVecVV)*(iJ-1) + (i-1)*nVir(jSym)
               Call dCopy_(nVir(jSym),Cab(1+iOff1+iOffVV),1,
     &                               Cpq(1+iOff3+iOffBB),1)
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
     &                 1.0d0,CMO_inv(1+iOffCInv(jSym)),nOrb(jSym),
     &                       Cpq(1+iOff1), nOrb(jSym),
     &                 0.0d0,Cpn, nBas(jSym))
*
*           C_nq^J x  (C^-1) = C_nm^J
*
            iOff2 = nLRb(iSym)*(iJ-1) + iOffLRb(iSym,jSym)
            Call dGemm_('N','N',nBas(jSym),nBas(kSym), nOrb(kSym),
     &                 1.0d0,Cpn,nBas(jSym),
     &                       CMO_inv(1+iOffCInv(kSym)),nOrb(kSym),
     &                 0.0d0,Cmn(1+iOff2), nBas(jSym))
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
     &                 1.0d0,Cmn(1+iOff),nBas(jSym),
     &                       MP2Density(1+iOff0), nBas(kSym),
     &                 0.0d0,Ukn(1+iOff), nBas(jSym))
*
*           D(HF)_kn x C_nm^J = V_km^J, Eq. (42)
*
            iOff0= iOffD(jSym)
            Call dGemm_('N','N',nBas(jSym),nBas(kSym), nBas(jSym),
     &                 1.0d0,SCFDensity(1+iOff0), nBas(jSym),
     &                       Cmn(1+iOff), nBas(jSym),
     &                 0.0d0,Vkn(1+iOff), nBas(jSym))
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
         Call dDaFile(LuLVec,iWr,Cmn,lTot,iAdr)
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuXVec,iWr,Ukn,lTot,iAdr)
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuYVec,iWr,Vkn,lTot,iAdr)
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
            Call dDaFile(lUnit_F(lSym,iTypR),iRd,Ria,lTot,iAdr)

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
     &              1.0d0,CMO_v,nBas(jSym),
     &                    Ria(1+iOff1), nVir(jSym),
     &              0.0d0,Rin(1+iOff2), nBas(jSym))
*
*              R_ni^K x C_mi^T = R_nm^K
*
               iOff3 = iOff_Rmn
               iOff_Rmn = iOff_Rmn + nBas(jSym)*nBas(kSym)
               Call dGemm_('N','T',nBas(jSym),nBas(kSym), nOcc(kSym),
     &                 1.0d0,Rin(1+iOff2),nBas(jSym),
     &                       CMO_o, nBas(kSym),
     &                 0.0d0,Rmn(1+iOff3), nBas(jSym))
*
               End Do ! jSym
            End Do    ! iK
*
*           Write  R_mn^K to disk
*
            lTot = nLRb(iSym)*nK
            iAdr = iAdrR(iSym) +  nLRb(iSym)*(iiK-1)
            Call dDaFile(LuRVec,iWr,Rmn,lTot,iAdr)
*                                                                      *
************************************************************************
*                                                                      *
*           Sum_mn (R_mn^K)^T x C_mn^J = W_KJ  (Eq. 38)
*
            If (iSym.eq.lSym) Then
            iOffZ = iOff_WJL(iSym) + nMP2Vec(iSym)*(iiJ-1) + iiK-1
            Call dGemm_('T','N',nK,nJ,nLRb(iSym),
     &                 1.0d0,Rmn,nLRb(iSym),
     &                       Cmn,nLRb(iSym),
     &                 0.0d0,WJL(1+iOffZ), nMP2Vec(iSym))
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
     &                        1.0d0,Cmn(1+iOffL),nBas(jSym),
     &                              Rin(1+iOffR),nBas(kSym),
     &                        0.0d0,WmjKJ,nBas(jSym))
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
     &                       1.0d0,Rmn(1+iOffR),nBas(mSym2),
     &                             WmjKJ, nBas(jSym),
     &                       Fac,  B3jl(1+iOffB),nBas(mSym2))
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
     &                -8.0d0,B3jl(1+iOffB1) , nBas(jSym),
     &                       CMO_o(1+iOff),nBas(kSym),
     &                 0.0d0,B3kl(1+iOffB2),nBas(jSym))
*
            End Do ! jSym
         End Do    ! iJ
*
*        Write matrix to file. Note that it has the same structure as
*        the R and L matrices.
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuBTmp,iWr,B3kl,lTot,iAdr)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      End Do ! iiJ
      End Do ! iSym
*
      Call mma_deallocate(B3jl)
      Call mma_deallocate(Rin)
      Call mma_deallocate(Cpn)
      Call mma_deallocate(Ria)

      Call mma_deallocate(Cpq)
      Call mma_deallocate(Cab)
      Call mma_deallocate(Cai)
      Call mma_deallocate(CaK)
      Call mma_deallocate(Cia)
      Call mma_deallocate(Cij)
      Call mma_deallocate(CiK)
      Call mma_deallocate(CKa)
      Call mma_deallocate(CKi)
      Call mma_deallocate(CJK)

      Call mma_deallocate(SCFDensity)
*
      iSym = 1
*
      lB1kl = nBas(iSym)*nBas(iSym)*nVec
      Call mma_allocate(B1kl,lB1kl,Label='B1kl')
*
      lB3kl_s = nBas(iSym)*nBas(iSym)*nVec
      Call mma_allocate(B3kl_s,lB3kl_s,Label='B3kl_s')
*
      lA1 = NumCho(iSym)*NumCho(iSym)
      Call mma_allocate(A1,lA1,Label='A1')
      A1(:)=Zero
*
      lA2 = NumCho(iSym)*NumCho(iSym)
      Call mma_allocate(A2,lA2,Label='A2')
      A2(:)=Zero
*
      lB2kl = nBas(iSym)*nBas(iSym)*nVec
      Call mma_allocate(B2kl,lB2kl,Label='B2kl')
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
         Call dDaFile(LuYVec,iRd,Vkn,lTot,iAdr)
*
         Do iiK = 1, NumCho(iSym), nVec
            nK = Min(nVec, NumCho(iSym) - (iiK-1))
*
*           Read U_km_K
*
            lTot = nLRb(iSym)*nK
            iAdr = iAdrB(iSym) + nLRb(iSym)*(iiK-1)
            Call dDaFile(LuXVec,iRd,Ukn,lTot,iAdr)
*
*           U_km^K x V_km^J, 1st term RHS eq. 41

            iOffA = iiK-1 + (iiJ-1)*NumCho(iSym) + iOff_A
            Fact = 1.0D0
            Call dGemm_('T','N', nK, nJ,nLRb(iSym),
     &                 Fact,Ukn,nLRb(iSym),
     &                      Vkn,nLRb(iSym),
     &                0.0d0,A1(1+iOffA),NumCho(iSym))
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
            Call dDaFile(LuRVec,iRd,Rmn,lTot,iAdr)
*
*           R_mn_K x W_JK, 2nd RHS term Eq. 35
*
            Fac = 1.0D0
            If (iiK.eq.1) Fac = 0.0D0
            iOffZ = iOff_WJL(iSym) + iiK-1 + (iiJ-1)*nMp2Vec(iSym)
            Call dGemm_('N','N',nLRb(iSym),nJ,nK,
     &                 -8.0d0,Rmn, nLRb(iSym),
     &                        WJL(1+iOffZ),nMP2Vec(iSym),
     &                  Fac,  B2kl,nLRb(iSym))
         End Do
*
*        Write to disk
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuB(2),iWr,B2kl,lTot,iAdr)
*                                                                      *
************************************************************************
*                                                                      *
*        Read 3rd RHS term in Eq. 40
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuBTmp,iRd,B3kl,lTot,iAdr)
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
                  B3kl_s(1+iOffL+kl_s-1) =
     &               ( B3kl(1+iOffL+kl-1) +
     &                 B3kl(1+iOffL+lk-1) )/2
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
         Call dDaFile(LuBTmp,iWr,B3kl_s,lTot,iAdr)
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
     &                 1.0d0,Vkn(1+iOff1), nBas(jSym),
     &                       MP2Density(1+iOff2), nBas(kSym),
     &                 0.0d0,B1kl(1+iOff1), nBas(jSym))
            End Do
         End Do
*
*        Compound 2nd and 3rd RHS term in Eq. 40.

         Call DaXpY_(nLRb(iSym)*nJ,1.0d0,B3kl,1,B1kl,1)
*
*        Write compounded terms to disk
*
         lTot = nLRb(iSym)*nJ
         iAdr = iAdrB(iSym) + nLRb(iSym)*(iiJ-1)
         Call dDaFile(LuB(1),iWr,B1kl,lTot,iAdr)
*                                                                      *
************************************************************************
*                                                                      *
      End Do ! iiJ
         iOff_A = iOff_A + NumCho(iSym)**2
      End Do ! iSym
*
      Call mma_deallocate(MP2Density)
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
     &            Fact ,WJL(1+iOff1), nMP2Vec(iSym),!-2
     &                  WJL(1+iOff1), nMP2Vec(iSym),
     &            0.0d0,A2(1+iOff_A2),NumCho(iSym))
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
         Call dDaFile(LuBTmp,iRd,B3kl_s,lTot,iAdr)
*
         Do iiK = 1, NumCho(iSym), nVec
            nK = Min(nVec,NumCho(iSym)-(iiK-1))
*
            lTot = nLRb(iSym)*nK
            iAdr = iAdrB(iSym) + nLRb(iSym)*(iiK-1)
            Call dDaFile(LuLVec,iRd,Cmn,lTot,iAdr)
*
*           B_nm^L x C_nm^K
*
*           This is the second term of the RHS of Eq. (41)
*
            iOffA = iiJ-1 + (iiK-1)*NumCho(iSym) + iOff_A
            Fact = 1.0D0
            Call dGemm_('T', 'N', nJ, nK, nLRb(iSym),
     &                 Fact ,B3kl_s,nLRb(iSym),
     &                       Cmn, nLRb(iSym),
     &                 1.0d0,A1(1+iOffA),NumCho(iSym))
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
         Call dDaFile(LuA(1),iWr,A1(1+iOff),lTot,iAdr1)
         Call dDaFile(LuA(2),iWr,A2(1+iOff),lTot,iAdr2)
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
      Call mma_deallocate(Cmn)
      Call mma_deallocate(Rmn)
      Call mma_deallocate(Ukn)
      Call mma_deallocate(Vkn)
*
      Call mma_deallocate(B1kl)
      Call mma_deallocate(A1)
      Call mma_deallocate(B2kl)
      Call mma_deallocate(A2)
      Call mma_deallocate(B3kl)
      Call mma_deallocate(B3kl_s)
      Call mma_deallocate(WmjKJ)
      Call mma_deallocate(WJL)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(CMO_v)
      Call mma_deallocate(CMO_o)
      Call mma_deallocate(CMO_Inv)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
