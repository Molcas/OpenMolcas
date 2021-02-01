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
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Gamma_new()
*                                                                      *
************************************************************************
*                                                                      *
      Implicit Real*8 (a-h,o-z)
*
#include "orbinf2.fh"
#include "corbinf.fh"
#include "WrkSpc.fh"
#include "mp2grad.fh"
#include "files_mbpt2.fh"
      Integer nTOrb(nSym),iOffCMO(nSym),iOffCMO_o(nSym),iOffCMO_v(nSym)
      Logical Triangular, Done, LoadZeros
      Logical NonZeroSym(4)
*                                                                      *
************************************************************************
*                                                                      *
*     Some statement functions
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
      iTable(i,j)=iWork(ipiTab-1+(j-1)*6+i)
      iInt1(i,j,iSyI,iSyJ) =
     &             (nFro(iSyJ) + j-1)*nTOrb(iSyI) +
     &              nFro(iSyI) + nOcc(iSyI) + i-1
      iBinOff(iX,iY,iBin) =  iY-1 + (iX)*2 + (iBin-1)*iBinLength
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
      iAdrGam = 1
      LoadZeros = .false.
*
      iMaxBas = 0
      Do i = 1, nSym
         nTOrb(i) = nOrb(i) + nDel(i)
         iMaxBas = Max(iMaxBas,nBas(i))
      End Do
*
*     Setup one CMO for occupied and one for virtual orbitals together
*     with a corresponding symmetryoffset.
*
      iOffCMO(1) = 0
      iOffCMO_o(1) = 0
      iOffCMO_v(1) = 0
      Do iSym = 2 , nSym
         iOffCMO(iSym)   = iOffCMO(iSym-1)+nTOrb(iSym-1)*nBas(iSym-1)
         iOffCMO_o(iSym) = iOffCMO_o(iSym-1)+nOcc(iSym-1)*nBas(iSym-1)
         iOffCMO_v(iSym) = iOffCMO_v(iSym-1)+nExt(iSym-1)*nBas(iSym-1)
      End Do
*
*     Calculate total length of CMO-matrices in all irreps.
*
      lCMO_o = 0
      lCMO_v = 0
      Do i = 1, nSym
         lCMO_o = lCMO_o + nOcc(i)*nBas(i)
         lCMO_v = lCMO_v + nExt(i)*nBas(i)
      End Do
*
*     Allocate memory for a virtual and an occupied CMO-matrix
*
      Call GetMem('CMO_o','Allo','Real',ipCMO_o,lCMO_o)
      Call GetMem('CMO_v','Allo','Real',ipCMO_v,lCMO_v)
*
*     Copy CMO to CMO_o and CMO_v
*
      iOff = 0
      Do iSym = 1, nSym
         iOff = nBas(iSym)*nFro(iSym)
         nNO=nBas(iSym)*nOcc(iSym)
         nNV=nBas(iSym)*nExt(iSym)
*
         Call dCopy_(nNO,
     &        Work(ipCMO + iOffCMO(iSym)+iOff),1,
     &        Work(ipCMO_o + iOffCMO_o(iSym)),1)
*
         iOff = iOff + nNO
*
         Call dCopy_(nNV,
     &        Work(ipCMO + iOffCMO(iSym)+iOff),1,
     &        Work(ipCMO_v + iOffCMO_v(iSym)),1)
      End Do
*
#ifdef _DEBUGPRINT_
*     Print the elements of the Full CMO-matrices as well as CMO_o and
*     CMO_v.
      Do iSym = 1, nSym
         Call RecPrt('Full CMO',' ',Work(ipCMO + iOffCMO(iSym)),
     &               nBas(iSym),nTOrb(iSym))
      End Do
      Do iSym = 1, nSym
         Call RecPrt('Occupied CMO',' ',Work(ipCMO_o + iOffCMO_o(iSym)),
     &               nBas(iSym),nOcc(iSym))
      End Do
      Do iSym = 1, nSym
         Call RecPrt('Virtual CMO',' ',Work(ipCMO_v + iOffCMO_v(iSym)),
     &               nBas(iSym),nExt(iSym))
      End Do
#endif
*
*     Setup indeces for symmetryblocks to make it compatible with
*     the aces-routine used for reading the gammas later.
*
      If (nSym.eq.8) nBlocks=106
      If (nSym.eq.4) nBlocks= 19
      If (nSym.eq.2) nBlocks=  4
      If (nSym.eq.1) nBlocks=  1
*
*     Construct a table for block information.
*
      Call GetMem('iTable','Allo','Inte',ipiTab,6*nBlocks)
      Call Gamma_Blocks(iWork(ipiTab),nBlocks,nSym)
*
*     Open a file for writing gammas.
*
      iSeed = 10
      LuGam = IsFreeUnit(iSeed)
      FnGam='LuGam   '
      Call DaName_MF_WA(LuGam,FnGam)
*
*     Open a file to store bins of half-transformed gammas.
*
      iSeed = 11
      LuBin = IsFreeUnit(iSeed)
      FnBin='TmpBin  '
      Call DaName_MF_WA(LuBin,FnBin)
*
*     SETUP OF AMPLITUDE BATCHES
*     (The handling is somewhat simple when nBlocks=1 since we only have
*      blocks with one symmetry)
*
*     Check how large the occ*vir-products and bas**2-products may be.
*
      iMaxBasProd = 0
      iMaxOccVir = 0
      Do iSym1 = 1, nSym
         Do iSym2 = 1, iSym1
            iMaxBasProd = Max(iMaxBasProd,nBas(iSym1)*nBas(iSym2))
            iMaxOccVir  = Max(iMaxOccVir,nOcc(iSym1)*nExt(iSym2))
            iMaxOccVir  = Max(iMaxOccVir,nExt(iSym1)*nOcc(iSym2))
         End Do
      End Do
*
*     Setup two pairs of temporary vectors to use for transformation
*     and symmetrization
*
      lTemp = iMaxBasProd
      Call GetMem('Temp1','Allo','Real',ipTemp1,lTemp)
      Call GetMem('Temp2','Allo','Real',ipTemp2,lTemp)
      If(nBlocks .ne. 1) Then
         Call GetMem('Temp1_1','Allo','Real',ipTemp1_2,lTemp)
         Call GetMem('Temp2_2','Allo','Real',ipTemp2_2,lTemp)
      End If
*
*     Check max available memory
*
      Call GetMem('Max','Max','Real',ipMax,lMax)
*
*     Calculate how much space is needed to store all full bins
*     on disk.
*
      If(nBlocks .eq. 1) Then
         iMemNeeded = 2*iMaxBasProd*(iMaxOccVir+1)
      Else
         iMemNeeded = 4*iMaxBasProd*(iMaxOccVir+1)
      End If
*
*     Take maximum one third of the memory for each bin-complex
*     and leave one third for other tasks.
*
      If(nBlocks .eq. 1) Then
         iMemAvail = 2*lMax/3
      Else
         iMemAvail = lMax/3
      End If
*
*     Make the bins as large as possible but not larger than is
*     needed to keep all bins in memory.
*
      lAllBins = min(iMemAvail,iMemNeeded)
*
*     Allocate two bins
*
      Call GetMem('Bins1','Allo','Real',ipBin,lAllBins)
      If(nBlocks .ne. 1) Then
         Call GetMem('Bins2','Allo','Real',ipBin2,lAllBins)
      End If
*
      iBinLength = lAllBins/iMaxBasProd
*
*     The construction and transformation of Tiajb are now made
*     looping over symmetryblocks
*
      Do iBlock = 1, nBlocks
         iType = iTable(1,iBlock)
*
         If((iType.eq.1) .or. (iType.eq.2)) Then
            Triangular = .True.
         Else
            Triangular = .False.
         End If
*
*        nBas, nOcc etc. are defined as nBas(1:8) here in MBPT2
*        and as nBas(0:7) in integral_util and alaska. (which might not
*        be optimal...)
*        The weird order here is because I accidently constructed
*        (kap lam|mu nu) when I needed (mu nu|kap lam), with this order
*        it will be compatible with src/integral_util/read_blocks.f
         iSym_A=iTable(4,iBlock)+1
         iSym_B=iTable(5,iBlock)+1
         iSym_C=iTable(2,iBlock)+1
         iSym_D=iTable(3,iBlock)+1
*
*        Number of orbitals in (ia|jb)
         nI = nOcc(iSym_A)
         nA = nExt(iSym_B)
         nJ = nOcc(iSym_C)
         nB = nExt(iSym_D)
*
*        Number of orbitals in (ai|bj) (nX and nX1 can be
*        combined to get (ia|bj) and (ai|jb) which are needed for
*        full symmetrization)
         nI2 = nOcc(iSym_B)
         nA2 = nExt(iSym_A)
         nJ2 = nOcc(iSym_D)
         nB2 = nExt(iSym_C)
*
*        Initialize adress for storing Bins on disk.
         iAdrBin = 1
*        Setup the number of bins needed.
         If(Triangular) Then
            nBins = nBas(iSym_C)*(nBas(iSym_C)+1)/2
         Else
            nBins = nBas(iSym_C)*nBas(iSym_D)
         End If
*
*        Initialize the bins to have length 0 and adress -1 to
*        next element
         Do i = 1, nBins
            Work(ipBin+iBinOff(0,1,i)) = 0.0d0
            Work(ipBin+iBinOff(0,2,i)) = -1.0d0
            If(nBlocks .ne. 1) Then
               Work(ipBin2+iBinOff(0,1,i)) = 0.0d0
               Work(ipBin2+iBinOff(0,2,i)) = -1.0d0
            End If
         End Do
*
*        Check which symmetry combinations that are present.
         Do i = 1,4
            NonZeroSym(i) = .True.
         End Do
         If(nI*nA*nJ*nB .eq. 0) NonZeroSym(1) = .False.
         If(nI*nA*nJ2*nB2 .eq. 0) NonZeroSym(2) = .False.
         If(nI2*nA2*nJ*nB .eq. 0) NonZeroSym(3) = .False.
         If(nI2*nA2*nJ2*nB2 .eq. 0) NonZeroSym(4) = .False.
*
*        Skip cases where nABCD is zero. If nABCD is nonzero but
*        nIJAB (different A and B) is zero we jump to the end and load
*        some zeros.
*
         If(Triangular) Then
            nABCD = (nBas(iSym_A)*(nBas(iSym_A)+1)/2)*
     &              (nBas(iSym_C)*(nBas(iSym_C)+1)/2)
         Else
            nABCD = nBas(iSym_A)*nBas(iSym_B)*nBas(iSym_C)*nBas(iSym_D)
         End If
*
         If(nABCD .eq. 0) Go To 100
         If(((.not. NonZeroSym(1)) .and. Triangular) .or.
     &       (.not. (NonZeroSym(1) .or. NonZeroSym(2) .or.
     &               NonZeroSym(3) .or. NonZeroSym(4)))) Then
            LoadZeros = .True.
            Go To 500
         End If
*
*        Start the loop to construct Tiajb.
         Do iI = 1, nI
            Do iA = 1, nA
               If(NonZeroSym(1)) Then
                  Call Exch(iSym_D,iSym_A,iSym_C,iSym_B,
     &                 iI+nFro(iSym_A), iA+nFro(iSym_B)+nOcc(iSym_B),
     &                 Work(ipInt1),Work(ipScr1))
                  Call Coul(iSym_D,iSym_C,iSym_A,iSym_B,
     &                 iI+nFro(iSym_A),iA+nFro(iSym_B)+nOcc(iSym_B),
     &                 Work(ipInt2),Work(ipScr1))
#ifdef _DEBUGPRINT_
                  Write(6,*) ' *  I,A = ',iI, iA
                  Call RecPrt('Int1:','(8F10.6)',Work(ipInt1),
     &                 nOrb(iSym_D)+nDel(iSym_D),
     &                 nOrb(iSym_C)+nDel(iSym_C))
                  Write(6,*) ' *  I,A = ',iI, iA
                  Call RecPrt('Int2:','(8F10.6)',Work(ipInt2),
     &                 nOrb(iSym_D)+nDel(iSym_D),
     &                 nOrb(iSym_C)+nDel(iSym_C))
#endif
               End If
*
*              When Symmetry of I and A differs we need to load both
*              possible (IA|xx) to have (mu nu| xx) and (nu mu|xx)
*              available at the same time for symmetrization.
               If(NonZeroSym(2) .and. (.not. Triangular)) Then
                  Call Exch(iSym_C,iSym_A,iSym_D,iSym_B,
     &                 iI+nFro(iSym_A), iA+nFro(iSym_B)+nOcc(iSym_B),
     &                 Work(ipInt1_2),Work(ipScr1))
                  Call Coul(iSym_C,iSym_D,iSym_A,iSym_B,
     &                 iI+nFro(iSym_A),iA+nFro(iSym_B)+nOcc(iSym_B),
     &                 Work(ipInt2_2),Work(ipScr1))
#ifdef _DEBUGPRINT_
                  Write(6,*) ' *  I,A = ',iI, iA
                  Call RecPrt('Int1:','(8F10.6)',Work(ipInt1_2),
     &                 nOrb(iSym_C)+nDel(iSym_C),
     &                 nOrb(iSym_D)+nDel(iSym_D))
                  Write(6,*) ' *  I,A = ',iI, iA
                  Call RecPrt('Int2:','(8F10.6)',Work(ipInt2_2),
     &                 nOrb(iSym_C)+nDel(iSym_C),
     &                 nOrb(iSym_D)+nDel(iSym_D))
#endif
               End If
*
*     Construct Tiajb for a specific (i,a)-pair
               If(NonZeroSym(1)) Then
                  Do iJ = 1, nJ
                     Do iB = 1, nB
                        xiajb = Work(ipInt2 +iInt1(iB,iJ,iSym_D,iSym_C))
                        xibja = Work(ipInt1 +iInt1(iB,iJ,iSym_D,iSym_C))
                        EDenom =(work(mAdOcc(iSym_A)+iI-1)
     &                       + work(mAdOcc(iSym_C)+iJ-1)
     &                       - work(mAdVir(iSym_B)+iA-1)
     &                       - work(mAdVir(iSym_D)+iB-1)  )
                        Tiajb = (2.0d0*xiajb - xibja)/EDenom
                        Work(ipTemp1 + (iJ-1)*nB + iB-1) = Tiajb
                     End Do     !iSymB
                  End Do        !iSymJ
#ifdef _DEBUGPRINT_
               Call RecPrt('Tiajb',' ',Work(ipTemp1),nExt(iSym_D),
     &                     nOcc(iSym_C))
#endif
               End If
*
               If(NonZeroSym(2) .and. (.not. Triangular)) Then
                  Do iJ = 1, nJ2
                     Do iB = 1, nB2
                        xiajb =Work(ipInt2_2+iInt1(iB,iJ,iSym_C,iSym_D))
                        xibja =Work(ipInt1_2+iInt1(iB,iJ,iSym_C,iSym_D))
                        EDenom =(work(mAdOcc(iSym_A)+iI-1)
     &                         + work(mAdOcc(iSym_D)+iJ-1)
     &                         - work(mAdVir(iSym_B)+iA-1)
     &                         - work(mAdVir(iSym_C)+iB-1)  )
                        Tiajb = (2.0d0*xiajb - xibja)/EDenom
                        Work(ipTemp1_2 + (iJ-1)*nB2 + iB-1) = Tiajb
                     End Do     !iSymB
                  End Do        !iSymJ
#ifdef _DEBUGPRINT_
               Call RecPrt('Tiajb',' ',Work(ipTemp1),nExt(iSym_D),
     &                     nOcc(iSym_C))
#endif
               End If
*
*              Do first halftransformation.
               If(NonZeroSym(1)) Then
                  Call dGemm_('N','N',nBas(iSym_D),nOcc(iSym_C),
     &                 nExt(iSym_D),1.0d0,
     &                 Work(ipCMO_v+iOffCMO_v(iSym_D)),nBas(iSym_D),
     &                 Work(ipTemp1),nExt(iSym_D),0.0d0,
     &                 Work(ipTemp2),nBas(iSym_D))
                  Call dGemm_('N','T',nBas(iSym_D),nBas(iSym_C),
     &                 nOcc(iSym_C),1.0d0,
     &                 Work(ipTemp2),nBas(iSym_D),
     &                 Work(ipCMO_o+iOffCMO_o(iSym_C)),nBas(iSym_C),
     &                 0.0d0,Work(ipTemp1),nBas(iSym_D))
               End If
*
               If(NonZeroSym(2) .and. (.not. Triangular)) Then
                  Call dGemm_('N','N',nBas(iSym_C),nOcc(iSym_D),
     &                 nExt(iSym_C),1.0d0,
     &                 Work(ipCMO_v+iOffCMO_v(iSym_C)),nBas(iSym_C),
     &                 Work(ipTemp1_2),nExt(iSym_C),0.0d0,
     &                 Work(ipTemp2_2),nBas(iSym_C))
                  Call dGemm_('N','T',nBas(iSym_C),nBas(iSym_D),
     &                 nOcc(iSym_D),1.0d0,
     &                 Work(ipTemp2_2),nBas(iSym_C),
     &                 Work(ipCMO_o+iOffCMO_o(iSym_D)),nBas(iSym_D),
     &                 0.0d0,Work(ipTemp1_2),nBas(iSym_C))
               End If
*     Do a halfsymmetrization with respect to the two AO-indeces if they are
*     the same symmetry.
               If(Triangular) Then
                  Do iKap = 1,nBas(iSym_D)
                     Do iLam = 1,nBas(iSym_D)
                        iLamKap1 = iLam-1 + (iKap-1)*nBas(iSym_D)
                        iLamKap2 = iKap-1 + (iLam-1)*nBas(iSym_C)
                        Work(ipTemp2+iLamKap1) =
     &                       (Work(ipTemp1+iLamKap1)
     &                       + Work(ipTemp1+iLamKap2))/2
                     End Do
                  End Do
*     Do a halfsymmetrization with respect to the two AO-indeces if they are
*     different symmetries.
               Else
                  Do iKap = 1,nBas(iSym_C)
                     Do iLam = 1,nBas(iSym_D)
                        iLamKap1 = iLam-1 + (iKap-1)*nBas(iSym_D)
                        iLamKap2 = iKap-1 + (iLam-1)*nBas(iSym_C)
                        If(NonZeroSym(1) .and. NonZeroSym(2)) Then
                           Work(ipTemp2 + iLamKap1) =
     &                          (Work(ipTemp1+iLamKap1)
     &                          +  Work(ipTemp1_2 + iLamKap2))/2
                        Else If(NonZeroSym(1)) Then
                           Work(ipTemp2 + iLamKap1) =
     &                          Work(ipTemp1 + iLamKap1)/2
                        Else If(NonZeroSym(2)) Then
                           Work(ipTemp2 + iLamKap1) =
     &                          Work(ipTemp1_2 + iLamKap2)/2
                        End If
                     End Do
                  End Do
               End If
*
*     Place the result in a bin. A bin has fixed Lambda and Kappa and
*     includes Ti,a,lam,kap as well as i,a-adress
               iIA = iA-1 + (iI-1)*nA
               Do iKap = 1, nBas(iSym_C)
                  nLam = nBas(iSym_D)
                  If(Triangular) nLam = iKap
                  Do iLam = 1, nLam
                     iRec = iLam-1 + (iKap-1)*nBas(iSym_D)
                     If(Triangular) Then
                        iBin = iTri(iKap,iLam)
                     Else
                        iBin = iRec+1
                     End If
*                    Write (*,*) 'iKap,iLam,iBin=',iKap,iLam,iBin
                     iNextX = nInt(Work(ipBin+iBinOff(0,1,iBin)))
                     iNextX = iNextX + 1
                     iBinSize = iNextX*2 + 2
                     If(iBinSize .gt. iBinLength) Then
                        iLastAdr = iAdrBin
                        Call dDaFile(LuBin,1,
     &                               Work(ipBin+iBinOff(0,1,iBin)),
     &                               iBinLength,iAdrBin)
                        iNextX = 1
                        Work(ipBin+iBinOff(0,1,iBin)) = 0.0d0
                        Work(ipBin+iBinOff(0,2,iBin)) = dble(iLastAdr)
                     End If
*                    Write (*,*) Work(ipTemp2 + iRec), Dble(iIA), iBin
                     Work(ipBin + iBinOff(iNextX,2,iBin)) =
     &                    Work(ipTemp2 + iRec)
                     Work(ipBin + iBinOff(iNextX,1,iBin)) =
     &                    Dble(iIA)
                     Work(ipBin+iBinOff(0,1,iBin)) = dble(iNextX)
                  End Do
               End Do
            End Do              !iA
         End Do                 !iI
*
*     Now there is a second loop over everything to fix the mu-nu-symmetrization
*
         If(Triangular) Go To 300
         Do iI = 1, nI2
            Do iA =  1, nA2
               If(NonZeroSym(3)) Then
                  Call Exch(iSym_D,iSym_B,iSym_C,iSym_A,
     &                 iI+nFro(iSym_B), iA+nFro(iSym_A)+nOcc(iSym_A),
     &                 Work(ipInt1),Work(ipScr1))
                  Call Coul(iSym_D,iSym_C,iSym_B,iSym_A,
     &                 iI+nFro(iSym_B),iA+nFro(iSym_A)+nOcc(iSym_A),
     &                 Work(ipInt2),Work(ipScr1))
#ifdef _DEBUGPRINT_
                  Write(6,*) ' *  I,A = ',iI, iA
                  Call RecPrt('Int1_lap2:','(8F10.6)',Work(ipInt1),
     &                 nOrb(iSym_D)+nDel(iSym_D),
     &                 nOrb(iSym_C)+nDel(iSym_C))
                  Write(6,*) ' *  I,A = ',iI, iA
                  Call RecPrt('Int2_lap2:','(8F10.6)',Work(ipInt2),
     &                 nOrb(iSym_D)+nDel(iSym_D),
     &                 nOrb(iSym_C)+nDel(iSym_C))
#endif
               End If
*
               If(NonZeroSym(4)) Then
                  Call Exch(iSym_C,iSym_B,iSym_D,iSym_A,
     &                 iI+nFro(iSym_B), iA+nFro(iSym_A)+nOcc(iSym_A),
     &                 Work(ipInt1_2),Work(ipScr1))
                  Call Coul(iSym_C,iSym_D,iSym_B,iSym_A,
     &                 iI+nFro(iSym_B),iA+nFro(iSym_A)+nOcc(iSym_A),
     &                 Work(ipInt2_2),Work(ipScr1))
#ifdef _DEBUGPRINT_
                  Write(6,*) ' *  I,A = ',iI, iA
                  Call RecPrt('Int1_lap2:','(8F10.6)',Work(ipInt1_2),
     &                 nOrb(iSym_C)+nDel(iSym_C),
     &                 nOrb(iSym_D)+nDel(iSym_D))
                  Write(6,*) ' *  I,A = ',iI, iA
                  Call RecPrt('Int2_lap:','(8F10.6)',Work(ipInt2_2),
     &                 nOrb(iSym_C)+nDel(iSym_C),
     &                 nOrb(iSym_D)+nDel(iSym_D))
#endif
               End If
*     Construct Tiajb for a specific (i,a)-pair
               If(NonZeroSym(3)) Then
                  Do iJ = 1, nJ
                     Do iB = 1, nB
                        xiajb = Work(ipInt2 +iInt1(iB,iJ,iSym_D,iSym_C))
                        xibja = Work(ipInt1 +iInt1(iB,iJ,iSym_D,iSym_C))
                        EDenom =(work(mAdOcc(iSym_B)+iI-1)
     &                       +   work(mAdOcc(iSym_C)+iJ-1)
     &                       -   work(mAdVir(iSym_A)+iA-1)
     &                       -   work(mAdVir(iSym_D)+iB-1)  )
                        Tiajb = (2.0d0*xiajb - xibja)/EDenom
                        Work(ipTemp1 + (iJ-1)*nB + iB-1) = Tiajb
                     End Do     !iSymB
                  End Do        !iSymJ
               End If
*
               If(NonZeroSym(4)) Then
                  Do iJ = 1, nJ2
                     Do iB = 1, nB2
                        xiajb =Work(ipInt2_2+iInt1(iB,iJ,iSym_C,iSym_D))
                        xibja =Work(ipInt1_2+iInt1(iB,iJ,iSym_C,iSym_D))
                        EDenom =(work(mAdOcc(iSym_B)+iI-1)
     &                         + work(mAdOcc(iSym_D)+iJ-1)
     &                         - work(mAdVir(iSym_A)+iA-1)
     &                         - work(mAdVir(iSym_C)+iB-1)  )
                        Tiajb = (2.0d0*xiajb - xibja)/EDenom

                        Work(ipTemp1_2 + (iJ-1)*nB2 + iB-1) = Tiajb
                     End Do     !iSymB
                  End Do        !iSymJ
               End If
*     Do first halftransformation.
               If(NonZeroSym(3)) Then
                  Call dGemm_('N','N',nBas(iSym_D),nOcc(iSym_C),
     &                 nExt(iSym_D),1.0d0,
     &                 Work(ipCMO_v+iOffCMO_v(iSym_D)),nBas(iSym_D),
     &                 Work(ipTemp1),nExt(iSym_D),0.0d0,
     &                 Work(ipTemp2),nBas(iSym_D))
                  Call dGemm_('N','T',nBas(iSym_D),nBas(iSym_C),
     &                 nOcc(iSym_C),1.0d0,
     &                 Work(ipTemp2),nBas(iSym_D),
     &                 Work(ipCMO_o+iOffCMO_o(iSym_C)),nBas(iSym_C),
     &                 0.0d0,Work(ipTemp1),nBas(iSym_D))
               End If
*
               If(NonZeroSym(4)) Then
                  Call dGemm_('N','N',nBas(iSym_C),nOcc(iSym_D),
     &                 nExt(iSym_C),1.0d0,
     &                 Work(ipCMO_v+iOffCMO_v(iSym_C)),nBas(iSym_C),
     &                 Work(ipTemp1_2),nExt(iSym_C),0.0d0,
     &                 Work(ipTemp2_2),nBas(iSym_C))
                  Call dGemm_('N','T',nBas(iSym_C),nBas(iSym_D),
     &                 nOcc(iSym_D),1.0d0,
     &                 Work(ipTemp2_2),nBas(iSym_C),
     &                 Work(ipCMO_o+iOffCMO_o(iSym_D)),nBas(iSym_D),
     &                 0.0d0,Work(ipTemp1_2),nBas(iSym_C))
               End If
*     Do a halfsymmetrization with respect to the two AO-indeces if they are
*     the same symmetry.
               Do iKap = 1,nBas(iSym_C)
                  Do iLam = 1,nBas(iSym_D)
                     iLamKap1 = iLam-1 + (iKap-1)*nBas(iSym_D)
                     iLamKap2 = iKap-1 + (iLam-1)*nBas(iSym_C)
                     If(NonZeroSym(3) .and. NonZeroSym(4)) Then
                        Work(ipTemp2 + iLamKap1) =
     &                       (Work(ipTemp1+iLamKap1)
     &                       +  Work(ipTemp1_2 + iLamKap2))/2
                     Else If(NonZeroSym(3)) Then
                        Work(ipTemp2 + iLamKap1) =
     &                       Work(ipTemp1 + iLamKap1)/2
                     Else If(NonZeroSym(4)) Then
                        Work(ipTemp2 + iLamKap1) =
     &                       Work(ipTemp1_2 + iLamKap2)/2
                     End If
                  End Do
               End Do
*
               iIA = iA-1 + (iI-1)*nA2
               Do iKap = 1, nBas(iSym_C)
                  Do iLam = 1, nBas(iSym_D)
                     iRec = iLam-1 + (iKap-1)*nBas(iSym_D)
                     iBin = iRec+1
*                    Write (*,*) 'iKap,iLam,iBin=',iKap,iLam,iBin
                     iNextX = nInt(Work(ipBin2+iBinOff(0,1,iBin)))
                     iNextX = iNextX + 1
*
                     iBinSize = iNextX*2 + 2
                     If(iBinSize .gt. iBinLength) Then
                        Call dDaFile(LuBin,1,
     &                               Work(ipBin2+iBinOff(0,1,iBin)),
     &                               iBinLength,iAdrBin)
                        iNextX = 1
                        Work(ipBin2+iBinOff(0,1,i)) = 0.0d0
                        Work(ipBin2+iBinOff(0,2,i)) = dble(iAdrBin)
                     End If
*
*                    Write (*,*) Work(ipTemp2 + iRec), Dble(iIA), iBin
                     Work(ipBin2 + iBinOff(iNextX,2,iBin)) =
     &                    Work(ipTemp2 + iRec)
                     Work(ipBin2 + iBinOff(iNextX,1,iBin)) =
     &                    Dble(iIA)
                     Work(ipBin2+iBinOff(0,1,iBin)) = dble(iNextX)
                  End Do
               End Do
            End Do              !iA
         End Do                 !iI


 300     Continue

*     Read the halftransformed integrals from one bin at the time
*     and do the two remaining transformations.

*        Write (*,*) 'Do the second half transformation'
 500     Do iBin = 1, nBins
            If(NonZeroSym(1) .or. NonZeroSym(2)) Then
               Done = .false.
 200           If(.not.Done) Then
                  iLen = nInt(Work(ipBin+iBinOff(0,1,iBin)))
                  Do i = 1, iLen
                     iIA = nInt(Work(ipBin+iBinOff(i,1,iBin)))
                     Work(ipTemp1+ iIA) = Work(ipBin+iBinOff(i,2,iBin))
                  End Do

                  iNextAdr = nInt(Work(ipBin+iBinOff(0,2,iBin)))
                  If(iNextAdr.eq.-1) Then
                     Done = .True.
                  Else
                     iAdrRdBin = iNextAdr
                     Call dDaFile(LuBin,2,Work(ipBin+iBinOff(0,1,iBin)),
     &                            iBinLength,iAdrRdBin)
                  End If
                  Go To 200
               End If
*     Do second halftransformation
               Call dGemm_('N','N',nBas(iSym_B),nOcc(iSym_A),
     &              nExt(iSym_B),1.0d0,
     &              Work(ipCMO_v+iOffCMO_v(iSym_B)),nBas(iSym_B),
     &              Work(ipTemp1),nExt(iSym_B),0.0d0,
     &              Work(ipTemp2),nBas(iSym_B))
               Call dGemm_('N','T',nBas(iSym_B),nBas(iSym_A),
     &              nOcc(iSym_A),1.0d0,
     &              Work(ipTemp2),nBas(iSym_B),
     &              Work(ipCMO_o+iOffCMO_o(iSym_A)),nBas(iSym_A),
     &              0.0d0,Work(ipTemp1),nBas(iSym_B))
*              Call RecPrt('(m,n,i,k)',' ',Work(ipTemp1),
*    &                     nBas(iSym_B),nBas(iSym_A))
            End If
*
            If((NonZeroSym(3) .or. NonZeroSym(4)) .and.
     &          (.not.Triangular)) Then
               Done = .false.
 400           If(.not.Done) Then
                  iLen = nInt(Work(ipBin2+iBinOff(0,1,iBin)))
                  Do i = 1, iLen
                     iIA = nInt(Work(ipBin2+iBinOff(i,1,iBin)))
                     Work(ipTemp1_2+ iIA) =
     &                               Work(ipBin2+iBinOff(i,2,iBin))
                  End Do
*
                  iNextAdr = nInt(Work(ipBin2+iBinOff(0,2,iBin)))
                  If(iNextAdr.eq.-1) Then
                     Done = .True.
                  Else
                     iAdrRdBin = iNextAdr
                     Call dDaFile(LuBin,2,
     &                            Work(ipBin2+iBinOff(0,1,iBin)),
     &                            iBinLength,iAdrRdBin)
                  End If
                  Go To 400
               End If
               Call dGemm_('N','N',nBas(iSym_A),nOcc(iSym_B),
     &              nExt(iSym_A),1.0d0,
     &              Work(ipCMO_v+iOffCMO_v(iSym_A)),nBas(iSym_A),
     &              Work(ipTemp1_2),nExt(iSym_A),0.0d0,
     &              Work(ipTemp2_2),nBas(iSym_A))
               Call dGemm_('N','T',nBas(iSym_A),nBas(iSym_B),
     &              nOcc(iSym_B),1.0d0,
     &              Work(ipTemp2_2),nBas(iSym_A),
     &              Work(ipCMO_o+iOffCMO_o(iSym_B)),nBas(iSym_B),
     &              0.0d0,Work(ipTemp1_2),nBas(iSym_A))
            End If
*     Do the second symmetrization. (of the mu and nu-index)
            If(Triangular) Then
               Do iMu = 1,nBas(iSym_A)
                  nNu = nBas(iSym_B)
                  If(Triangular) nNu = iMu
                  Do iNu = 1,nNu
                     iTriMuNu = iTri(iNu,iMu) - 1
                     iMuNu1 = iNu-1 + (iMu-1)*nBas(iSym_B)
                     iMuNu2 = iMu-1 + (iNu-1)*nBas(iSym_A)
                     Work(ipTemp2+iTriMuNu) =
     &                    (Work(ipTemp1+iMuNu1)
     &                    + Work(ipTemp1+iMuNu2))/2
                  End Do
               End Do
            Else
               Do iMu = 1, nBas(iSym_A)
                  Do iNu = 1, nBas(iSym_B)
                     iMuNu1 = iNu-1 + (iMu-1)*nBas(iSym_B)
                     iMuNu2 = iMu-1 + (iNu-1)*nBas(iSym_A)
                     If((NonZeroSym(1) .or. NonZeroSym(2)) .and.
     &                  (NonZeroSym(3) .or. NonZeroSym(4))) Then
                        Work(ipTemp2 + iMuNu1) =
     &                       (Work(ipTemp1+iMuNu1)
     &                       + Work(ipTemp1_2+iMuNu2))/2
                     Else If(NonZeroSym(1) .or. NonZeroSym(2)) Then
                        Work(ipTemp2 + iMuNu1) =
     &                       Work(ipTemp1 + iMuNu1)/2
                     Else If(NonZeroSym(3) .or. NonZeroSym(4)) Then
                        Work(ipTemp2 + iMuNu1) =
     &                       Work(ipTemp1_2 + iMuNu2)/2
                     End If
                  End Do
               End Do
            End If
*           Fix the special case with empty iajb:s but "nonempty"
*           (mu nu|kap lam)
            If(LoadZeros) Then
               Call FZero(Work(ipTemp2), ltemp)
            End If
*
*     Store the result on disk.
            If(Triangular) Then
               iSize = nBas(iSym_A)*(nBas(iSym_A)+1)/2
#ifdef _DEBUGPRINT_
               Call TriPrt('(Gamma)',' ',Work(ipTemp2),
     &                     nBas(iSym_A))
#endif
            Else
               iSize = nBas(iSym_A)*nBas(iSym_B)
#ifdef _DEBUGPRINT_
               Call RecPrt('(Gamma)',' ',Work(ipTemp2),
     &                     nBas(iSym_A),nBas(iSym_B))
#endif
            End If
            Call dDaFile(LuGam,1,Work(ipTemp2),iSize,
     &                   iAdrGam)
         End Do
         LoadZeros = .false.
 100     Continue
      End Do
*
      Call DaClos(LuGam)
      Call DaClos(LuBin)
*
*     Deallocate memory
*
      Call GetMem('iTable','Free','Inte',ipiTab,6*nBlocks)
      Call GetMem('CMO_occ','Free','Real',ipCMO_o,lCMO_o)
      Call GetMem('CMO_vir','Free','Real',ipCMO_v,lCMO_v)
*
      Call GetMem('Temp1','Free','Real',ipTemp1,lTemp)
      Call GetMem('Temp2','Free','Real',ipTemp2,lTemp)
      Call GetMem('Bins','Free','Real',ipBin,lAllBins)
*
      If(nBlocks .ne. 1) Then
         Call GetMem('Temp1_1','Free','Real',ipTemp1_2,lTemp)
         Call GetMem('Temp2_2','Free','Real',ipTemp2_2,lTemp)
         Call GetMem('Bins2','Free','Real',ipBin2,lAllBins)
      End If
*
      Return
      End
