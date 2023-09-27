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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
      SUBROUTINE OLagNS2(iSym,DPT2C,T2AO)
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "caspt2_grad.fh"
C
      DIMENSION DPT2C(*),T2AO(*)
C
      !! orbital Lagrangian from the T-amplitude
      !! See the loop structure in rhs_mp2.f
      !! and the helper subroutine in rhs_mp2_help1/2
C
      nMaxOrb=0
      Do jSym = 1, nSym
        nMaxOrb = Max(nMaxOrb,nBas(jSym))
      End Do
C     write(6,*) "nmaxorb = ", nmaxorb
C
      lInt = nMaxOrb*nMaxOrb
C
      Call GetMem('Int1','Allo','Real',ipInt1,lInt) !! for (ia|jb)
      Call GetMem('Int2','Allo','Real',ipInt2,lInt) !! for (ib|ja)
      Call GetMem('Scr1','Allo','Real',ipScr1,lInt) !! work space

      Call GetMem('Amp1','Allo','Real',ipAmp1,lInt) !! for amplitude
C

      !! (ia|jb)
      Do iSymI = 1, nSym !! Symmetry of occupied (docc+act) orbitals
        !! Check, in particular nFro
        If (nFro(iSymI)+nIsh(iSymI)+nAsh(iSymI).eq.0) Cycle
        Do iSymJ = 1, iSymI
          If (nFro(iSymJ)+nIsh(iSymJ)+nAsh(iSymJ).eq.0) Cycle
          iSymIJ = 1 + iEor(iSymI-1,iSymJ-1)
          Do iSymA = 1, nSym !! Symmetry of non-filled (act+virt) orbs
            If (nAsh(iSymA)+nSsh(iSymA)+nDel(iSymA).eq.0) Cycle
            Do iSymB = 1, iSymA
              If (nAsh(iSymB)+nSsh(iSymB)+nDel(iSymB).eq.0) Cycle
              iSymAB = 1 + iEor(iSymA-1,iSymB-1)
              iSymIJAB = 1 + iEor(iSymIJ-1,iSymAB-1)
              If (iSym.ne.iSymIJAB) Cycle
              Do iCase = 1, 13
C               if (icase.ne.12.and.icase.ne.13) cycle
                Call OLagNs_Hel2(iCase,iSym,iSymA,iSymB,iSymI,iSymJ,
     *                           nMaxOrb,Work(ipInt1),Work(ipInt2),
     *                           Work(ipAmp1),Work(ipScr1),
     *                           DPT2C,T2AO)
              End Do
            End Do
          End Do
        End Do
      End Do
C     write(6,*) "olag before VVVO"
C     call sqprt(work(ipolag),12)
C
      Call DScal_(nBasSq,1.0D+00/DBLE(MAX(1,NACTEL)),DPT2C,1)
C
C
      Call GetMem('Int1','Free','Real',ipInt1,lInt)
      Call GetMem('Int2','Free','Real',ipInt2,lInt)
      Call GetMem('Scr1','Free','Real',ipScr1,lInt)

      Call GetMem('Amp1','Free','Real',ipAmp1,lInt)
C
      END SUBROUTINE OLagNS2
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OLagNS_Hel2(iCase,iSym,iSymA,iSymB,iSymI,iSymJ,nMaxOrb,
     *                       ERI1,ERI2,Amp1,Scr,DPT2C,T2AO)
      USE SUPERINDEX
      USE iSD_data
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "caspt2_grad.fh"
      DIMENSION ERI1(*),ERI2(*),Amp1(nMaxOrb,nMaxOrb),
     *          Scr(nMaxOrb,nMaxOrb)
      DIMENSION DPT2C(*),T2AO(*)
C
      LOGICAL   PM
      DIMENSION IOFF1(8),IOFF2(8)
      !! just to avoid the unused ... of ERI2
      if (.false.) write (6,*) eri2(1)
C
C     DMNS_{ijkl}*d(ij|kl)/dx -> (pj|kl)*D_{qjkl} + (ip|kl)*D_{iqkl}
C                              + (ij|pl)*D_{ijql} + (ij|kp)*D_{ijkq}

C     Integrals needed:
C     (OO|OO), (VO|OO), (VV|OO), (VO|VO), (VV|VO)
C     -> <OO|OO>, <VO|OO>, <VO|VO>, <VV|OO>, <VV|VO>
C     Here, O is occupied (doubly and partially) orbitals
C           V is not filled (partially and virtual) orbitals
C
C     <**|VV> and <**|OO> are fetched by Exch
C     <V*|V*> and <O*|O*> are fetched by Exch
C
C     However, (*O|*O) is split into (*C|*C), (*C|*A), (*A|*A)
C              (V*|V*) is split into (A*|A*), (A*|V*), (V*|V*)
C              (*O|V*) is split into (*A|A*), (*A|V*)
C
C     EXCH(ISYP,ISYI,ISYQ,ISYJ,II,IJ,ERI,SCR)
C
C     rhs_mp2_help1.f
C
      !! The amplitude is in the IC (internally contracted) basis, so
      !! the active orbital index (indices) must be transformed to the
      !! (quasi-)canonical MO (or contravatiant) basis.
      !! IC = SR (why?), contravariant = C
C     write(6,*) "icase = ", icase
C
      PM = .false.
      If (iCase.eq.2.or.iCase.eq.6.or.iCase.eq.8.or.
     *    iCase.eq.10.or.iCase.eq.12) PM = .true.
      If (iCase.eq.3.or.iCase.eq.7.or.iCase.eq.9.or.
     *    iCase.eq.11.or.iCase.eq.13) Return
C
      SQ2    = SQRT(2.0D+00)
      SQI2   = 1.0D+00/SQ2
      SQ3    = SQRT(3.0D+00)
      ! iVec   = iVecX
      IO1=0
      IO2=0
      DO iSymK = 1, nSym
        IOFF1(iSymK) = IO1
        IOFF2(iSymK) = IO2
        iSymAB       = MUL(iSymK,iSym)
        IO1          = IO1+nIsh(iSymK)*nAgeB(iSymAB)
        IO2          = IO2+nIsh(iSymK)*nAgtB(iSymAB)
      END DO
C
      !! Some setup
      !! Read T-amplitude, hopefully in contravariant form
      nINP = 0
      nINM = 0
      nIN  = 0
      If (PM) Then
        nINP = nINDEP(iSym,iCase)
        nASP = nASup(iSym,iCase)
        nISP = nISup(iSym,iCase)
        ! If (nINP.ne.0) Then
        !   nVec = nINP*nISP
        ! End If
        nINM = nINDEP(iSym,iCase+1)
        nASM = nASup(iSym,iCase+1)
        nISM = nISup(iSym,iCase+1)
        ! If (nINM.ne.0) Then
        !   nVec = nINM*nISM
        ! End If
        If (nASP*nISP.ne.0) Then
          Call RHS_ALLO(nASP,nISP,ipTCP)
          CALL RHS_READ_C(ipTCP,iCase,iSym,iVecC2)
        End If
        If (nASM*nISM.ne.0) Then
          Call RHS_ALLO(nASM,nISM,ipTCM)
          CALL RHS_READ_C(ipTCM,iCase+1,iSym,iVecC2)
        End If
      Else
        nIN = nINDEP(iSym,iCase)
        nAS = nASup(iSym,iCase)
        nIS = nISup(iSym,iCase)
        ! If (nIN.ne.0) Then
        !   nVec = nIN*nIS
        ! End If
        If (nAS*nIS.ne.0) Then
          Call RHS_ALLO(nAS,nIS,ipTC)
          CALL RHS_READ_C(ipTC,iCase,iSym,iVecC2)
        End If
      End If
C
      If (iCase.ne.12.and.iCase.ne.13) Then
        If (PM) Then
          If (nASP*nISP.ne.0) Then
            CALL GETMEM('LSTP','ALLO','REAL',LSTP,nASP*nINP)
            idST = idSTMAT(iSym,iCase)
            CALL dDaFile(LUSBT,2,Work(LSTP),nASP*nINP,idST)
          End If
          If (nASM*nISM.ne.0) Then
            CALL GETMEM('LSTM','ALLO','REAL',LSTM,nASM*nINM)
            idST = idSTMAT(iSym,iCase+1)
            CALL dDaFile(LUSBT,2,Work(LSTM),nASM*nINM,idST)
          End If
        Else
          If (nAS*nIN.ne.0) Then
            CALL GETMEM('LST','ALLO','REAL',LST,nAS*nIN)
            idST = idSTMAT(iSym,iCase)
            CALL dDaFile(LUSBT,2,Work(LST),nAS*nIN,idST)
          End If
        End If
      End If
C
C     If (PM) Then
C       If (nINP*nISP.eq.0.and.nINM*nISM.eq.0) GO TO 1
C     Else
C       If (nIN*nIS) GO TO 1
C     End If
C
      Call GetMem('WRK1','Allo','Real',ipWRK1,nBasT*nBasT)
      Call GetMem('WRK2','Allo','Real',ipWRK2,nBasT*nBasT)
C
      nFroI = nFro(iSymI)
      nFroJ = nFro(iSymJ)
      nFroA = nFro(iSymA)
      nFroB = nFro(iSymB)
      nIshI = nIsh(iSymI)
      nIshJ = nIsh(iSymJ)
      nIshA = nIsh(iSymA)
      nIshB = nIsh(iSymB)
      nAshI = nAsh(iSymI)
      nAshJ = nAsh(iSymJ)
      nAshA = nAsh(iSymA)
      nAshB = nAsh(iSymB)
      nSshA = nSsh(iSymA)
      nSshB = nSsh(iSymB)
      nBasI = nBas(iSymI)
      nBasJ = nBas(iSymJ)
      nBasA = nBas(iSymA)
      nBasB = nBas(iSymB)
C
      nCorI = nFroI+nIshI
      nCorJ = nFroJ+nIshJ
      nCorA = nFroA+nIshA
      nCorB = nFroB+nIshB
      nOccA = nCorA+nAshA
      nOccB = nCorB+nAshB
      nOccA2= nOccA-nFroA
      nOccB2= nOccB-nFroB
      nOrbA = nOccA+nSshA
C
      ! nOcc  = nFro(iSym)+nIsh(iSym)+nAsh(iSym)
C
      !! active+virtual part for the right index
C     nJ = nFro(iSymJ)+nIsh(iSymJ)+nAsh(iSymJ)
C
      If (iCase.eq.1) Then
        Call OLagNS_A(Amp1)
      Else If (iCase.eq.2.or.iCase.eq.3) Then
        Call OLagNS_B(Amp1)
      Else If (iCase.eq.4) Then
        Call OLagNS_C(Amp1)
      Else If (iCase.eq.5) Then
        Call OLagNS_D(Amp1)
      Else If (iCase.eq.6.or.iCase.eq.7) Then
        Call OLagNS_E(Amp1)
      Else If (iCase.eq.8.or.iCase.eq.9) Then
        Call OLagNS_F(Amp1)
      Else If (iCase.eq.10.or.iCase.eq.11) Then
        Call OLagNS_G(Amp1)
      Else If (iCase.eq.12.or.iCase.eq.13) Then
        Call OLagNS_H(Amp1)
      End If
C
      Call GetMem('WRK1','Free','Real',ipWRK1,nBasT*nBasT)
      Call GetMem('WRK2','Free','Real',ipWRK2,nBasT*nBasT)
C
C   1 CONTINUE
      If (PM) Then
        If (nASP*nISP.ne.0) Call RHS_FREE(nASP,nISP,ipTCP)
        If (nASM*nISM.ne.0) Call RHS_FREE(nASM,nISM,ipTCM)
      Else
        If (nAS*nIS.ne.0) Call RHS_FREE(nAS,nIS,ipTC)
      End If
C
      If (iCase.ne.12.and.iCase.ne.13) Then
        If (PM) Then
          If (nASP*nISP.ne.0)
     *      CALL GETMEM('LSTP','FREE','REAL',LSTP,nASP*nINP)
          If (nASM*nISM.ne.0)
     *      CALL GETMEM('LSTM','FREE','REAL',LSTM,nASM*nINM)
        Else
          If (nAS*nIN.ne.0) CALL GETMEM('LST','FREE','REAL',LST,nAS*nIN)
        End If
      End If
C
      Return
C
      Contains
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_A(AmpL1)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension AmpL1(nAshA,nAshB)
C
      If (nAshI.eq.0.or.nIshJ.eq.0.or.nAshA.eq.0.or.nAshB.eq.0) Return
C
      ! nJ = nIshJ
      Do iI = 1, nAshI
        iIabs = iI + nIshI + nAes(iSymI)
        ! iItot = iI + nCorI
C       If (iSymI.eq.iSymJ) nJ = iI
        Do iJ = 1, nIshJ
          iJabs = iJ + nIes(iSymJ)
          iJtot = iJ + nFroJ
C         If (iIabs.lt.iJabs) Cycle
          Fac = 1.0D+00
C         If ((iI.ne.iJ).and.(iSymI.eq.iSymJ)) Fac = 2.0d+00
          If (iSymI.ne.iSymJ) Fac = 2.0d+00
C
          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     *              iI+nCorI,iJ+nFroJ,
     *              ERI1,Scr)
C         If ((iI.ne.iJ).or.(iSymI.ne.iSymJ)) then
C           Call Exch(iSymA,iSymJ,iSymB,iSymI,
C    *                iJ+nFroJ,iI+nFroI,
C    *                ERI2,Scr)
C         End If
C
          Call DCopy_(nAshA*nAshB,[0.0D+00],0,AmpL1,1)
C
          Do iA = 1, nAshA
            iAabs = iA + nAes(iSymA)
            ! iAtot = iA + nCorA
            Do iB = 1, nAshB
              iBabs = iB + nAes(iSymB)
              iBtot = iB + nCorB
C
              iTabs = iBabs
              iUabs = iI + nAes(iSymI)
              iVabs = iAabs
              IW1 = kTUV(iTabs,iUabs,iVabs) - nTUVes(iSym)
C             ValA = 0.0D+00
C             Do iICB = 1, nIN
C               iVA  = iICB + nIN*(iJabs-1)
C               ValA  = ValA
C    *        + Work(ipT+iVA-1)*Work(LST+IW1-1+nAS*(iICB-1))
C             End Do
C             ValA = ValA*2.0D+00
              iIS = iJabs
              iAS = IW1
              ValA = Work(ipTC+iAS-1+nAS*(iIS-1))*2.0D+00
C
              If (iUabs.eq.iVabs) Then
                !! For FIMO derivative
                DPT2C(iBtot+nOrbA*(iJtot-1))
     *            = DPT2C(iBtot+nOrbA*(iJtot-1)) + ValA
              End If
C
              AmpL1(iA,iB) = AmpL1(iA,iB) + ValA
            End Do
          End Do
C
          Call DScal_(nAshA*nAshB,Fac,AmpL1,1)
C
          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(1,1,ERI1,AmpL1)
C
          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nAshA,nAshB,nCorA,nCorB,
     *                      AmpL1,Work(ipWRK2))
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,Work(ipWRK2))
        End Do
      End Do
C
      End Subroutine OLagNS_A
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_B(AmpL1)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension AmpL1(nAshA,nAshB)
C
      If (nIshI.eq.0.or.nIshJ.eq.0.or.nAshA.eq.0.or.nAshB.eq.0) Return
C
      nJ = nIshJ
      Do iI = 1, nIshI
        iIabs = iI + nIes(iSymI)
        ! iItot = iI + nFroI
        If (iSymI.eq.iSymJ) nJ = iI
        Do iJ = 1, nJ
          iJabs = iJ + nIes(iSymJ)
          ! iJtot = iJ + nFroJ
          If (iIabs.lt.iJabs) Cycle
          Fac = 1.0D+00
C         If ((iI.ne.iJ).and.(iSymI.eq.iSymJ)) Fac = 2.0d+00
          If (iSymI.ne.iSymJ) Fac = 2.0d+00
C
          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     *              iI+nFroI,iJ+nFroJ,
     *              ERI1,Scr)
C         If ((iI.ne.iJ).or.(iSymI.ne.iSymJ)) then
C           Call Exch(iSymA,iSymJ,iSymB,iSymI,
C    *                iJ+nFroJ,iI+nFroI,
C    *                ERI2,Scr)
C         End If
C
          Call DCopy_(nAshA*nAshB,[0.0D+00],0,AmpL1,1)
C
          Do iA = 1, nAshA
            iAabs = iA + nAes(iSymA)
            ! iAtot = iA + nCorA
            Do iB = 1, nAshB
              iBabs = iB + nAes(iSymB)
              ! iBtot = iB + nCorB
C
              if (iaabs.gt.ibabs) then
                iTabs = iAabs
                iUabs = iBabs
              else
                iTabs = iBabs
                iUabs = iAabs
              end if
              iViP  = kIgeJ(iIabs,iJabs)-nIgeJes(iSym) ! inactive
              iVaP  = kTgeU(iTabs,iUabs)-nTgeUes(iSym) !   active
              iViM  = kIgtJ(iIabs,iJabs)-nIgtJes(iSym)
              iVaM  = kTgtU(iTabs,iUabs)-nTgtUes(iSym)
              !! transform internally contracted (SR)
              !!        to contravariant (C)
              ValBP = 0.0D+00
              ValBM = 0.0D+00
C             Do iICB = 1, nINP
C               iVP  = iICB + nINP*(iViP-1)
C               ValBP = ValBP
C    *        + Work(ipTP+iVP-1)*Work(LSTP+iVaP-1+nASP*(iICB-1))
C             End Do
              iIS = iViP
              iAS = iVaP
              ValBP = Work(ipTCP+iAS-1+nASP*(iIS-1))
              If (iAabs.ne.iBabs.and.iIabs.ne.iJabs) Then
                If (iIabs.ne.iJabs) Then
C                 Do iICB = 1, nINM
C                   iVM  = iICB + nINM*(iViM-1)
C                   ValBM = ValBM
C    *          + Work(ipTM+iVM-1)*Work(LSTM+iVaM-1+nASM*(iICB-1))
C                 End Do
                  iIS = iViM
                  iAS = iVaM
                  ValBM = Work(ipTCM+iAS-1+nASM*(iIS-1))
                End If
                !! permutated
                If (iAabs.lt.iBabs) ValBM = -ValBM
              End If
              If (iIabs.eq.iJabs) ValBP = ValBP*SQI2
C
              AmpL1(iA,iB) = AmpL1(iA,iB) + ValBP + ValBM
            End Do
          End Do
C
          Call DScal_(nAshA*nAshB,Fac,AmpL1,1)
C
          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(1,1,ERI1,AmpL1)
C
          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nAshA,nAshB,nCorA,nCorB,
     *                      AmpL1,Work(ipWRK2))
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,Work(ipWRK2))
        End Do
      End Do
C
      End Subroutine OLagNS_B
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_C(AmpL1)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension AmpL1(nAshA+nSshA,nAshB+nSshB)
C
      If (nAshI.eq.0.or.nAshJ.eq.0.or.nSshA.eq.0.or.nAshB.eq.0) Return
C
      ! nJ = nIshJ
      Do iI = 1, nAshI
        iIabs = iI + nIshI + nAes(iSymI)
        ! iItot = iI + nCorI
C       If (iSymI.eq.iSymJ) nJ = iI
        Do iJ = 1, nAshJ
          iJabs = iJ + nIshJ + nAes(iSymJ)
          ! iJtot = iJ + nCorJ
          If (iIabs.lt.iJabs) Cycle
          Fac = 1.0D+00
C         If ((iI.ne.iJ).and.(iSymI.eq.iSymJ)) Fac = 2.0d+00
          If (iSymI.ne.iSymJ) Fac = 2.0d+00
C
          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     *              iI+nCorI,iJ+nCorJ,
     *              ERI1,Scr)
C         If ((iI.ne.iJ).or.(iSymI.ne.iSymJ)) then
C           Call Exch(iSymA,iSymJ,iSymB,iSymI,
C    *                iJ+nFroJ,iI+nFroI,
C    *                ERI2,Scr)
C         End If
C
          Call DCopy_((nAshA+nSshA)*(nAshB+nSshB),[0.0D+00],0,AmpL1,1)
C
          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            iAtot = iA + nFroA + nIshA + nAshA
            Do iB = 1, nAshB
              iBabs = iB + nAes(iSymB)
              iBtot = iB + nCorB
C
              iTabs = iI + nAes(iSymI)
              iUabs = iBabs
              iVabs = iJ + nAes(iSymJ)
C             write(6,*) itabs,iuabs,ivabs
              !! (at|uv) -> (ai|bj) -> (at|uv)
              ! IW1 = kTUV(iTabs,iUabs,iVabs) - nTUVes(iSym)
              ValC1 = 0.0D+00
              ValC2 = 0.0D+00
C             Do iICB = 1, nIN
C               iV  = iICB + nIN*(iAabs-1)
C               ValC1 = ValC1
C    *        + Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
C             End Do
C             If (iIabs.ne.iJabs) Then
C               IW1 = kTUV(iVabs,iUabs,iTabs) - nTUVes(iSym)
C               Do iICB = 1, nIN
C                 iV  = iICB + nIN*(iAabs-1)
C                 ValC2 = ValC2
C    *          + Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
C               End Do
C             End If
C
C             ValC1 = ValC1*2.0D+00
C             ValC2 = ValC2*2.0D+00
C
              iIS = iAabs
              iAS = kTUV(iTabs,iUabs,iVabs) - nTUVes(iSym)
              ValC1 = Work(ipTC+iAS-1+nAS*(iIS-1))*2.0D+00
              If (iIabs.ne.iJabs) Then
                iAS = kTUV(iVabs,iUabs,iTabs) - nTUVes(iSym)
                ValC2 = Work(ipTC+iAS-1+nAS*(iIS-1))*2.0D+00
              End If
C
              iTabs = iBabs
              iUabs = iI + nAes(iSymI)
              iVabs = iJ + nAes(iSymJ)
              If (iUabs.eq.iVabs) Then
                !! For FIMO derivative
                ONEADD = 0.0D+00
C               IW1 = kTUV(iTabs,iUabs,iVabs) - nTUVes(iSym)
C               Do iICB = 1, nIN
C                 iV  = iICB + nIN*(iAabs-1)
C                 ONEADD = ONEADD
C    *          + Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
C               End Do
C               ONEADD = ONEADD*2.0D+00
                iAS = kTUV(iTabs,iUabs,iVabs) - nTUVes(iSym)
                ONEADD = Work(ipTC+iAS-1+nAS*(iIS-1))*2.0D+00
                DPT2C(iAtot+nOrbA*(iBtot-1))
     *            = DPT2C(iAtot+nOrbA*(iBtot-1)) + ONEADD
C
                !! For -sum(y)(ay,yt) -> (ay,ty) derivative
                !! It is correct, but should be rewritten
C               ONEADD = 0.0D+00
C               Do iXabs = 1, nAshI !?
C                 IW1 = kTUV(iTabs,iXabs,iXabs) - nTUVes(iSym)
C                 Do iICB = 1, nIN
C                   iV  = iICB + nIN*(iAabs-1)
C                   ONEADD = ONEADD
C    *            + Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
C                 End Do
C               End Do
C               ONEADD = 2.0D+00*ONEADD/DBLE(MAX(1,NACTEL))
                ONEADD = 0.0D+00
                Do iXabs = 1, nAshI !?
                  iAS = kTUV(iTabs,iXabs,iXabs) - nTUVes(iSym)
                  ONEADD = ONEADD + Work(ipTC+iAS-1+nAS*(iIS-1))
                End Do
                ONEADD = 2.0D+00*ONEADD/DBLE(MAX(1,NACTEL))
                AmpL1(iAtot-nCorA,iBtot-nCorA)
     *            = AmpL1(iAtot-nCorA,iBtot-nCorA) - ONEADD
              End If
              AmpL1(iAtot-nCorA,iBtot-nCorB)
     *          = AmpL1(iAtot-nCorA,iBtot-nCorB) + ValC1
              AmpL1(iBtot-nCorB,iAtot-nCorA)
     *          = AmpL1(iBtot-nCorB,iAtot-nCorA) + ValC2
            End Do
          End Do
C
          Call DScal_((nAshA+nSshA)*(nAshB+nSshB),Fac,AmpL1,1)
C
          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(3,3,ERI1,AmpL1)
C
          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nAshA+nSshA,nAshB+nSshB,nCorA,nCorB,
     *                      AmpL1,Work(ipWRK2))
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,Work(ipWRK2))
        End Do
      End Do
C
      End Subroutine OLagNS_C
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_D(AmpL1)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension AmpL1(nAshA+nSshA,nAshB+nSshB)
C
      If (nAshI.eq.0.or.nIshJ.eq.0.or.nSshA.eq.0.or.nAshB.eq.0) Return
C
      ! nJ = nIshJ
      Do iI = 1, nAshI
        iIabs = iI + nIshI + nAes(iSymI)
        iItot = iI + nCorI
C       If (iSymI.eq.iSymJ) nJ = iI
        Do iJ = 1, nIshJ
          iJabs = iJ + nIes(iSymJ)
          iJtot = iJ + nFroJ
C         If (iIabs.lt.iJabs) Cycle
          Fac = 1.0D+00
C         If ((iI.ne.iJ).and.(iSymI.eq.iSymJ)) Fac = 2.0d+00
          If (iSymI.ne.iSymJ) Fac = 2.0d+00
C
          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     *              iI+nCorI,iJ+nFroJ,
     *              ERI1,Scr)
C         If ((iI.ne.iJ).or.(iSymI.ne.iSymJ)) then
C           Call Exch(iSymA,iSymJ,iSymB,iSymI,
C    *                iJ+nFroJ,iI+nFroI,
C    *                ERI2,Scr)
C         End If
C
          Call DCopy_((nAshA+nSshA)*(nAshB+nSshB),[0.0D+00],0,AmpL1,1)
C
          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            iAtot = iA + nFroA + nIshA + nAshA
            Do iB = 1, nAshB
              ! iBabs = iB + nAes(iSymB)
              iBtot = iB + nCorB
C
C             iVi   = iJabs + nIshA*(iAabs-1)+IOFF1(iSymA)
C             iVa1  = iB + nAshB*(iI+nAes(iSymI)-1)+IOFF1(iSymB)
C             iVa2  = iB + nAshB*(iI+nAes(iSymI)-1)+IOFF1(iSymB)
C    *                   + nAshT*nAshT
C             !! transform internally contracted (SR)
C             !!        to contravariant (C)
C             ValD1 = 0.0D+00
C             ValD2 = 0.0D+00
C             Do iICB = 1, nIN
C               iVD   = iICB + nIN*(iVi-1)
C               ValD1 = ValD1
C    *        + Work(ipT+iVD-1)*Work(LST+iVa1-1+nAS*(iICB-1))
C               ValD2 = ValD2
C    *        + Work(ipT+iVD-1)*Work(LST+iVa2-1+nAS*(iICB-1))
C             End Do
C             ValD1 = ValD1*2.0D+00
C             ValD2 = ValD2*2.0D+00
              iIS = iJabs + nIshA*(iAabs-1)+iOFF1(iSymA)
              iAS = kTU(iB,iI)-nTUes(iSymA)
              ValD1 = Work(ipTC+iAS-1+nAS*(iIS-1))*2.0D+00
              iAS = iAS + nTU(iSymA)
              ValD2 = Work(ipTC+iAS-1+nAS*(iIS-1))*2.0D+00
C
              !! Fock contributions from the inactive density
              If (iItot.eq.iBtot) Then
                DPT2C(iAtot+nOrbA*(iJtot-1))
     *            = DPT2C(iAtot+nOrbA*(iJtot-1)) + ValD1
              End If
C
              AmpL1(iA+nAshA,iB) = AmpL1(iA+nAshA,iB) + ValD2
              AmpL1(iB,iA+nAshB) = AmpL1(iB,iA+nAshB) + ValD1
            End Do
          End Do
C
          Call DScal_((nAshA+nSshA)*(nAshB+nSshB),Fac,AmpL1,1)
C
          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(3,3,ERI1,AmpL1)
C
          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nAshA+nSshA,nAshB+nSshB,nCorA,nCorB,
     *                      AmpL1,Work(ipWRK2))
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,Work(ipWRK2))
        End Do
      End Do
C
      End Subroutine OLagNS_D
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_E(AmpL1)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension AmpL1(nAshA+nSshA,nAshB+nSshB)
C
      If (nIshI.eq.0.or.nIshJ.eq.0.or.nSshA.eq.0.or.nAshB.eq.0) Return
C
      nJ = nIshJ
      Do iI = 1, nIshI
        iIabs = iI + nIes(iSymI)
        ! iItot = iI
        If (iSymI.eq.iSymJ) nJ = iI
        Do iJ = 1, nJ
          iJabs = iJ + nIes(iSymJ)
          ! iJtot = iJ
          If (iIabs.lt.iJabs) Cycle
          Fac = 1.0D+00
C         If ((iI.ne.iJ).and.(iSymI.eq.iSymJ)) Fac = 2.0d+00
          If (iSymI.ne.iSymJ) Fac = 2.0d+00
C
          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     *              iI+nFroI,iJ+nFroJ,
     *              ERI1,Scr)
C         If ((iI.ne.iJ).or.(iSymI.ne.iSymJ)) then
C           Call Exch(iSymA,iSymJ,iSymB,iSymI,
C    *                iJ+nFroJ,iI+nFroI,
C    *                ERI2,Scr)
C         End If
C
          Call DCopy_((nAshA+nSshA)*(nAshB+nSshB),[0.0D+00],0,AmpL1,1)
C
          IgeJ  = kIgeJ(iIabs,iJabs) - nIgeJes(iSym) ! iSymIJ
          IgtJ  = kIgtJ(iIabs,iJabs) - nIgtJes(iSym) ! iSymIJ
          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            ! iAtot = iA + nIshA + nAshA
            Do iB = 1, nAshB
              iBabs = iB + nAes(iSymB)
              ! iBtot = iB + nIshB
C
              iASP  = iBabs
              iISP  = iAabs + nSshA*(IgeJ-1)+iOFF1(iSymA)
              ValEP = Work(ipTCP+iASP-1+nASP*(iISP-1))
              ValEM = 0.0D+00
              If (iIabs.gt.iJabs) Then
                ValEP = ValEP * SQ2
                iASM  = iBabs
                iISM  = iAabs + nSshA*(IgtJ-1)+iOFF1(iSymA)
                ValEM = Work(ipTCM+iASM-1+nASM*(iISM-1))*SQ2*SQ3
              Else
              End If
C
              AmpL1(iA+nAshA,iB) = AmpL1(iA+nAshA,iB) + ValEP + ValEM
              AmpL1(iB,iA+nAshB) = AmpL1(iB,iA+nAshB) + ValEP - ValEM
            End Do
          End Do
C
          Call DScal_((nAshA+nSshA)*(nAshB+nSshB),Fac,AmpL1,1)
C
          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(3,3,ERI1,AmpL1)
C
          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nAshA+nSshA,nAshB+nSshB,nCorA,nCorB,
     *                      AmpL1,Work(ipWRK2))
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,Work(ipWRK2))
        End Do
      End Do
C
      End Subroutine OLagNS_E
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_F(AmpL1)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension AmpL1(nSshA,nSshB)
C
      If (nAshI.eq.0.or.nAshJ.eq.0.or.nSshA.eq.0.or.nSshB.eq.0) Return
C
      ! nJ = nIshJ
      Do iI = 1, nAshI
        iIabs = iI + nIshI + nAes(iSymI)
        ! iItot = iI + nCorI
C       If (iSymI.eq.iSymJ) nJ = iI
        Do iJ = 1, nAshI
          iJabs = iJ + nIshJ + nIes(iSymJ)
          ! iJtot = iJ + nCorJ
          If (iIabs.lt.iJabs) Cycle
          Fac = 1.0D+00
C         If ((iI.ne.iJ).and.(iSymI.eq.iSymJ)) Fac = 2.0d+00
          If (iSymI.ne.iSymJ) Fac = 2.0d+00
C
          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     *              iI+nCorI,iJ+nCorJ,
     *              ERI1,Scr)
C         If ((iI.ne.iJ).or.(iSymI.ne.iSymJ)) then
C           Call Exch(iSymA,iSymJ,iSymB,iSymI,
C    *                iJ+nFroJ,iI+nFroI,
C    *                ERI2,Scr)
C         End If
C
          Call DCopy_(nSshA*nSshB,[0.0D+00],0,AmpL1,1)
C
          iTabs = iI + nAes(iSymI)
          iUabs = iJ + nAes(iSymJ)
          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            ! iAtot = iA + nIsh(iSymA) + nAsh(iSymA)
            Do iB = 1, nSshB
              iBabs = iB + nSes(iSymB)
              If (iAabs.lt.iBabs) Cycle
              ! iBtot = iB + nIsh(iSymB) + nAsh(iSymB)
C
              iASP  = kTgeU(iTabs,iUabs)-nTgeUes(iSym)
              iISP  = kAgeB(iAabs,iBabs)-nAgeBes(iSym)
              ValFP = Work(ipTCP+iASP-1+nASP*(iISP-1))
              If (iIabs.eq.iJabs) ValFP = ValFP*0.5D+00
              ValFM = 0.0D+00
              If (iAabs.ne.iBabs) Then
                If (iTabs.ne.iUabs) Then
                  iASM  = kTgtU(iTabs,iUabs)-nTgtUes(iSym)
                  iISM  = kAgtB(iAabs,iBabs)-nAgtBes(iSym)
                  ValFM = Work(ipTCM+iASM-1+nASM*(iISM-1))
                End If
              Else
                ValFP = ValFP * SQI2
              End If
              VALFM = -VALFM !! why?
C
              AmpL1(iA,iB) = AmpL1(iA,iB) + ValFP + ValFM
              AmpL1(iB,iA) = AmpL1(iB,iA) + ValFP - ValFM
            End Do
          End Do
C
          Call DScal_(nSshA*nSshB,Fac,AmpL1,1)
C
          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(2,2,ERI1,AmpL1)
C
          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nSshA,nSshB,nOccA,nOccB,AmpL1,Work(ipWRK2))
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,Work(ipWRK2))
        End Do
      End Do
C
      End Subroutine OLagNS_F
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_G(AmpL1)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension AmpL1(nSshA,nSshB)
C
      If (nAshI.eq.0.or.nIshJ.eq.0.or.nSshA.eq.0.or.nSshB.eq.0) Return
C
      ! nJ = nIshJ
      Do iI = 1, nAshI
        iIabs = iI + nIshI + nAes(iSymI)
        ! iItot = iI + nCorI
C       If (iSymI.eq.iSymJ) nJ = iI
        Do iJ = 1, nIshI
          iJabs = iJ + nIes(iSymJ)
          ! iJtot = iJ
C         If (iIabs.lt.iJabs) Cycle
          Fac = 1.0D+00
C         If ((iI.ne.iJ).and.(iSymI.eq.iSymJ)) Fac = 2.0d+00
          If (iSymI.ne.iSymJ) Fac = 2.0d+00
C
          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     *              iI+nCorI,iJ+nFroJ,
     *              ERI1,Scr)
C      write(6,*) "integral",ii+nfroi,ij+nfroj
C      call sqprt(eri1,12)
C         If ((iI.ne.iJ).or.(iSymI.ne.iSymJ)) then
C           Call Exch(iSymA,iSymJ,iSymB,iSymI,
C    *                iJ+nFroJ,iI+nFroI,
C    *                ERI2,Scr)
C         End If
C
          Call DCopy_(nSshA*nSshB,[0.0D+00],0,AmpL1,1)
C
          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            ! iAtot = iA + nIsh(iSymA) + nAsh(iSymA)
            Do iB = 1, nSshB
              iBabs = iB + nSes(iSymB)
              If (iAabs.lt.iBabs) Cycle
              ! iBtot = iB + nIsh(iSymB) + nAsh(iSymB)
C
              iAgeB = kAgeB(iAabs,iBabs)-nAgeBes(iSym) !! iSymAB
              iVjP  = iJ + nIsh(iSymJ)*(iAgeB-1)+IOFF1(iSymJ)
              ValGP = Work(ipTCP+iI-1+nASP*(iVjP-1))
              ValGM = 0.0D+00
              If (iAabs.ne.iBabs) Then
                ValGP = ValGP * SQ2
                iAgtB = kAgtB(iAabs,iBabs) - nAgtBes(iSym) !! iSymAB
                iVjM  = iJ + nIsh(iSymJ)*(iAgtB-1)+IOFF2(iSymJ)
                ValGM = Work(ipTCM+iI-1+nASM*(iVjM-1))*SQ2*SQ3
              End If
C
              AmpL1(iA,iB) = AmpL1(iA,iB) + ValGP + ValGM
              AmpL1(iB,iA) = AmpL1(iB,iA) + ValGP - ValGM
            End Do
          End Do
C
          Call DScal_(nSshA*nSshB,Fac,AmpL1,1)
C
          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(2,2,ERI1,AmpL1)
C         Call DGEMM_('N','T',nOrbA,nSshA,nSshB,
C    *                1.0D+00,ERI1(1+nOrbA*nOccB),nOrbA,
C    *                        AmpL1,nSshA,
C    *                1.0D+00,Work(ipOLAG+nOrbA*nOccB),nOrbA)
C         Call DGEMM_('T','N',nOrbA,nSshA,nSshB,
C    *                1.0D+00,ERI1(nOccA+1),nOrbA,
C    *                        AmpL1,nSshA,
C    *                1.0D+00,Work(ipOLAG+nOrbA*nOccB),nOrbA)
C
          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu} for all ij
          Call OLagNS_post2(nSshA,nSshB,nOccA,nOccB,AmpL1,Work(ipWRK2))
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,Work(ipWRK2))
        End Do
      End Do
C
      End Subroutine OLagNS_G
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_H(AmpL1)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension AmpL1(nSshA,nSshB)
C
      If (nIshI.eq.0.or.nIshJ.eq.0.or.nSshA.eq.0.or.nSshB.eq.0) Return
C
      nJ = nIshJ
      Do iI = 1, nIshI
        iIabs = iI + nIes(iSymI)
        ! iItot = iI
        If (iSymI.eq.iSymJ) nJ = iI
        Do iJ = 1, nJ
          iJabs = iJ + nIes(iSymJ)
          ! iJtot = iJ
          If (iIabs.lt.iJabs) Cycle
          Fac = 1.0D+00
C         If ((iI.ne.iJ).and.(iSymI.eq.iSymJ)) Fac = 2.0d+00
          If (iSymI.ne.iSymJ) Fac = 2.0d+00
C
          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     *              iI+nFroI,iJ+nFroJ,
     *              ERI1,Scr)
C         If ((iI.ne.iJ).or.(iSymI.ne.iSymJ)) then
C           Call Exch(iSymA,iSymJ,iSymB,iSymI,
C    *                iJ+nFroJ,iI+nFroI,
C    *                ERI2,Scr)
C         End If
C
          Call DCopy_(nSshA*nSshB,[0.0D+00],0,AmpL1,1)
C
          iViHP0= kIgeJ(iIabs,iJabs) - nIgeJes(iSym)
          iViHP = nAgeB(iSym)*(iViHP0-1)
          iViHM0= kIgtJ(iIabs,iJabs) - nIgtJes(iSym)
          iViHM = nAgtB(iSym)*(iViHM0-1)
          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            ! iAtot = iA + nIsh(iSymA) + nAsh(iSymA)
            Do iB = 1, nSshB
              iBabs = iB + nSes(iSymB)
              If (iAabs.lt.iBabs) Cycle
              ! iBtot = iB + nIsh(iSymB) + nAsh(iSymB)
              iVaHP = kAgeB(iAabs,iBabs) - nAgeBes(iSym)
              iVHP  = iVaHP + iViHP !! nAgeB(iSym)*(iViP-1)
C
              ValHP = Work(ipTCP+iVHP-1)
              ValHM = 0.0D+00
              If (iIabs.ne.iJabs) Then
                If (iAabs.ne.iBabs) Then
                  ValHP = ValHP * 2.0D+00
                  iVaHM = kAgtB(iAabs,iBabs) - nAgtBes(iSym)
                  iVHM  = iVaHM + iViHM !! nAgtB(iSym)*(iViM-1)
                  ValHM = Work(ipTCM+iVHM-1) * 2.0D+00*SQ3
                Else
                  ValHP = ValHP * SQ2
                End If
              Else
                If (iAabs.ne.iBabs) ValHP = ValHP * SQ2
              End If
C
C     write(6,'(2i3,2f20.10)')ia,ib,valhp,valhm
              AmpL1(iA,iB) = AmpL1(iA,iB) + ValHP + ValHM
              AmpL1(iB,iA) = AmpL1(iB,iA) + ValHP - ValHM
            End Do
          End Do
C
          Call DScal_(nSshA*nSshB,Fac,AmpL1,1)
C
          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(2,2,ERI1,AmpL1)
C         Call DGEMM_('N','T',nOrbA,nSshA,nSshB,
C    *                1.0D+00,ERI1(1+nOrbA*nOccB),nOrbA,
C    *                        AmpL1,nSshA,
C    *                1.0D+00,Work(ipOLAG+nOrbA*nOccB),nOrbA)
C         Call DGEMM_('T','N',nOrbA,nSshA,nSshB,
C    *                1.0D+00,ERI1(nOccA+1),nOrbA,
C    *                        AmpL1,nSshA,
C    *                1.0D+00,Work(ipOLAG+nOrbA*nOccB),nOrbA)
C
          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu} for all ij
          Call OLagNS_post2(nSshA,nSshB,nOccA,nOccB,AmpL1,Work(ipWRK2))
C         Call DGEMM_('N','N',nBasT,nSshA,nSshB,
C    *                1.0D+00,Work(LCMOPT2+nBasT*nOccB),nBasT,
C    *                        AmpL1,nSshA,
C    *                0.0D+00,Work(ipWRK1),nBasT)
C         Call DGEMM_('N','T',nBasT,nBasT,nSshB,
C    *                1.0D+00,Work(ipWRK1),nBasT,
C    *                        Work(LCMOPT2+nBasT*nOccA),nBasT,
C    *                0.0D+00,Work(ipWRK2),nBasT)
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,Work(ipWRK2))
C         Do iBas = 1, nBasT
C           Do jBas = 1, nBasT
C             loc1 = iJ-1 + (jBas-1)*nOccA2
C    *             + (iI-1)*nOccA2*nBasT + (iBas-1)*nOccA2*nBasT*nOccA2
C             loc2 = iI-1 + (jBas-1)*nOccA2
C    *             + (iJ-1)*nOccA2*nBasT + (iBas-1)*nOccA2*nBasT*nOccA2
C             loc3 = iBas-1 + (jBas-1)*nBasT
C             loc4 = jBas-1 + (iBas-1)*nBasT
C             T2AO(1+loc1) = T2AO(1+loc1) + Work(ipWRK2+loc3)
C             T2AO(1+loc2) = T2AO(1+loc2) + Work(ipWRK2+loc4)
C           End Do
C         End Do
        End Do
      End Do
C
      End Subroutine OLagNS_H
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_post1(iLeft,iRight,ERI,AmpMO)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension ERI(*),AmpMO(*)
C
      nSkpA = 0
      nSkpB = 0
      nDimA = 0
      nDimB = 0
      If (iLeft .eq.1) Then
        nDimA = nAshA
        nSkpA = nCorA
      Else If (iLeft.eq.2) Then
        nDimA = nSshA
        nSkpA = nOccA
      Else If (iLeft.eq.3) Then
        nDimA = nAshA+nSshA
        nSkpA = nCorA
      End If
      If (iRight.eq.1) Then
        nDimB = nAshB
        nSkpB = nCorB
      Else If (iRight.eq.2) Then
        nDimB = nSshB
        nSkpB = nOccB
      Else If (iRight.eq.3) Then
        nDimB = nAshB+nSshB
        nSkpB = nCorB
      End If
C
      Call DGEMM_('N','T',nOrbA,nDimA,nDimB,
     *            1.0D+00,ERI(1+nOrbA*nSkpB),nOrbA,
     *                    AmpMO,nDimA,
     *            1.0D+00,Work(ipOLAG+nOrbA*nSkpB),nOrbA)
      Call DGEMM_('T','N',nOrbA,nDimA,nDimB,
     *            1.0D+00,ERI(nSkpA+1),nOrbA,
     *                    AmpMO,nDimA,
     *            1.0D+00,Work(ipOLAG+nOrbA*nSkpB),nOrbA)
C
      End Subroutine OLagNS_post1
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_post2(nDimA,nDimB,nSkpA,nSkpB,AmpMO,AmpAO)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension AmpMO(nDimA,nDimB),AmpAO(nBasA,nBasB)
C
       Call DGEMM_('N','N',nBasA,nDimB,nDimA,
     *             1.0D+00,Work(LCMOPT2+nBasA*nSkpA),nBasA,
     *                     AmpMO,nDimA,
     *             0.0D+00,Work(ipWRK1),nBasA)
       Call DGEMM_('N','T',nBasA,nBasB,nDimB,
     *             1.0D+00,Work(ipWRK1),nBasA,
     *                     Work(LCMOPT2+nBasB*nSkpB),nBasA,
     *             0.0D+00,AmpAO,nBasA)
C
      End Subroutine OLagNS_post2
C
C-----------------------------------------------------------------------
C
      Subroutine OLagNS_post3(iIabs,iJabs,TampAO,TampIJ)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension TampAO(nOccA2,nBasI,nOccB2,nBasJ),TampIJ(nBasI,nBasJ)
C
      Do iBas = 1, nBasI
        Do jBas = 1, nBasJ
          TampAO(iJabs,jBas,iIabs,iBas)
     *      = TampAO(iJabs,jBas,iIabs,iBas) + TampIJ(iBas,jBas)
          TampAO(iIabs,jBas,iJabs,iBas)
     *      = TampAO(iIabs,jBas,iJabs,iBas) + TampIJ(jBas,iBas)
        End Do
      End Do
C     Do iBas = 1, nBasT
C       Do jBas = 1, nBasT
C         loc1 = iJ-1 + (jBas-1)*nOccA2
C    *         + (iI-1)*nOccA2*nBasT + (iBas-1)*nOccA2*nBasT*nOccA2
C         loc2 = iI-1 + (jBas-1)*nOccA2
C    *         + (iJ-1)*nOccA2*nBasT + (iBas-1)*nOccA2*nBasT*nOccA2
C         loc3 = iBas-1 + (jBas-1)*nBasT
C         loc4 = jBas-1 + (iBas-1)*nBasT
C         T2AO(1+loc1) = T2AO(1+loc1) + Work(ipWRK2+loc3)
C         T2AO(1+loc2) = T2AO(1+loc2) + Work(ipWRK2+loc4)
C       End Do
C     End Do
C
      End Subroutine OLagNS_post3
C
      END SUBROUTINE OLagNS_Hel2
C
C-----------------------------------------------------------------------
C
! MO->AO or AO->MO transformation of 1-RDM
      Subroutine OLagTrf(mode,iSym,CMO,DPT2,DPT2AO,WRK)
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
C
      Dimension CMO(*),DPT2(*),DPT2AO(*),WRK(*)
C
      !! Mode = 1: MO -> AO transformation
      !! Mode = 2: AO -> MO transformation
      iCMO =1
      iAO = 1
      iMO = 1
      Do jSym = 1, iSym-1
        iCMO = iCMO + nBas(jSym)**2 !! ??
        iAO  = iAO  + nBas(jSym)**2
        iMO  = iMO  + (nOrb(jSym)+nFro(jSym))**2
      End Do
C
      If (nOrb(iSym)+nFro(iSym).GT.0) Then
        nBasI = nBas(iSym)
        nOrbI = nBas(iSym)-nDel(iSym)
        If (Mode.eq.1) Then
          !! MO -> AO
          CALL DGEMM_('N','N',nBasI,nOrbI,nOrbI,
     *                1.0D+00,CMO(iCMO),nBasI,DPT2(iMO),nOrbI,
     *                0.0D+00,WRK,nBasI)
          CALL DGEMM_('N','T',nBasI,nBasI,nOrbI,
     *                1.0D+00,WRK,nBasI,CMO(iCMO),nBasI,
     *                0.0D+00,DPT2AO(iAO),nBasI)
          !! Symmetrize, just in case
          Do iBas = 1, nBasI
            Do jBas = 1, iBas-1
              Val =(DPT2AO(iAO+iBas-1+nBasI*(jBas-1))
     *            + DPT2AO(iAO+jBas-1+nBasI*(iBas-1)))*0.5D+00
              DPT2AO(iAO+iBas-1+nBasI*(jBas-1)) = Val
              DPT2AO(iAO+jBas-1+nBasI*(iBas-1)) = Val
            End Do
          End Do
        Else If (Mode.eq.2) Then
          !! AO -> MO
          CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,
     *                1.0D+00,CMO(iCMO),nBasI,DPT2AO(iAO),nBasI,
     *                0.0D+00,WRK,nOrbI)
          CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,
     *                1.0D+00,WRK,nOrbI,CMO(iCMO),nBasI,
     *                0.0D+00,DPT2(iMO),nOrbI)
        End If
      END IF
C
      Return
C
      End Subroutine OLagTrf
