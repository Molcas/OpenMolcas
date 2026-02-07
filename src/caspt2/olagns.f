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

      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: NSYM, NACTEL, NFRO, NISH, NASH, NSSH,
     &                         NDEL, NBAS, NBSQT
      use Constants, only: One

      implicit none

      integer(kind=iwp), intent(in) :: iSym
      real(kind=wp), intent(inout) :: DPT2C(*), T2AO(*)

      real(kind=wp), allocatable :: Int1(:), Scr1(:), Amp1(:)
      integer(kind=iwp) :: nMaxOrb, jSym, lInt, iSymI, iSymJ, iSymIJ,
     &  iSymA, iSymB, iSymAB, iSymIJAB, iCase

      !! orbital Lagrangian from the T-amplitude
      !! See the loop structure in rhs_mp2.f
      !! and the helper subroutine in rhs_mp2_help1/2

      nMaxOrb=0
      Do jSym = 1, nSym
        nMaxOrb = Max(nMaxOrb,nBas(jSym))
      End Do
      lInt = nMaxOrb*nMaxOrb

      call mma_allocate(Int1,lInt,Label='Int1') !! for (ia|jb)
      call mma_allocate(Scr1,lInt,Label='Scr1') !! work space
      call mma_allocate(Amp1,lInt,Label='Amp1') !! for amplitude

      !! (ia|jb)
      Do iSymI = 1, nSym !! Symmetry of occupied (docc+act) orbitals
        !! Check, in particular nFro
        If (nFro(iSymI)+nIsh(iSymI)+nAsh(iSymI) == 0) Cycle
        Do iSymJ = 1, iSymI
          If (nFro(iSymJ)+nIsh(iSymJ)+nAsh(iSymJ) == 0) Cycle
          iSymIJ = 1 + iEor(iSymI-1,iSymJ-1)
          Do iSymA = 1, nSym !! Symmetry of non-filled (act+virt) orbs
            If (nAsh(iSymA)+nSsh(iSymA)+nDel(iSymA) == 0) Cycle
            Do iSymB = 1, iSymA
              If (nAsh(iSymB)+nSsh(iSymB)+nDel(iSymB) == 0) Cycle
              iSymAB = 1 + iEor(iSymA-1,iSymB-1)
              iSymIJAB = 1 + iEor(iSymIJ-1,iSymAB-1)
              If (iSym /= iSymIJAB) Cycle
              Do iCase = 1, 13
                Call OLagNs_Hel2(iCase,iSym,iSymA,iSymB,iSymI,iSymJ,
     &                           nMaxOrb,Int1,Amp1,Scr1,DPT2C,T2AO)
              End Do
            End Do
          End Do
        End Do
      End Do

      DPT2C(1:NBSQT) = DPT2C(1:NBSQT)/real(max(1,NACTEL),kind=wp)

      call mma_deallocate(Int1)
      call mma_deallocate(Scr1)
      call mma_deallocate(Amp1)

      END SUBROUTINE OLagNS2
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OLagNS_Hel2(iCase,iSym,iSymA,iSymB,iSymI,iSymJ,nMaxOrb,
     &                       ERI1,Amp1,Scr,DPT2C,T2AO)

      USE SUPERINDEX, only: KTU, KTUV, KTGEU, KTGTU, KAGEB, KAGTB,
     &                      KIGEJ, KIGTJ
      use caspt2_global, only: OLag
      use EQSOLV, only: IVECC2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use fake_GA, only: GA_Arrays
      use caspt2_module, only: NACTEL, NSYM, NFRO, NISH, NIES, NASH,
     &                         NAES, NSSH, NSES, NBAS, NBAST, MUL, NTU,
     &                         NAGEB, NAGTB, NTUVES, NTUES, NTGEUES,
     &                         NTGTUES, NIGEJES, NIGTJES, NAGEBES,
     &                         NAGTBES, NASUP, NISUP
      use Constants, only: Zero, One, Half, Two, Three

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: iCase, iSym, iSymA, iSymB, iSymI,
     &  iSymJ, nMaxOrb
      real(kind=wp), intent(_OUT_) :: ERI1(*)
      real(kind=wp), intent(out) :: Amp1(nMaxOrb,nMaxOrb),
     &  Scr(nMaxOrb,nMaxOrb)
      real(kind=wp), intent(inout) :: DPT2C(*), T2AO(*)

      real(kind=wp), allocatable :: WRK1(:), WRK2(:)

      logical(kind=iwp) :: PM
      integer(kind=iwp) :: IOFF1(8), IOFF2(8), IO1, IO2, iSymK, iSymAB,
     &  nASP, nISP, nASM, nISM, ipTCP, ipTCM, nAS, nIS, ipTC, nFroI,
     &  nFroJ, nFroA, nFroB, nIshI, nIshJ, nIshA, nIshB, nAshI, nAshJ,
     &  nAshA, nAshB, nSshA, nSshB, nBasI, nBasJ, nBasA, nBasB,
     &  nCorI, nCorJ, nCorA, nCorB, nOccA, nOccB, nOccA2, nOccB2,
     &  nOrbA
      integer(kind=iwp) :: iI, iIabs, iJ, iJabs, iJtot, iA, iAabs,
     &  iB, iBabs, iBtot, iTabs, iUabs, iVabs, IW1, iIS, iAS, nJ,
     &  iViP, iVaP, iViM, iVaM, iAtot, iXabs, iItot, IgeJ, IgtJ, iASP,
     &  iISP, iASM, iISM, iAgeB, iVjP, iAgtB, iVjM, iViHP0, iViHP,
     &  iViHM0, iViHM, iVaHP, iVHP, iVaHM, iVHM
      real(kind=wp) :: SQ2, SQI2, SQ3, Fac
      real(kind=wp) :: ValA, ValBP, ValBM, ValC1, ValC2, ONEADD, ValD1,
     &  ValD2, ValEP, ValEM, ValFP, ValFM, ValGP, ValGM, ValHP, ValHM
!
!     DMNS_{ijkl}*d(ij|kl)/dx -> (pj|kl)*D_{qjkl} + (ip|kl)*D_{iqkl}
!                              + (ij|pl)*D_{ijql} + (ij|kp)*D_{ijkq}

!     Integrals needed:
!     (OO|OO), (VO|OO), (VV|OO), (VO|VO), (VV|VO)
!     -> <OO|OO>, <VO|OO>, <VO|VO>, <VV|OO>, <VV|VO>
!     Here, O is occupied (doubly and partially) orbitals
!           V is not filled (partially and virtual) orbitals
!
!     <**|VV> and <**|OO> are fetched by Exch
!     <V*|V*> and <O*|O*> are fetched by Exch
!
!     However, (*O|*O) is split into (*C|*C), (*C|*A), (*A|*A)
!              (V*|V*) is split into (A*|A*), (A*|V*), (V*|V*)
!              (*O|V*) is split into (*A|A*), (*A|V*)
!
!     EXCH(ISYP,ISYI,ISYQ,ISYJ,II,IJ,ERI,SCR)
!
!     rhs_mp2_help1.f
!
      !! The amplitude is in the IC (internally contracted) basis, so
      !! the active orbital index (indices) must be transformed to the
      !! (quasi-)canonical MO (or contravatiant) basis.
      !! IC = SR (why?), contravariant = C
!     write(u6,*) 'icase = ', icase

      PM = .false.
      If (iCase == 2 .or. iCase == 6 .or .iCase == 8 .or.
     &    iCase == 10 .or. iCase == 12) PM = .true.
      If (iCase == 3 .or. iCase == 7 .or .iCase == 9 .or.
     &    iCase == 11 .or. iCase == 13) Return

      SQ2    = SQRT(Two)
      SQI2   = One/SQ2
      SQ3    = SQRT(Three)
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

      !! Some setup
      !! Read T-amplitude, hopefully in contravariant form
!     nINP = 0
!     nINM = 0
!     nIN  = 0
      If (PM) Then
!       nINP = nINDEP(iSym,iCase)
        nASP = nASup(iSym,iCase)
        nISP = nISup(iSym,iCase)
        ! If (nINP /= 0) Then
        !   nVec = nINP*nISP
        ! End If
!       nINM = nINDEP(iSym,iCase+1)
        nASM = nASup(iSym,iCase+1)
        nISM = nISup(iSym,iCase+1)
        ! If (nINM /= 0) Then
        !   nVec = nINM*nISM
        ! End If
        If (nASP*nISP /= 0) Then
          Call RHS_ALLO(nASP,nISP,ipTCP)
          CALL RHS_READ_C(ipTCP,iCase,iSym,iVecC2)
        End If
        If (nASM*nISM /= 0) Then
          Call RHS_ALLO(nASM,nISM,ipTCM)
          CALL RHS_READ_C(ipTCM,iCase+1,iSym,iVecC2)
        End If
      Else
!       nIN = nINDEP(iSym,iCase)
        nAS = nASup(iSym,iCase)
        nIS = nISup(iSym,iCase)
        ! If (nIN /= 0) Then
        !   nVec = nIN*nIS
        ! End If
        If (nAS*nIS /= 0) Then
          Call RHS_ALLO(nAS,nIS,ipTC)
          CALL RHS_READ_C(ipTC,iCase,iSym,iVecC2)
        End If
      End If

      ! finish if no contributions
      ! should we use the number of independent vectors?
      if ((PM .and. nASP*nISP == 0 .and. nASM*nISM == 0)
     &  .or. (.not.PM .and. nAS*nIS == 0)) return

      call mma_allocate(WRK1,nBasT*nBasT,Label='WRK1')
      call mma_allocate(WRK2,nBasT*nBasT,Label='WRK2')

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

      nCorI = nFroI+nIshI
      nCorJ = nFroJ+nIshJ
      nCorA = nFroA+nIshA
      nCorB = nFroB+nIshB
      nOccA = nCorA+nAshA
      nOccB = nCorB+nAshB
      nOccA2= nOccA-nFroA
      nOccB2= nOccB-nFroB
      nOrbA = nOccA+nSshA

      ! nOcc  = nFro(iSym)+nIsh(iSym)+nAsh(iSym)

      !! active+virtual part for the right index
!     nJ = nFro(iSymJ)+nIsh(iSymJ)+nAsh(iSymJ)

      If (iCase == 1) Then
        Call OLagNS_A(Amp1)
      Else If (iCase == 2.or.iCase == 3) Then
        Call OLagNS_B(Amp1)
      Else If (iCase == 4) Then
        Call OLagNS_C(Amp1)
      Else If (iCase == 5) Then
        Call OLagNS_D(Amp1)
      Else If (iCase == 6.or.iCase == 7) Then
        Call OLagNS_E(Amp1)
      Else If (iCase == 8.or.iCase == 9) Then
        Call OLagNS_F(Amp1)
      Else If (iCase == 10.or.iCase == 11) Then
        Call OLagNS_G(Amp1)
      Else If (iCase == 12.or.iCase == 13) Then
        Call OLagNS_H(Amp1)
      End If

      call mma_deallocate(WRK1)
      call mma_deallocate(WRK2)

      If (PM) Then
        If (nASP*nISP /= 0) Call RHS_FREE(ipTCP)
        If (nASM*nISM /= 0) Call RHS_FREE(ipTCM)
      Else
        If (nAS*nIS /= 0) Call RHS_FREE(ipTC)
      End If

      Return

      Contains
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_A(AmpL1)

      implicit none

      real(kind=wp), intent(out) :: AmpL1(nAshA,nAshB)

      if (nAshI*nIshJ*nAshA*nAshB == 0) return

      ! nJ = nIshJ
      Do iI = 1, nAshI
        iIabs = iI + nIshI + nAes(iSymI)
        ! iItot = iI + nCorI
!       If (iSymI == iSymJ) nJ = iI
        Do iJ = 1, nIshJ
          iJabs = iJ + nIes(iSymJ)
          iJtot = iJ + nFroJ
!         If (iIabs < iJabs) Cycle
          Fac = One
!         If ((iI /= iJ).and.(iSymI == iSymJ)) Fac = Two
          If (iSymI /= iSymJ) Fac = Two

          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     &              iI+nCorI,iJ+nFroJ,
     &              ERI1,Scr)

          AmpL1(:,:) = Zero

          Do iA = 1, nAshA
            iAabs = iA + nAes(iSymA)
            ! iAtot = iA + nCorA
            Do iB = 1, nAshB
              iBabs = iB + nAes(iSymB)
              iBtot = iB + nCorB

              iTabs = iBabs
              iUabs = iI + nAes(iSymI)
              iVabs = iAabs
              IW1 = kTUV(iTabs,iUabs,iVabs) - nTUVes(iSym)
!             ValA = Zero
!             Do iICB = 1, nIN
!               iVA  = iICB + nIN*(iJabs-1)
!               ValA  = ValA
!    *        + Work(ipT+iVA-1)*Work(LST+IW1-1+nAS*(iICB-1))
!             End Do
!             ValA = ValA*Two
              iIS = iJabs
              iAS = IW1
              ValA = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two

              If (iUabs == iVabs) Then
                !! For FIMO derivative
                DPT2C(iBtot+nOrbA*(iJtot-1))
     &            = DPT2C(iBtot+nOrbA*(iJtot-1)) + ValA
              End If

              AmpL1(iA,iB) = AmpL1(iA,iB) + ValA
            End Do
          End Do

          AmpL1(:,:) = AmpL1(:,:)*Fac

          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(1,1,ERI1,AmpL1)

          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nAshA,nAshB,nCorA,nCorB,AmpL1,WRK2)
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
        End Do
      End Do

      End Subroutine OLagNS_A
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_B(AmpL1)

      implicit none

      real(kind=wp), intent(out) :: AmpL1(nAshA,nAshB)

      if (nIshI*nIshJ*nAshA*nAshB == 0) return

      nJ = nIshJ
      Do iI = 1, nIshI
        iIabs = iI + nIes(iSymI)
        ! iItot = iI + nFroI
        If (iSymI == iSymJ) nJ = iI
        Do iJ = 1, nJ
          iJabs = iJ + nIes(iSymJ)
          ! iJtot = iJ + nFroJ
          If (iIabs < iJabs) Cycle
          Fac = One
!         If ((iI /= iJ).and.(iSymI == iSymJ)) Fac = Two
          If (iSymI /= iSymJ) Fac = Two

          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     &              iI+nFroI,iJ+nFroJ,
     &              ERI1,Scr)

          AmpL1(:,:) = Zero

          Do iA = 1, nAshA
            iAabs = iA + nAes(iSymA)
            ! iAtot = iA + nCorA
            Do iB = 1, nAshB
              iBabs = iB + nAes(iSymB)
              ! iBtot = iB + nCorB

              if (iaabs > ibabs) then
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
              ValBP = Zero
              ValBM = Zero
!             Do iICB = 1, nINP
!               iVP  = iICB + nINP*(iViP-1)
!               ValBP = ValBP
!    *        + Work(ipTP+iVP-1)*Work(LSTP+iVaP-1+nASP*(iICB-1))
!             End Do
              iIS = iViP
              iAS = iVaP
              ValBP = GA_Arrays(ipTCP)%A(iAS+nASP*(iIS-1))
              If (iAabs /= iBabs .and. iIabs /= iJabs) Then
                If (iIabs /= iJabs) Then
!                 Do iICB = 1, nINM
!                   iVM  = iICB + nINM*(iViM-1)
!                   ValBM = ValBM
!    *          + Work(ipTM+iVM-1)*Work(LSTM+iVaM-1+nASM*(iICB-1))
!                 End Do
                  iIS = iViM
                  iAS = iVaM
                  ValBM = GA_Arrays(ipTCM)%A(iAS+nASM*(iIS-1))
                End If
                !! permutated
                If (iAabs < iBabs) ValBM = -ValBM
              End If
              If (iIabs == iJabs) ValBP = ValBP*SQI2
!
              AmpL1(iA,iB) = AmpL1(iA,iB) + ValBP + ValBM
            End Do
          End Do

          AmpL1(:,:) = AmpL1(:,:)*Fac

          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(1,1,ERI1,AmpL1)

          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nAshA,nAshB,nCorA,nCorB,AmpL1,WRK2)
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
        End Do
      End Do

      End Subroutine OLagNS_B
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_C(AmpL1)

      implicit none

      real(kind=wp), intent(out) :: AmpL1(nAshA+nSshA,nAshB+nSshB)

      if (nAshI*nAshJ*nSshA*nAshB == 0) return

      ! nJ = nIshJ
      Do iI = 1, nAshI
        iIabs = iI + nIshI + nAes(iSymI)
        ! iItot = iI + nCorI
!       If (iSymI == iSymJ) nJ = iI
        Do iJ = 1, nAshJ
          iJabs = iJ + nIshJ + nAes(iSymJ)
          ! iJtot = iJ + nCorJ
          If (iIabs < iJabs) Cycle
          Fac = One
!         If ((iI /= iJ).and.(iSymI == iSymJ)) Fac = Two
          If (iSymI /= iSymJ) Fac = Two

          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     &              iI+nCorI,iJ+nCorJ,
     &              ERI1,Scr)

          AmpL1(:,:) = Zero

          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            iAtot = iA + nFroA + nIshA + nAshA
            Do iB = 1, nAshB
              iBabs = iB + nAes(iSymB)
              iBtot = iB + nCorB

              iTabs = iI + nAes(iSymI)
              iUabs = iBabs
              iVabs = iJ + nAes(iSymJ)
!             write(u6,*) itabs,iuabs,ivabs
              !! (at|uv) -> (ai|bj) -> (at|uv)
              ! IW1 = kTUV(iTabs,iUabs,iVabs) - nTUVes(iSym)
              ValC1 = Zero
              ValC2 = Zero
!             Do iICB = 1, nIN
!               iV  = iICB + nIN*(iAabs-1)
!               ValC1 = ValC1
!    *        + Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
!             End Do
!             If (iIabs /= iJabs) Then
!               IW1 = kTUV(iVabs,iUabs,iTabs) - nTUVes(iSym)
!               Do iICB = 1, nIN
!                 iV  = iICB + nIN*(iAabs-1)
!                 ValC2 = ValC2
!    *          + Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
!               End Do
!             End If
!
!             ValC1 = ValC1*Two
!             ValC2 = ValC2*Two
!
              iIS = iAabs
              iAS = kTUV(iTabs,iUabs,iVabs) - nTUVes(iSym)
              ValC1 = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two
              If (iIabs /= iJabs) Then
                iAS = kTUV(iVabs,iUabs,iTabs) - nTUVes(iSym)
                ValC2 = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two
              End If

              iTabs = iBabs
              iUabs = iI + nAes(iSymI)
              iVabs = iJ + nAes(iSymJ)
              If (iUabs == iVabs) Then
                !! For FIMO derivative
                ONEADD = Zero
!               IW1 = kTUV(iTabs,iUabs,iVabs) - nTUVes(iSym)
!               Do iICB = 1, nIN
!                 iV  = iICB + nIN*(iAabs-1)
!                 ONEADD = ONEADD
!    *          + Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
!               End Do
!               ONEADD = ONEADD*Two
                iAS = kTUV(iTabs,iUabs,iVabs) - nTUVes(iSym)
                ONEADD = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two
                DPT2C(iAtot+nOrbA*(iBtot-1))
     &            = DPT2C(iAtot+nOrbA*(iBtot-1)) + ONEADD

                !! For -sum(y)(ay,yt) -> (ay,ty) derivative
                !! It is correct, but should be rewritten
!               ONEADD = Zero
!               Do iXabs = 1, nAshI !?
!                 IW1 = kTUV(iTabs,iXabs,iXabs) - nTUVes(iSym)
!                 Do iICB = 1, nIN
!                   iV  = iICB + nIN*(iAabs-1)
!                   ONEADD = ONEADD
!    *            + Work(ipT+iV-1)*Work(LST+IW1-1+nAS*(iICB-1))
!                 End Do
!               End Do
!               ONEADD = Two*ONEADD/DBLE(MAX(1,NACTEL))
                ONEADD = Zero
                Do iXabs = 1, nAshI !?
                  iAS = kTUV(iTabs,iXabs,iXabs) - nTUVes(iSym)
                  ONEADD = ONEADD
     &                   + GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))
                End Do
                ONEADD = Two*ONEADD/real(MAX(1,NACTEL),kind=wp)
                AmpL1(iAtot-nCorA,iBtot-nCorA)
     &            = AmpL1(iAtot-nCorA,iBtot-nCorA) - ONEADD
              End If
              AmpL1(iAtot-nCorA,iBtot-nCorB)
     &          = AmpL1(iAtot-nCorA,iBtot-nCorB) + ValC1
              AmpL1(iBtot-nCorB,iAtot-nCorA)
     &          = AmpL1(iBtot-nCorB,iAtot-nCorA) + ValC2
            End Do
          End Do

          AmpL1(:,:) = AmpL1(:,:)*Fac

          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(3,3,ERI1,AmpL1)

          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nAshA+nSshA,nAshB+nSshB,nCorA,nCorB,
     &                      AmpL1,WRK2)
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
        End Do
      End Do

      End Subroutine OLagNS_C
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_D(AmpL1)

      implicit none

      real(kind=wp), intent(out) :: AmpL1(nAshA+nSshA,nAshB+nSshB)

      if (nAshI*nIshJ*nSshA*nAshB == 0) return

      ! nJ = nIshJ
      Do iI = 1, nAshI
        iIabs = iI + nIshI + nAes(iSymI)
        iItot = iI + nCorI
!       If (iSymI == iSymJ) nJ = iI
        Do iJ = 1, nIshJ
          iJabs = iJ + nIes(iSymJ)
          iJtot = iJ + nFroJ
!         If (iIabs < iJabs) Cycle
          Fac = One
!         If ((iI /= iJ).and.(iSymI == iSymJ)) Fac = Two
          If (iSymI /= iSymJ) Fac = Two

          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     &              iI+nCorI,iJ+nFroJ,
     &              ERI1,Scr)

          AmpL1(:,:) = Zero

          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            iAtot = iA + nFroA + nIshA + nAshA
            Do iB = 1, nAshB
              ! iBabs = iB + nAes(iSymB)
              iBtot = iB + nCorB
!
!             iVi   = iJabs + nIshA*(iAabs-1)+IOFF1(iSymA)
!             iVa1  = iB + nAshB*(iI+nAes(iSymI)-1)+IOFF1(iSymB)
!             iVa2  = iB + nAshB*(iI+nAes(iSymI)-1)+IOFF1(iSymB)
!    *                   + nAshT*nAshT
!             !! transform internally contracted (SR)
!             !!        to contravariant (C)
!             ValD1 = Zero
!             ValD2 = Zero
!             Do iICB = 1, nIN
!               iVD   = iICB + nIN*(iVi-1)
!               ValD1 = ValD1
!    *        + Work(ipT+iVD-1)*Work(LST+iVa1-1+nAS*(iICB-1))
!               ValD2 = ValD2
!    *        + Work(ipT+iVD-1)*Work(LST+iVa2-1+nAS*(iICB-1))
!             End Do
!             ValD1 = ValD1*Two
!             ValD2 = ValD2*Two
              iIS = iJabs + nIshA*(iAabs-1)+iOFF1(iSymA)
              iAS = kTU(iB,iI)-nTUes(iSymA)
              ValD1 = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two
              iAS = iAS + nTU(iSymA)
              ValD2 = GA_Arrays(ipTC)%A(iAS+nAS*(iIS-1))*Two

              !! Fock contributions from the inactive density
              If (iItot == iBtot) Then
                DPT2C(iAtot+nOrbA*(iJtot-1))
     &            = DPT2C(iAtot+nOrbA*(iJtot-1)) + ValD1
              End If

              AmpL1(iA+nAshA,iB) = AmpL1(iA+nAshA,iB) + ValD2
              AmpL1(iB,iA+nAshB) = AmpL1(iB,iA+nAshB) + ValD1
            End Do
          End Do

          AmpL1(:,:) = AmpL1(:,:)*Fac

          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(3,3,ERI1,AmpL1)

          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nAshA+nSshA,nAshB+nSshB,nCorA,nCorB,
     &                      AmpL1,WRK2)
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
        End Do
      End Do

      End Subroutine OLagNS_D
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_E(AmpL1)

      implicit none

      real(kind=wp), intent(out) :: AmpL1(nAshA+nSshA,nAshB+nSshB)

      if (nIshI*nIshJ*nSshA*nAshB == 0) return

      nJ = nIshJ
      Do iI = 1, nIshI
        iIabs = iI + nIes(iSymI)
        ! iItot = iI
        If (iSymI == iSymJ) nJ = iI
        Do iJ = 1, nJ
          iJabs = iJ + nIes(iSymJ)
          ! iJtot = iJ
          If (iIabs < iJabs) Cycle
          Fac = One
!         If ((iI /= iJ).and.(iSymI == iSymJ)) Fac = Two
          If (iSymI /= iSymJ) Fac = Two

          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     &              iI+nFroI,iJ+nFroJ,
     &              ERI1,Scr)

          AmpL1(:,:) = Zero

          IgeJ  = kIgeJ(iIabs,iJabs) - nIgeJes(iSym) ! iSymIJ
          IgtJ  = kIgtJ(iIabs,iJabs) - nIgtJes(iSym) ! iSymIJ
          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            ! iAtot = iA + nIshA + nAshA
            Do iB = 1, nAshB
              iBabs = iB + nAes(iSymB)
              ! iBtot = iB + nIshB

              iASP  = iBabs
              iISP  = iAabs + nSshA*(IgeJ-1)+iOFF1(iSymA)
              ValEP = GA_Arrays(ipTCP)%A(iASP+nASP*(iISP-1))
              ValEM = Zero
              If (iIabs > iJabs) Then
                ValEP = ValEP * SQ2
                iASM  = iBabs
                iISM  = iAabs + nSshA*(IgtJ-1)+iOFF1(iSymA)
                ValEM = GA_Arrays(ipTCM)%A(iASM+nASM*(iISM-1))
     &                *SQ2*SQ3
              Else
              End If

              AmpL1(iA+nAshA,iB) = AmpL1(iA+nAshA,iB) + ValEP + ValEM
              AmpL1(iB,iA+nAshB) = AmpL1(iB,iA+nAshB) + ValEP - ValEM
            End Do
          End Do

          AmpL1(:,:) = AmpL1(:,:)*Fac

          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(3,3,ERI1,AmpL1)

          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nAshA+nSshA,nAshB+nSshB,nCorA,nCorB,
     &                      AmpL1,WRK2)
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
        End Do
      End Do

      End Subroutine OLagNS_E
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_F(AmpL1)

      implicit none

      real(kind=wp), intent(out) :: AmpL1(nSshA,nSshB)

      if (nAshI*nAshJ*nSshA*nSshB == 0) return

      ! nJ = nIshJ
      Do iI = 1, nAshI
        iIabs = iI + nIshI + nAes(iSymI)
        ! iItot = iI + nCorI
!       If (iSymI == iSymJ) nJ = iI
        Do iJ = 1, nAshI
          iJabs = iJ + nIshJ + nIes(iSymJ)
          ! iJtot = iJ + nCorJ
          If (iIabs < iJabs) Cycle
          Fac = One
!         If ((iI /= iJ).and.(iSymI == iSymJ)) Fac = Two
          If (iSymI /= iSymJ) Fac = Two

          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     &              iI+nCorI,iJ+nCorJ,
     &              ERI1,Scr)

          AmpL1(:,:) = Zero

          iTabs = iI + nAes(iSymI)
          iUabs = iJ + nAes(iSymJ)
          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            ! iAtot = iA + nIsh(iSymA) + nAsh(iSymA)
            Do iB = 1, nSshB
              iBabs = iB + nSes(iSymB)
              If (iAabs < iBabs) Cycle
              ! iBtot = iB + nIsh(iSymB) + nAsh(iSymB)

              iASP  = kTgeU(iTabs,iUabs)-nTgeUes(iSym)
              iISP  = kAgeB(iAabs,iBabs)-nAgeBes(iSym)
              ValFP = GA_Arrays(ipTCP)%A(iASP+nASP*(iISP-1))
              If (iIabs == iJabs) ValFP = ValFP*Half
              ValFM = Zero
              If (iAabs /= iBabs) Then
                If (iTabs /= iUabs) Then
                  iASM  = kTgtU(iTabs,iUabs)-nTgtUes(iSym)
                  iISM  = kAgtB(iAabs,iBabs)-nAgtBes(iSym)
                  ValFM = GA_Arrays(ipTCM)%A(iASM+nASM*(iISM-1))
                End If
              Else
                ValFP = ValFP * SQI2
              End If
              VALFM = -VALFM !! why?

              AmpL1(iA,iB) = AmpL1(iA,iB) + ValFP + ValFM
              AmpL1(iB,iA) = AmpL1(iB,iA) + ValFP - ValFM
            End Do
          End Do

          AmpL1(:,:) = AmpL1(:,:)*Fac

          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(2,2,ERI1,AmpL1)

          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu}
          Call OLagNS_post2(nSshA,nSshB,nOccA,nOccB,AmpL1,WRK2)
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
        End Do
      End Do

      End Subroutine OLagNS_F
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_G(AmpL1)

      implicit none

      real(kind=wp), intent(out) :: AmpL1(nSshA,nSshB)

      if (nAshI*nIshJ*nSshA*nSshB == 0) return

      ! nJ = nIshJ
      Do iI = 1, nAshI
        iIabs = iI + nIshI + nAes(iSymI)
        ! iItot = iI + nCorI
!       If (iSymI == iSymJ) nJ = iI
        Do iJ = 1, nIshI
          iJabs = iJ + nIes(iSymJ)
          ! iJtot = iJ
!         If (iIabs < iJabs) Cycle
          Fac = One
!         If ((iI /= iJ).and.(iSymI == iSymJ)) Fac = Two
          If (iSymI /= iSymJ) Fac = Two

          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     &              iI+nCorI,iJ+nFroJ,
     &              ERI1,Scr)

          AmpL1(:,:) = Zero

          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            ! iAtot = iA + nIsh(iSymA) + nAsh(iSymA)
            Do iB = 1, nSshB
              iBabs = iB + nSes(iSymB)
              If (iAabs < iBabs) Cycle
              ! iBtot = iB + nIsh(iSymB) + nAsh(iSymB)

              iAgeB = kAgeB(iAabs,iBabs)-nAgeBes(iSym) !! iSymAB
              iVjP  = iJ + nIsh(iSymJ)*(iAgeB-1)+IOFF1(iSymJ)
              ValGP = GA_Arrays(ipTCP)%A(iI+nASP*(iVjP-1))
              ValGM = Zero
              If (iAabs /= iBabs) Then
                ValGP = ValGP * SQ2
                iAgtB = kAgtB(iAabs,iBabs) - nAgtBes(iSym) !! iSymAB
                iVjM  = iJ + nIsh(iSymJ)*(iAgtB-1)+IOFF2(iSymJ)
                ValGM = GA_Arrays(ipTCM)%A(iI+nASM*(iVjM-1))*SQ2*SQ3
              End If

              AmpL1(iA,iB) = AmpL1(iA,iB) + ValGP + ValGM
              AmpL1(iB,iA) = AmpL1(iB,iA) + ValGP - ValGM
            End Do
          End Do

          AmpL1(:,:) = AmpL1(:,:)*Fac

          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(2,2,ERI1,AmpL1)
!         Call DGEMM_('N','T',nOrbA,nSshA,nSshB,
!    *                One,ERI1(1+nOrbA*nOccB),nOrbA,
!    *                        AmpL1,nSshA,
!    *                One,OLAG(nOrbA*nOccB+1),nOrbA)
!         Call DGEMM_('T','N',nOrbA,nSshA,nSshB,
!    *                One,ERI1(nOccA+1),nOrbA,
!    *                        AmpL1,nSshA,
!    *                One,OLAG(nOrbA*nOccB+1),nOrbA)

          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu} for all ij
          Call OLagNS_post2(nSshA,nSshB,nOccA,nOccB,AmpL1,WRK2)
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
        End Do
      End Do

      End Subroutine OLagNS_G
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_H(AmpL1)

      implicit none

      real(kind=wp), intent(out) :: AmpL1(nSshA,nSshB)

      if (nIshI*nIshJ*nSshA*nSshB == 0) return

      nJ = nIshJ
      Do iI = 1, nIshI
        iIabs = iI + nIes(iSymI)
        ! iItot = iI
        If (iSymI == iSymJ) nJ = iI
        Do iJ = 1, nJ
          iJabs = iJ + nIes(iSymJ)
          ! iJtot = iJ
          If (iIabs < iJabs) Cycle
          Fac = One
!         If ((iI /= iJ).and.(iSymI == iSymJ)) Fac = Two
          If (iSymI /= iSymJ) Fac = Two

          Call Exch(iSymA,iSymI,iSymB,iSymJ,
     &              iI+nFroI,iJ+nFroJ,
     &              ERI1,Scr)

          AmpL1(:,:) = Zero

          iViHP0= kIgeJ(iIabs,iJabs) - nIgeJes(iSym)
          iViHP = nAgeB(iSym)*(iViHP0-1)
          iViHM0= kIgtJ(iIabs,iJabs) - nIgtJes(iSym)
          iViHM = nAgtB(iSym)*(iViHM0-1)
          Do iA = 1, nSshA
            iAabs = iA + nSes(iSymA)
            ! iAtot = iA + nIsh(iSymA) + nAsh(iSymA)
            Do iB = 1, nSshB
              iBabs = iB + nSes(iSymB)
              If (iAabs < iBabs) Cycle
              ! iBtot = iB + nIsh(iSymB) + nAsh(iSymB)
              iVaHP = kAgeB(iAabs,iBabs) - nAgeBes(iSym)
              iVHP  = iVaHP + iViHP !! nAgeB(iSym)*(iViP-1)

              ValHP = GA_Arrays(ipTCP)%A(iVHP)
              ValHM = Zero
              If (iIabs /= iJabs) Then
                If (iAabs /= iBabs) Then
                  ValHP = ValHP * Two
                  iVaHM = kAgtB(iAabs,iBabs) - nAgtBes(iSym)
                  iVHM  = iVaHM + iViHM !! nAgtB(iSym)*(iViM-1)
                  ValHM = GA_Arrays(ipTCM)%A(iVHM) * Two*SQ3
                Else
                  ValHP = ValHP * SQ2
                End If
              Else
                If (iAabs /= iBabs) ValHP = ValHP * SQ2
              End If

!     write(u6,'(2i3,2f20.10)')ia,ib,valhp,valhm
              AmpL1(iA,iB) = AmpL1(iA,iB) + ValHP + ValHM
              AmpL1(iB,iA) = AmpL1(iB,iA) + ValHP - ValHM
            End Do
          End Do

          AmpL1(:,:) = AmpL1(:,:)*Fac

          !! Calculate the actual contributions
          !! L_{pq} = sum_{j,ab} (pa|jb) * T_{qj}^{ab}
          Call OLagNS_post1(2,2,ERI1,AmpL1)
!         Call DGEMM_('N','T',nOrbA,nSshA,nSshB,
!    *                One,ERI1(1+nOrbA*nOccB),nOrbA,
!    *                        AmpL1,nSshA,
!    *                One,OLAG(nOrbA*nOccB+1),nOrbA)
!         Call DGEMM_('T','N',nOrbA,nSshA,nSshB,
!    *                One,ERI1(nOccA+1),nOrbA,
!    *                        AmpL1,nSshA,
!    *                One,OLAG(nOrbA*nOccB+1),nOrbA)
!
          !! Prepare for implicit (VV|VO) integrals
          !! T_{ij}^{ab} -> T_{ij}^{mu nu} back-transformation
          !! 1) T_{ij}^{ab} -> T_{ij}^{mu nu} for all ij
          Call OLagNS_post2(nSshA,nSshB,nOccA,nOccB,AmpL1,WRK2)
!         Call DGEMM_('N','N',nBasT,nSshA,nSshB,
!    *                One,CMOPT2(1+nBasT*nOccB),nBasT,
!    *                        AmpL1,nSshA,
!    *                Zero,WRK1,nBasT)
!         Call DGEMM_('N','T',nBasT,nBasT,nSshB,
!    *                One,WRK1,nBasT,
!    *                        CMOPT2(1+nBasT*nOccA),nBasT,
!    *                Zero,WRK2,nBasT)
          !! Reorder T_{ij}^{rho sigma} to T2AO(j,sigma,i,rho)
          Call OLagNS_post3(iIabs,iJabs,T2AO,WRK2)
!         Do iBas = 1, nBasT
!           Do jBas = 1, nBasT
!             loc1 = iJ-1 + (jBas-1)*nOccA2
!    *             + (iI-1)*nOccA2*nBasT + (iBas-1)*nOccA2*nBasT*nOccA2
!             loc2 = iI-1 + (jBas-1)*nOccA2
!    *             + (iJ-1)*nOccA2*nBasT + (iBas-1)*nOccA2*nBasT*nOccA2
!             loc3 = iBas-1 + (jBas-1)*nBasT
!             loc4 = jBas-1 + (iBas-1)*nBasT
!             T2AO(1+loc1) = T2AO(1+loc1) + WRK2(1+loc3)
!             T2AO(1+loc2) = T2AO(1+loc2) + WRK2(1+loc4)
!           End Do
!         End Do
        End Do
      End Do

      End Subroutine OLagNS_H
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_post1(iLeft,iRight,ERI,AmpMO)

      implicit none

      integer(kind=iwp), intent(in) :: iLeft, iRight
      real(kind=wp), intent(in) :: ERI(*), AmpMO(*)

      integer(kind=iwp) :: nSkpA, nSkpB, nDimA, nDimB

      nSkpA = 0
      nSkpB = 0
      nDimA = 0
      nDimB = 0
      If (iLeft == 1) Then
        nDimA = nAshA
        nSkpA = nCorA
      Else If (iLeft == 2) Then
        nDimA = nSshA
        nSkpA = nOccA
      Else If (iLeft == 3) Then
        nDimA = nAshA+nSshA
        nSkpA = nCorA
      End If
      If (iRight == 1) Then
        nDimB = nAshB
        nSkpB = nCorB
      Else If (iRight == 2) Then
        nDimB = nSshB
        nSkpB = nOccB
      Else If (iRight == 3) Then
        nDimB = nAshB+nSshB
        nSkpB = nCorB
      End If

      Call DGEMM_('N','T',nOrbA,nDimA,nDimB,
     &            One,ERI(1+nOrbA*nSkpB),nOrbA,AmpMO,nDimA,
     &            One,OLAG(nOrbA*nSkpB+1),nOrbA)
      Call DGEMM_('T','N',nOrbA,nDimA,nDimB,
     &            One,ERI(nSkpA+1),nOrbA,AmpMO,nDimA,
     &            One,OLAG(nOrbA*nSkpB+1),nOrbA)

      End Subroutine OLagNS_post1
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_post2(nDimA,nDimB,nSkpA,nSkpB,AmpMO,AmpAO)

      use caspt2_global, only: CMOPT2

      implicit none

      integer(kind=iwp), intent(in) :: nDimA, nDimB, nSkpA, nSkpB
      real(kind=wp), intent(in) :: AmpMO(nDimA,nDimB)
      real(kind=wp), intent(out) :: AmpAO(nBasA,nBasB)

       Call DGEMM_('N','N',nBasA,nDimB,nDimA,
     &             One,CMOPT2(1+nBasA*nSkpA),nBasA,AmpMO,nDimA,
     &             Zero,WRK1,nBasA)
       Call DGEMM_('N','T',nBasA,nBasB,nDimB,
     &             One,WRK1,nBasA,CMOPT2(1+nBasB*nSkpB),nBasA,
     &             Zero,AmpAO,nBasA)

      End Subroutine OLagNS_post2
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_post3(iIabs,iJabs,TampAO,TampIJ)

      implicit none

      integer(kind=iwp), intent(in) :: iIabs, iJabs
      real(kind=wp), intent(inout) :: TampAO(nOccA2,nBasI,nOccB2,nBasJ)
      real(kind=wp), intent(in) :: TampIJ(nBasI,nBasJ)

      integer(kind=iwp) :: iBas, jBas

      Do iBas = 1, nBasI
        Do jBas = 1, nBasJ
          TampAO(iJabs,jBas,iIabs,iBas)
     &      = TampAO(iJabs,jBas,iIabs,iBas) + TampIJ(iBas,jBas)
          TampAO(iIabs,jBas,iJabs,iBas)
     &      = TampAO(iIabs,jBas,iJabs,iBas) + TampIJ(jBas,iBas)
        End Do
      End Do

      End Subroutine OLagNS_post3

      END SUBROUTINE OLagNS_Hel2
!
!-----------------------------------------------------------------------
!
! MO->AO or AO->MO transformation of 1-RDM
      Subroutine OLagTrf(mode,iSym,CMO,DPT2,DPT2AO,WRK)

      use caspt2_module, only: NFRO, NORB, NDEL, NBAS
      use Constants, only: Zero, One, Half
      use definitions, only: wp, iwp

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: mode, iSym
      real(kind=wp), intent(in) :: CMO(*)
      real(kind=wp), intent(inout) :: DPT2(*), DPT2AO(*)
      real(kind=wp), intent(_OUT_) :: WRK(*)

      real(kind=wp) :: Val
      integer(kind=iwp) :: iCMO, iAO, iMO, jSym, nBasI, nOrbI, iBas,
     &  jBas

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

      If (nOrb(iSym)+nFro(iSym) > 0) Then
        nBasI = nBas(iSym)
        nOrbI = nBas(iSym)-nDel(iSym)
        If (Mode == 1) Then
          !! MO -> AO
          CALL DGEMM_('N','N',nBasI,nOrbI,nOrbI,
     &                One,CMO(iCMO),nBasI,DPT2(iMO),nOrbI,
     &                Zero,WRK,nBasI)
          CALL DGEMM_('N','T',nBasI,nBasI,nOrbI,
     &                One,WRK,nBasI,CMO(iCMO),nBasI,
     &                Zero,DPT2AO(iAO),nBasI)
          !! Symmetrize, just in case
          Do iBas = 1, nBasI
            Do jBas = 1, iBas-1
              Val =(DPT2AO(iAO+iBas-1+nBasI*(jBas-1))
     &            + DPT2AO(iAO+jBas-1+nBasI*(iBas-1)))*Half
              DPT2AO(iAO+iBas-1+nBasI*(jBas-1)) = Val
              DPT2AO(iAO+jBas-1+nBasI*(iBas-1)) = Val
            End Do
          End Do
        Else If (Mode == 2) Then
          !! AO -> MO
          CALL DGEMM_('T','N',nOrbI,nBasI,nBasI,
     &                One,CMO(iCMO),nBasI,DPT2AO(iAO),nBasI,
     &                Zero,WRK,nOrbI)
          CALL DGEMM_('N','N',nOrbI,nOrbI,nBasI,
     &                One,WRK,nOrbI,CMO(iCMO),nBasI,
     &                Zero,DPT2(iMO),nOrbI)
        End If
      END IF

      Return

      End Subroutine OLagTrf
