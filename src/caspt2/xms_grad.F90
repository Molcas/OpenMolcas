!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

      Subroutine XMS_Grad(H0,U0,UEFF,OMGDER)

      use caspt2_global, only: do_nac, do_csf, iRoot1, iRoot2,          &
     &                         nOLag, CLag,CLagFull,OLag,DPT2_tot,      &
     &                         FIFA_all,FIFASA_all
      use caspt2_global, only: FIFA, TORB, NDREF
      use caspt2_global, only: CMOPT2, if_equalW, weight
      use sguga, only: SGS
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp, iwp
      use caspt2_module, only: ENERGY, IFXMS, IFRMS, IFDW, STSYM, NCONF,&
     &                         NFRO, NISH, NASHT, NDEL, NBAS, NBAST,    &
     &                         NBSQT, NROOTS, NSTATE, ZETA, ORBIN, MXCI
      use Constants, only: Zero, One, Half, Two, Quart
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif

      implicit none

      real(kind=wp), intent(in) :: H0(nState,nState), U0(nState,nState),&
     &                             OMGDER(nState,nState)
      real(kind=wp), intent(inout) :: UEFF(nState,nState)

      real(kind=wp),allocatable :: CI1(:),CI2(:),SGM1(:),SGM2(:),TG1(:),&
     &                             TG2(:),DG1(:),DG2(:),DG3(:),G1(:,:), &
     &                             RDMEIG(:,:),DPT2(:),Trf(:),          &
     &                             RDMSA(:,:),WRK1(:),WRK2(:),DPT2_AO(:)

      real(kind=wp) :: SLag(nState*nState), Scal, OVL, EEI, EEJ, fact,  &
     &  Wgt, TRC, EINACT, EDIFF
      integer(kind=iwp) :: nLev, iStat, jStat, NTG1, NTG3, iState,      &
     &  kStat, lStat, nOrbI, nBasI, iAsh, jAsh, nCor, I, II
      real(kind=wp), external :: ddot_

      nLev = SGS%nLev
!
!     The XMS rotation applies to any variants: XMS-CASPT2, XDW-CASPT2,
!     and RMS-CASPT2.
!
      If (IFXMS .or. IFRMS) Then
        call mma_allocate(CI1,nConf,Label='CI1')
        call mma_allocate(CI2,nConf,Label='CI2')
        call mma_allocate(SGM1,MXCI,Label='SGM1')
        call mma_allocate(SGM2,MXCI,Label='SGM2')
        call mma_allocate(TG1,nAshT**2,Label='TG1')
        call mma_allocate(TG2,nAshT**4,Label='TG2')

        call mma_allocate(DG1,nAshT**2,Label='DG1')
        call mma_allocate(DG2,nAshT**4,Label='DG2')

!       ----- Construct pseudo-density matrix -----

        !! First, we need to consider the derivative of the rotation
        !! of the CASSCF energy (Heff[1]).

        !! Forward transformation of UEFF (CASSCF basis to XMS basis)
        !! SLag is used as a working array
        Call DGEMM_('T','N',nState,nState,nState,                       &
     &              One,U0,nState,UEFF,nState,                          &
     &              Zero,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)

        !! Then actual calculation
        CLag(1:nconf,1:nstate) = Zero
        Do iStat = 1, nState
          Call LoadCI_XMS('C',0,nConf,nState,CI1,iStat,U0)
          Do jStat = 1, nState
            Call LoadCI_XMS('C',0,nConf,nState,CI2,jStat,U0)

            Call Dens2T_RPT2(NLEV,NCONF,MXCI,CI1,CI2,SGM1,SGM2,TG1,TG2)
            TG1(:) = TG1(:)*Half
            TG2(:) = TG2(:)*Half
            DG1(:) = Zero
            DG2(:) = Zero
            Call CnstInt(1,DG1,DG2)

            If (do_nac) Then
              Scal =(UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)              &
     &             + UEFF(jStat,iRoot1)*UEFF(iStat,iRoot2))*Half
            Else
              Scal = UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)
            End If
            if (IFDW .and. zeta >= Zero) then
              scal = scal + OMGDER(iStat,jStat)
            end if
            DG1(:) = DG1(:)*Scal
            DG2(:) = DG2(:)*Scal

            call mma_allocate(DG3,nAshT**6,Label='DG3')
            DG3(:) = Zero

            NTG1 = NASHT**2
            NTG3 = (NTG1*(NTG1+1)*(NTG1+2))/6
            OVL = Zero
            CALL DERTG3(.False.,STSYM,STSYM,NCONF,NASHT,CI1,CI2,OVL,    &
     &                  DG1,DG2,NTG3,DG3,CLag(1,iStat),CLag(1,jStat))
            call mma_deallocate(DG3)
          End Do
        End Do
        call mma_deallocate(SGM2)
        call mma_deallocate(TG2)

        !! Back transformation of UEFF (XMS basis to CASSCF basis)
        !! SLag is used as a working array
        Call DGEMM_('N','N',nState,nState,nState,                       &
     &              One,U0,nState,UEFF,nState,                          &
     &              Zero,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)
        call mma_deallocate(DG1)
        call mma_deallocate(DG2)

      !! Add to the full CI Lagrangian
      !! CLag is in quasi-canonical basis, so transformation
      !! to natural basis is required.
        IF(ORBIN == 'TRANSFOR') Then
          Do iState = 1, nState
            Call CLagX_TrfCI(nConf,CLag(1,iState))
          End Do
        End If
        CLagFull(:,:) = CLagFull(:,:) + CLag(:,:)

        !! Compute the Lagrange multiplier for XMS
        !! The diagonal element is always zero.
        !! The code has an additional scaling with 0.5,
        !! because some contributions are doubled.
        SLag(:) = Zero
        Do iStat = 1, nState
          Call LoadCI_XMS('N',0,nConf,nState,CI1,iStat,U0)
          EEI = H0(iStat,iStat)
          Do jStat = 1, iStat
            If (iStat == jStat) Then
              SLag(iStat+nState*(jStat-1)) = Zero
            Else
              Call LoadCI_XMS('N',0,nConf,nState,CI2,jStat,U0)
              EEJ = H0(jStat,jStat)

              Scal = DDOT_(nConf,CI1,1,CLagFull(1,jStat),1)             &
     &             - DDOT_(nConf,CI2,1,CLagFull(1,iStat),1)
              If (do_csf) Then
                !! JCTC 2017, 13, 2561: eq.(66)
                !! iStat and jStat: XMS
                !! kStat and lStat: CASSCF
                fact = Zero
                Do kStat = 1, nState
                  Do lStat = 1, nState
                    fact = fact                                         &
     &                   + (UEFF(kStat,iRoot1)*UEFF(lStat,iRoot2)       &
     &                   -  UEFF(kStat,iRoot2)*UEFF(lStat,iRoot1))*Half &
     &                   * U0(kStat,iStat)*U0(lStat,jStat)
                  End Do
                End Do
                Scal = Scal+fact*(ENERGY(iRoot2)-ENERGY(iRoot1))*Two
              End If
              Scal = Quart*Scal/(EEJ-EEI)
              SLag(iStat+nState*(jStat-1)) = Scal
              SLag(jStat+nState*(iStat-1)) = Scal
            End If
          End Do
        End Do

      !! Either subtract the CI derivative of Heff(1)
      !! or add off-diagonal couplings in rhs_sa (Z-vector)
      !CLagFull = CLagFull - CLag

        !! Finally, construct the pseudo-density matrix
        !! d = \sum_{ST} w_{ST} * d_{ST}
        call mma_allocate(G1,nAshT,nAshT,Label='G1')
        G1(:,:) = Zero
        Do iStat = 1, nState
          Call LoadCI_XMS('C',0,nConf,nState,CI1,iStat,U0)
          Do jStat = 1, nState
            If (ABS(SLag(iStat+nState*(jStat-1))) <= 1.0e-10_wp) Cycle
            Call LoadCI_XMS('C',0,nConf,nState,CI2,jStat,U0)

            Call Dens1T_RPT2(CI1,CI2,                                   &
     &                       SGM1,TG1,nLev)
            Scal = SLag(iStat+nState*(jStat-1))*Two
            Call DaXpY_(nAshT**2,Scal,TG1,1,G1,1)
          End Do
        End Do

        call mma_deallocate(CI1)
        call mma_deallocate(CI2)
        call mma_deallocate(SGM1)
        call mma_deallocate(TG1)
!
!       ----- Calculate orbital derivatives -----
!
        call mma_allocate(RDMEIG,nAshT,nAshT,Label='RDMEIG')
        call mma_allocate(DPT2,NBSQT,Label='DPT2')
        DPT2(:) = Zero

        call mma_allocate(Trf,NBSQT,Label='TRFMAT')
        call mma_allocate(RDMSA,nAshT,nAshT,Label='RDMSA')
        call mma_allocate(WRK1,NBSQT,Label='WRK1')
        call mma_allocate(WRK2,NBSQT,Label='WRK2')

        !! Construct always state-averaged density; XMS basis is always
        !! generated with the equally-averaged density.
        WRK1(1:nDRef) = Zero
        call mma_allocate(CI1,nConf,Label='CI1')
        Wgt  = One/nState
        Do iState = 1, nState
          Call LoadCI(CI1,iState)
          call POLY1(CI1,nConf)
          call GETDREF(WRK2,nDRef)
          WRK1(1:nDRef) = WRK1(1:nDRef) + Wgt*WRK2(1:nDRef)
        End Do
        call mma_deallocate(CI1)
        Call SQUARE(WRK1,RDMSA,1,nAshT,nAshT)

        nOrbI = nBas(1) - nDel(1) !! nOrb(1)
        nBasI = nBas(1)
!       Call SQUARE(FIFA,FIFA_all,1,nOrbI,nOrbI)
        !! FIFASA_all is in natural orbital basis
        Trf(1:NBSQT) = Zero
        Call CnstTrf(NBSQT,TOrb,Trf)

        !! FIFA: natural -> quasi-canonical
        If (IFDW .or. IFRMS) Then
          Call DGemm_('T','N',nOrbI,nOrbI,nOrbI,                        &
     &                One,Trf,nBasI,FIFASA_all,nOrbI,                   &
     &                Zero,WRK1,nOrbI)
          Call DGemm_('N','N',nOrbI,nOrbI,nOrbI,                        &
     &                One,WRK1,nOrbI,Trf,nBasI,                         &
     &                Zero,FIFASA_all,nOrbI)
        End If
!       Call DCopy_(NBSQT,FIFASA_all,1,FIFA_all,1)
        FIFA_all(1:NBSQT) = FIFASA_all(1:NBSQT)

        !! Orbital derivatives of FIFA
        !! Both explicit and implicit orbital derivatives are computed
        !! Also, compute G(D) and put the active contribution to RDMEIG
        OLag(:) = Zero
        Call EigDer2(NBSQT,nAshT,RDMEIG,Trf,FIFA_all,RDMSA,G1,WRK1,WRK2)

        !! Add to PT2 density
        !! No inactive contributions. Correct as long as CASSCF CI
        !! vector are orthogonal.
        Call AddDEPSA(NBSQT,nAshT,DPT2,G1)
        Call DPT2_TrfStore(One,NBSQT,DPT2,DPT2_tot,Trf,WRK1)

        !! Finalize OLag (anti-symetrize) and construct WLag
        Call OLagFinal(nOLag,NBSQT,OLag,Trf)

        If (.not.if_equalW) then
          call mma_allocate(DPT2_AO,NBSQT,Label='DPT2_AO')

          !! AddDEPSA considers the frozen orbital, whereas DPT2_Trf
          !! does not. In any case, construct DPT2 again.
          DPT2(:) = Zero
          CALL DPT2_Trf(NBSQT,nAshT,DPT2,DPT2_AO,CMOPT2,G1,WRK1)
          !! Construct the SCF density
          WRK1(1:nDRef) = Zero
          call mma_allocate(CI1,nConf,Label='CI1')
          Do iState = 1, nState
            Call LoadCI_XMS('N',1,nConf,nState,CI1,iState,U0)
            call POLY1(CI1,nConf)
            call GETDREF(WRK2,nDRef)
            wgt = Weight(iState)
            WRK1(1:nDRef) = WRK1(1:nDRef) + Wgt*WRK2(1:nDRef)
          End Do
          call mma_deallocate(CI1)
          !! WRK2 is the SCF density (for nstate=nroots)
          Call SQUARE(WRK1,WRK2,1,nAshT,nAshT)
          Call DaXpY_(nAshT**2,-One,WRK2,1,RDMSA,1)
          !! Construct the SS minus SA density matrix in WRK1
          Call OLagFroD(NBSQT,nAshT,WRK1,WRK2,RDMSA,Trf)
          !! Subtract the inactive part
          Call DaXpY_(nBasT**2,-One,WRK2,1,WRK1,1)
          !! Save
          Call CnstAB_SSDM(NBSQT,DPT2_AO,WRK1)
          call mma_deallocate(DPT2_AO)
        End If

        call mma_deallocate(RDMSA)
        call mma_deallocate(WRK1)
        call mma_deallocate(WRK2)

        call mma_deallocate(DPT2)
!
!       ----- Calculate CI derivatives -----
!
        !! use quasi-canonical CSF rather than natural CSF
        CLag(:,:) = Zero

        !! 1) Explicit CI derivative
        !! a: Extract FIFA in the AS for explicit CI derivative
        !!    Note that this FIFA uses state-averaged density matrix
!       call sqprt(fifa_all,nbast)
        Do iAsh = 1, nAshT
          Do jAsh = 1, nAshT
            G1(iAsh,jAsh) = FIFA_all( nFro(1)                           &
     &          + nIsh(1) + iAsh + nBas(1)*(nFro(1)+nIsh(1)+jAsh-1))
          End Do
        End Do
!       call sqprt(g1),nasht)
        !! Transform quasi-canonical to natural
        nCor = nFro(1)+nIsh(1)
!     call sqprt(trf,nbast)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,                          &
     &              One,Trf(1+nBasT*nCor+nCor),nBasT,G1,nAshT,          &
     &              Zero,FIFA_all,nAshT)
!     call sqprt(FIFA_all,nasht)
        Call DGemm_('N','T',nAshT,nAshT,nAshT,                          &
     &              One,FIFA_all,nAshT,Trf(1+nBasT*nCor+nCor),nBasT,    &
     &              Zero,G1,nAshT)
!     call sqprt(g1,nasht)
        Call DGemm_('N','N',nAshT,nAshT,nAshT,                          &
     &              One,Trf(1+nBasT*nCor+nCor),nBasT,RDMEIG,nAshT,      &
     &              Zero,FIFA_all,nAshT)
        Call DGemm_('N','T',nAshT,nAshT,nAshT,                          &
     &              One,FIFA_all,nAshT,Trf(1+nBasT*nCor+nCor),nBasT,    &
     &              Zero,RDMEIG,nAshT)
        call mma_deallocate(Trf)

        !! b: FIFA in the inactive space
        TRC=Zero
!     DO ISYM=1,NSYM
        DO I=1,NISH(1) ! ISYM)
!         II=IOFF(ISYM)+(I*(I+1))/2
          II=(I*(I+1))/2
          TRC=TRC+FIFA(II)
        END DO
!     END DO
! Contribution from inactive orbitals:
        EINACT=Two*TRC
        !! This EINACT may be wrong. Perhaps, FIFA has to be
        !! back-transformed to natural orbital basis. However, this does
        !! not contribute to the final gradient as long as all the
        !! (internal) CI vectors are orthogonal.

        !! c: Finally, compute explicit CI derivative
        !! y_{I,T} = w_{ST} <I|f|S>
        !! Here, G1 is FIFA = ftu
        Call CLagEigT(CLag,G1,SLag,EINACT)

        !! 2) Implicit CI derivative
        Call CLagEig(.False.,.True.,nConf,nRoots,nState,nLev,CLag,      &
     &               RDMEIG)

#ifdef _MOLCAS_MPP_
        if (is_real_par()) then
          call GADGOP(CLag,nConf*nState,'+')
        end if
#endif

        call mma_deallocate(RDMEIG)
        call mma_deallocate(G1)

      End If

      If (do_csf) Then
        !! Eq (68)
        !! I think this contributes even in XMS-CASPT2.
        !! However, maybe I implemented other terms in a different way
        !! from BAGEL.
        call mma_allocate(CI1,nConf,Label='CI1')
        !! Do in the natural orbital in the XMS basis
        !! because CLag is like that
        !! CASSCF -> XMS
        Call DGEMM_('T','N',nState,nState,nState,                       &
     &              One,U0,nState,UEFF,nState,                          &
     &              Zero,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)

        EDIFF = ENERGY(iRoot2)-ENERGY(iRoot1)
        Do iStat = 1, nState
          Call LoadCI_XMS('N',0,nConf,nState,CI1,iStat,U0)
          Do jStat = 1, nState
            Scal = (UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)               &
     &           -  UEFF(iStat,iRoot2)*UEFF(jStat,iRoot1))*Half
            Scal = Scal * EDIFF
            Call DaXpY_(nConf,Scal,CI1,1,CLag(1,jStat),1)
          End Do
        End Do

        !! XMS -> CASSCF
        Call DGEMM_('N','N',nState,nState,nState,                       &
     &              One,U0,nState,UEFF,nState,                          &
     &              Zero,SLag,nState)
        Call DCopy_(nState**2,SLag,1,UEFF,1)

        call mma_deallocate(CI1)
      End If

      !! Finally, add to the total CI derivative array
      CLagFull(:,:) = CLagFull(:,:) + CLag(:,:)

      End Subroutine XMS_Grad
