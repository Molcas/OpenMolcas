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
      SUBROUTINE VVVOX(NSYM,NBAS,NFRO,KEEP,                             &
     &                 iSymI,iSymJ,iSymK,iSymL,NBMX,                    &
     &                 T2AO,vLag,CMO,nOcc,nBasT,LBUF,X1,X2,WRK,         &
     &                 DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,                   &
     &                 DIA,DI,FIFA,FIMO)

      use definitions, only: wp, iwp, u6
      use Constants, only: Zero, One, Half

      implicit none

      integer(kind=iwp), intent(in) :: NSYM, NBAS(8), NFRO(8), KEEP(8), &
     &  iSymI, iSymJ, iSymK, iSymL, NBMX, nOcc, nBasT, LBUF
      real(kind=wp), intent(in) :: T2AO(nOcc,nBasT,nOcc,nBasT),         &
     &  CMO(nBasT,nBasT), DPT2AO(nBasT**2), DPT2CAO(nBasT**2),          &
     &  DIA(nBasT**2), DI(nBasT**2)
      real(kind=wp), intent(inout) :: vLag(nBasT,nBasT), X1(LBUF),      &
     &  X2(NBMX*NBMX), WRK(nBasT**2), FPT2AO(nBasT**2),                 &
     &  FPT2CAO(nBasT**2), FIFA(nBasT**2), FIMO(nBasT**2)

      integer(kind=iwp) ::ISTLT(8), ISTSQ(8), iSym, nB, nB2, nB3, nFroT,&
     &  nBasI, KEEPI, nBasJ, KEEPJ, iSymIJ, nBasIJ, nBasK, KEEPK, iSMax,&
     &  iSymL_, nBasL, KEEPL, nBasKL, IOPT, LPQ, IPQ, NPQ, IP, JQ, IRC, &
     &  ISX
!
      ISTSQ(1)=0
      ISTLT(1)=0
      Do iSym = 2, nSym
        nB  = nBas(iSym-1)
        nB2 = nB*nB
        nB3 = (nB2+nB)/2
        ISTSQ(iSym) = ISTSQ(iSym-1) + nB2
        ISTLT(iSym) = ISTLT(iSym-1) + nB3
      End Do
      nFroT = 0
      Do iSym = 1, nSym
        nFroT = nFroT + nFro(iSym)
      End Do

      nBasI  = nBas(iSymI)
      KEEPI  = KEEP(iSymI)
      ! nAuxI  = nAux(iSymI)
      nBasJ  = nBas(iSymJ)
      KEEPJ  = KEEP(iSymJ)
      ! nAuxJ  = nAux(iSymJ)
      iSymIJ = 1+iEor(iSymI-1,iSymJ-1)
      nBasIJ = nBasI*nBasJ
      If (iSymI == iSymJ) nBasIJ = (nBasI*(nBasI+1))/2
      If (nBasIJ == 0) Return

      nBasK  = nBas(iSymK)
      KEEPK  = KEEP(iSymK)
      ! nAuxK  = nAux(iSymK)
      iSMax  = iSymK
      If (iSymK == iSymI) iSMax = iSymJ
      iSymL_ = 1+iEor(iSymIJ-1,iSymK-1)
      IF (iSymL_ > iSMax) Return !! should not
      nBasL  = nBas(iSymL_)
      KEEPL  = KEEP(iSymL_)
      ! nAuxL  = nAux(iSymL_)
      nBasKL = nBasK*nBasL
      IF (iSymK == iSymL_) nBasKL = (nBasK*(nBasK+1))/2
      If (nBasKL == 0) Return

      ! INTEGRAL BLOCK EXCLUDED BY SETTING KEEP PARAMETERS?
      IF (KEEPI+KEEPJ+KEEPK+KEEPL /= 0) Return
      !! This will not work when the number of the inactive orbital is 0
!     IF (nAuxI+nAuxJ+nAuxK+nAuxL == 0) Return ! frozen orbitals

      !! (ij|kl)
      If (iSymI == iSymJ .AND. iSymI == iSymK) Then
        IOPT=1
        LPQ=0
        IPQ=0
        NPQ=0
        DO IP = 1, nBasI
          DO JQ = 1, IP
            IPQ=IPQ+1
            LPQ=LPQ+1
            IF (IPQ > NPQ) THEN
              CALL RDORD(IRC,IOPT,iSymI,iSymJ,iSymK,iSymL,X1,LBUF,NPQ)
              IF(IRC > 1) then
                WRITE(u6,*)' Error return code IRC=',IRC
                WRITE(u6,*)' from RDORD call, in FTWOI.'
                CALL Abend()
              end if
              IOPT=2
              IPQ=1
            ENDIF
            ISX=(IPQ-1)*nBasKL+1
!           ISF=ISTLT(iSymI)+LPQ
!           ISD=ISTLT(iSymI)+1
            !! coulomb
!           TEMP=DDOT_(KLB,X1(ISX),1,DLT(ISD),1)
!           FLT(ISF)=FLT(ISF)+TEMP
            !! exchange
            !! f_{ik} = (ij|kl)*D_{jl}
            !! read (ij|kl) for given ij as (kl)
            !! first  dgemv: (kl)*D_{l,i} -> f_{ki}
            !! second dgemv: (kl)*D_{l,i} -> f_{ki}

            !! For (mu nu), read (mu nu|rho sigma) as I(rho sigma)
            !! (mu j|rho sigma) = sum_{nu} C_{nu j} (mu nu|rho sigma)

            !! L_{mu i} = T_{ij}^{mu nu} * (mu rho | j sigma)
            !!          = T_{ij}^{rho sigma} * C_{nu j}
            !!                 * (mu rho | nu sigma)

            !! IP = mu, JQ = rho
            !! X2 = (nu sigma) --> X2' = C_{nu j}*(nu sigma) = (j sigma)
            !! T_{ij}^{rho sigma} * (j sigma)
            !!   -> U_{i rho} for (mu, rho) pairs

            CALL SQUARE (X1(ISX),X2(1),1,nBasK,nBasL)
!           If (DoCholesky) Then
!           Else
            !! (mu(ip) rho(jq) | nu sigma) -> (mu(ip) rho(jq) | j sigma)
            call dgemm_('T','N',nOcc,nBasT,nBasT,                       &
     &                  One,CMO,nBasT,X2,nBasT,                         &
     &                  Zero,WRK,nOcc)
            !! wrk(j,sigma) for the given mu(ip), mu(jq)
            !! T2AO(j,sigma,i,rho) = T_{ij}^{rho sigma}
            !! rather than L_{mu i}, L_{i mu} is computed
            !! L_{i mu} = wrk(j,sigma)*(T2AO(j,sigma,i,rho)
            call dgemv_('t',nOcc*nBasT,nOcc,                            &
     &                  One,t2ao(1,1,1,jq),nOcc*nBasT,wrk,1,            &
     &                  One,vlag(1,ip),1)
            if (ip /= jq) then
              call dgemv_('t',nOcc*nBasT,nOcc,                          &
     &                    One,t2ao(1,1,1,ip),nOcc*nBasT,wrk,1,          &
     &                    One,vlag(1,jq),1)
!           call dgemm_('T','T',nOcc,nBasT,nBasT,
!    *                  One,CMO,nBasT,X2,nBasT,
!    *                  One,WRK,nOcc)
!           call dgemv_('n',nOcc*nBasT,nOcc,
!    *                  One,t2ao(1,ip,1,1),nOcc*nBasT,wrk,1,
!    *                  One,vlag(1,jq),1)
            end if
!           End If

            !! DPT2AO -> FPT2AO transformation
            !! FPT2 = G(DPT2)
            Call FDGTRF(DPT2AO,FPT2AO)
            Call FDGTRF(DPT2CAO,FPT2CAO)
            If (nFroT /= 0) Then
              Call FDGTRF(DIA,FIFA)
              Call FDGTRF(DI,FIMO)
            End If
            !! Coulomb
!           Val = DDot_(nBasK*nBasL,X2,1,DPT2AO(ISTSQ(iSymI)+1),1)
!           iSF = ISTSQ(iSYmI) + iP + nBasI*(jQ-1)
!           FPT2AO(iSF) = FPT2AO(iSF) + Val
!           !! Exchange
!           iSF = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!           iSD = ISTSQ(iSymI) + (iP-1)*nBasI+1
!           CALL DGEMV_('N',nBasK,nBasL,
!    *                  -Half,X2,nBasK,DPT2AO(iSD),1,
!    *                  One,FPT2AO(iSF),1)
!           If (iP /= jQ) Then
!             iSF = ISTSQ(iSymI) + jQ + nBasJ*(iP-1)
!             FPT2AO(iSF) = FPT2AO(iSF) + Val
!             iSF = ISTSQ(iSymI) + (iP-1)*nBasI+1
!             iSD = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!             CALL DGEMV_('N',nBasK,nBasL,
!    *                    -Half,X2,nBasK,DPT2AO(iSD),1,
!    *                    One,FPT2AO(iSF),1)
!           End If
!
!           !! DPT2CAO -> FPT2CAO transformation
!           !! Coulomb
!           Val = DDot_(nBasK*nBasL,X2,1,DPT2CAO(ISTSQ(iSymI)+1),1)
!           iSF = ISTSQ(iSYmI) + iP + nBasI*(jQ-1)
!           FPT2CAO(iSF) = FPT2CAO(iSF) + Val
!           !! Exchange
!           iSF = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!           iSD = ISTSQ(iSymI) + (iP-1)*nBasI+1
!           CALL DGEMV_('N',nBasK,nBasL,
!    *                  -Half,X2,nBasK,DPT2CAO(iSD),1,
!    *                  One,FPT2CAO(iSF),1)
!           If (iP /= jQ) Then
!             iSF = ISTSQ(iSymI) + jQ + nBasJ*(iP-1)
!             FPT2CAO(iSF) = FPT2CAO(iSF) + Val
!             iSF = ISTSQ(iSymI) + (iP-1)*nBasI+1
!             iSD = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
!             CALL DGEMV_('N',nBasK,nBasL,
!    *                    -Half,X2,nBasK,DPT2CAO(iSD),1,
!    *                    One,FPT2CAO(iSF),1)
!           End If
          End Do
        End Do
      End If

!     ELSE IF ( IS == JS .AND. IS /= KS ) THEN
!       ! CASE 2: Integrals are of symmetry type (II/JJ)
!       ! Coulomb terms need to be accumulated only
!     ELSE IF ( IS == KS .AND. JS == LS ) THEN
!       ! CASE 3: Integrals are of symmetry type (IJ/IJ)
!       ! Exchange terms need to be accumulated only
!         IOPT=1
!         LPQ=0
!         IPQ=0
!         NPQ=0
!         DO IP=1,IB
!           DO JQ=1,JB
!             IPQ=IPQ+1
!             LPQ=LPQ+1
!             IF ( IPQ > NPQ ) THEN
!               CALL RDORD(IRC,IOPT,IS,JS,KS,LS,X1,LBUF,NPQ)
!               IF(IRC > 1) then
!                 WRITE(u6,*)' Error return code IRC=',IRC
!                 WRITE(u6,*)' from RDORD call, in FTWOI.'
!                 CALL Abend()
!               end if
!               IOPT=2
!               IPQ=1
!             ENDIF
!             ISX=(IPQ-1)*KLB+1
!             IF ( NFI /= 0 ) THEN
!               ISD=ISTSQ(IS)+(IP-1)*IB+1
!               ISF=ISTSQ(JS)+(JQ-1)*JB+1
!               CALL DGEMV_('N',LB,KB,(-Half*ExFac),X1(ISX),LB,
!    &                       DSQ(ISD),1,One,FSQ(ISF),1)
!             ENDIF
!             IF ( NFJ /= 0 ) THEN
!               ISD=ISTSQ(JS)+(JQ-1)*JB+1
!               ISF=ISTSQ(IS)+(IP-1)*IB+1
!               CALL DGEMV_('T',LB,KB,(-Half*ExFac),X1(ISX),LB,
!    &                       DSQ(ISD),1,One,FSQ(ISF),1)
!             ENDIF
!           End Do
!         End Do
!     End If

      RETURN

      Contains

      Subroutine FDGTRF(DD,FF)

      implicit none

      real(kind=wp), intent(in) :: DD(nBasT**2)
      real(kind=wp), intent(inout) :: FF(nBasT**2)

      real(kind=wp) :: Val
      real(kind=wp), external :: ddot_
      integer(kind=iwp) :: iSF, iSD

      !! Coulomb
      Val = DDot_(nBasK*nBasL,X2,1,DD(ISTSQ(iSymI)+1),1)
      iSF = ISTSQ(iSYmI) + iP + nBasI*(jQ-1)
      FF(iSF) = FF(iSF) + Val

      !! Exchange
      iSF = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
      iSD = ISTSQ(iSymI) + (iP-1)*nBasI+1
      CALL DGEMV_('N',nBasK,nBasL,                                      &
     &           -Half,X2,nBasK,DD(iSD),1,                              &
     &            One,FF(iSF),1)
      If (iP /= jQ) Then
        iSF = ISTSQ(iSymI) + jQ + nBasJ*(iP-1)
        FF(iSF) = FF(iSF) + Val
        iSF = ISTSQ(iSymI) + (iP-1)*nBasI+1
        iSD = ISTSQ(iSymI) + (jQ-1)*nBasJ+1
        CALL DGEMV_('N',nBasK,nBasL,                                    &
     &             -Half,X2,nBasK,DD(iSD),1,                            &
     &              One,FF(iSF),1)
      End If

      End Subroutine FDGTRF

      End Subroutine VVVOX
