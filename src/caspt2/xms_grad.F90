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

subroutine XMS_Grad(H0,U0,UEFF,OMGDER)

use Index_Functions, only: nTri_Elem, nTri3_Elem
use sguga_states, only: SGS
use caspt2_global, only: CLag, CLagFull, CMOPT2, do_csf, do_nac, DPT2_tot, FIFA, FIFA_all, FIFASA_all, if_equalW, iRoot1, iRoot2, &
                         NDREF, nOLag, OLag, TORB, weight
use caspt2_module, only: ENERGY, IFDW, IFRMS, IFXMS, MXCI, NASHT, NBAS, NBAST, NBSQT, NCONF, NDEL, NFRO, NISH, NROOTS, NSTATE, &
                         ORBIN, STSYM, ZETA
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: H0(nState,nState), U0(nState,nState), OMGDER(nState,nState)
real(kind=wp), intent(inout) :: UEFF(nState,nState)
integer(kind=iwp) :: I, iStat, iState, jAsh, jStat, kStat, nBasI, nCor, nLev, nOrbI, NTG1, NTG3
real(kind=wp) :: EDIFF, EEI, EEJ, EINACT, fact, OVL, Scal, TRC, Wgt
real(kind=wp), allocatable :: CI1(:), CI2(:), DG1(:), DG2(:), DG3(:), DPT2(:), DPT2_AO(:), G1(:,:), RDMEIG(:,:), RDMSA(:,:), &
                              SGM1(:), SGM2(:), SLag(:,:), TG1(:,:), TG2(:), Trf(:), WRK1(:), WRK2(:)
real(kind=wp), external :: ddot_
integer(kind=iwp), parameter:: jstate=1

nLev = SGS(jstate)%nLev

call mma_allocate(SLag,nState,nState,Label='SLag')

! The XMS rotation applies to any variants: XMS-CASPT2, XDW-CASPT2,
! and RMS-CASPT2.

if (IFXMS .or. IFRMS) then
  call mma_allocate(CI1,nConf,Label='CI1')
  call mma_allocate(CI2,nConf,Label='CI2')
  call mma_allocate(SGM1,MXCI,Label='SGM1')
  call mma_allocate(SGM2,MXCI,Label='SGM2')
  call mma_allocate(TG1,nAshT,nAshT,Label='TG1')
  call mma_allocate(TG2,nAshT**4,Label='TG2')

  call mma_allocate(DG1,nAshT**2,Label='DG1')
  call mma_allocate(DG2,nAshT**4,Label='DG2')

  ! ----- Construct pseudo-density matrix -----

  !! First, we need to consider the derivative of the rotation
  !! of the CASSCF energy (Heff[1]).

  !! Forward transformation of UEFF (CASSCF basis to XMS basis)
  !! SLag is used as a working array
  call DGEMM_('T','N',nState,nState,nState,One,U0,nState,UEFF,nState,Zero,SLag,nState)
  UEFF(:,:) = SLag(:,:)

  !! Then actual calculation
  CLag(1:nconf,1:nstate) = Zero
  do iStat=1,nState
    call LoadCI_XMS('C',0,nConf,nState,CI1,iStat,U0)
    do jStat=1,nState
      call LoadCI_XMS('C',0,nConf,nState,CI2,jStat,U0)

      call Dens2T_RPT2(NLEV,NCONF,MXCI,CI1,CI2,SGM1,SGM2,TG1,TG2)
      TG1(:,:) = TG1(:,:)*Half
      TG2(:) = TG2(:)*Half
      DG1(:) = Zero
      DG2(:) = Zero
      call CnstInt(1,DG1,DG2)

      if (do_nac) then
        Scal = (UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)+UEFF(jStat,iRoot1)*UEFF(iStat,iRoot2))*Half
      else
        Scal = UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)
      end if
      if (IFDW .and. (zeta >= Zero)) scal = scal+OMGDER(iStat,jStat)
      DG1(:) = DG1(:)*Scal
      DG2(:) = DG2(:)*Scal

      call mma_allocate(DG3,nAshT**6,Label='DG3')
      DG3(:) = Zero

      NTG1 = NASHT**2
      NTG3 = nTri3_Elem(NTG1)
      OVL = Zero
      call DERTG3(.false.,STSYM,STSYM,NCONF,NASHT,CI1,CI2,OVL,DG1,DG2,NTG3,DG3,CLag(1,iStat),CLag(1,jStat))
      call mma_deallocate(DG3)
    end do
  end do
  call mma_deallocate(SGM2)
  call mma_deallocate(TG2)

  !! Back transformation of UEFF (XMS basis to CASSCF basis)
  !! SLag is used as a working array
  call DGEMM_('N','N',nState,nState,nState,One,U0,nState,UEFF,nState,Zero,SLag,nState)
  UEFF(:,:) = SLag(:,:)
  call mma_deallocate(DG1)
  call mma_deallocate(DG2)

  !! Add to the full CI Lagrangian
  !! CLag is in quasi-canonical basis, so transformation
  !! to natural basis is required.
  if (ORBIN == 'TRANSFOR') then
    do iState=1,nState
      call CLagX_TrfCI(nConf,CLag(1,iState))
    end do
  end if
  CLagFull(:,:) = CLagFull(:,:)+CLag(:,:)

  !! Compute the Lagrange multiplier for XMS
  !! The diagonal element is always zero.
  !! The code has an additional scaling with 0.5,
  !! because some contributions are doubled.
  SLag(:,:) = Zero
  do iStat=1,nState
    call LoadCI_XMS('N',0,nConf,nState,CI1,iStat,U0)
    EEI = H0(iStat,iStat)
    do jStat=1,iStat
      if (iStat == jStat) then
        SLag(iStat,jStat) = Zero
      else
        call LoadCI_XMS('N',0,nConf,nState,CI2,jStat,U0)
        EEJ = H0(jStat,jStat)

        Scal = DDOT_(nConf,CI1,1,CLagFull(1,jStat),1)-DDOT_(nConf,CI2,1,CLagFull(1,iStat),1)
        if (do_csf) then
          !! JCTC 2017, 13, 2561: eq.(66)
          !! iStat and jStat: XMS
          !! kStat and lStat: CASSCF
          fact = Zero
          do kStat=1,nState
            fact = fact+sum((UEFF(kStat,iRoot1)*UEFF(:,iRoot2)-UEFF(kStat,iRoot2)*UEFF(:,iRoot1))*Half*U0(kStat,iStat)*U0(:,jStat))
          end do
          Scal = Scal+fact*(ENERGY(iRoot2)-ENERGY(iRoot1))*Two
        end if
        Scal = Quart*Scal/(EEJ-EEI)
        SLag(iStat,jStat) = Scal
        SLag(jStat,iStat) = Scal
      end if
    end do
  end do

  !! Either subtract the CI derivative of Heff(1)
  !! or add off-diagonal couplings in rhs_sa (Z-vector)
  !CLagFull = CLagFull-CLag

  !! Finally, construct the pseudo-density matrix
  !! d = \sum_{ST} w_{ST} * d_{ST}
  call mma_allocate(G1,nAshT,nAshT,Label='G1')
  G1(:,:) = Zero
  do iStat=1,nState
    call LoadCI_XMS('C',0,nConf,nState,CI1,iStat,U0)
    do jStat=1,nState
      if (abs(SLag(iStat,jStat)) <= 1.0e-10_wp) cycle
      call LoadCI_XMS('C',0,nConf,nState,CI2,jStat,U0)

      call Dens1T_RPT2(CI1,CI2,SGM1,TG1,nLev)
      G1(:,:) = G1(:,:)+Two*SLag(iStat,jStat)*TG1(:,:)
    end do
  end do

  call mma_deallocate(CI1)
  call mma_deallocate(CI2)
  call mma_deallocate(SGM1)
  call mma_deallocate(TG1)

  ! ----- Calculate orbital derivatives -----

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
  Wgt = One/nState
  do iState=1,nState
    call LoadCI(CI1,iState)
    call POLY1(CI1,nConf)
    call GETDREF(WRK2,nDRef)
    WRK1(1:nDRef) = WRK1(1:nDRef)+Wgt*WRK2(1:nDRef)
  end do
  call mma_deallocate(CI1)
  call SQUARE(WRK1,RDMSA,1,nAshT,nAshT)

  nOrbI = nBas(1)-nDel(1) !! nOrb(1)
  nBasI = nBas(1)
  !call SQUARE(FIFA,FIFA_all,1,nOrbI,nOrbI)
  !! FIFASA_all is in natural orbital basis
  Trf(1:NBSQT) = Zero
  call CnstTrf(NBSQT,TOrb,Trf)

  !! FIFA: natural -> quasi-canonical
  if (IFDW .or. IFRMS) then
    call DGemm_('T','N',nOrbI,nOrbI,nOrbI,One,Trf,nBasI,FIFASA_all,nOrbI,Zero,WRK1,nOrbI)
    call DGemm_('N','N',nOrbI,nOrbI,nOrbI,One,WRK1,nOrbI,Trf,nBasI,Zero,FIFASA_all,nOrbI)
  end if
  FIFA_all(1:NBSQT) = FIFASA_all(1:NBSQT)

  !! Orbital derivatives of FIFA
  !! Both explicit and implicit orbital derivatives are computed
  !! Also, compute G(D) and put the active contribution to RDMEIG
  OLag(:) = Zero
  call EigDer2(NBSQT,nAshT,RDMEIG,Trf,FIFA_all,RDMSA,G1,WRK1,WRK2)

  !! Add to PT2 density
  !! No inactive contributions. Correct as long as CASSCF CI
  !! vector are orthogonal.
  call AddDEPSA(NBSQT,nAshT,DPT2,G1)
  call DPT2_TrfStore(One,NBSQT,DPT2,DPT2_tot,Trf,WRK1)

  !! Finalize OLag (anti-symetrize) and construct WLag
  call OLagFinal(nOLag,NBSQT,OLag,Trf)

  if (.not. if_equalW) then
    call mma_allocate(DPT2_AO,NBSQT,Label='DPT2_AO')

    !! AddDEPSA considers the frozen orbital, whereas DPT2_Trf
    !! does not. In any case, construct DPT2 again.
    DPT2(:) = Zero
    call DPT2_Trf(NBSQT,nAshT,DPT2,DPT2_AO,CMOPT2,G1,WRK1)
    !! Construct the SCF density
    WRK1(1:nDRef) = Zero
    call mma_allocate(CI1,nConf,Label='CI1')
    do iState=1,nState
      call LoadCI_XMS('N',1,nConf,nState,CI1,iState,U0)
      call POLY1(CI1,nConf)
      call GETDREF(WRK2,nDRef)
      WRK1(1:nDRef) = WRK1(1:nDRef)+Weight(iState)*WRK2(1:nDRef)
    end do
    call mma_deallocate(CI1)
    !! WRK2 is the SCF density (for nstate=nroots)
    call SQUARE(WRK1,WRK2,1,nAshT,nAshT)
    RDMSA(:,:) = RDMSA(:,:)-reshape(WRK2(1:nAshT**2),[nAshT,nAshT])
    !! Construct the SS minus SA density matrix in WRK1
    call OLagFroD(NBSQT,nAshT,WRK1,WRK2,RDMSA,Trf)
    !! Subtract the inactive part
    WRK1(1:nBasT**2) = WRK1(1:nBasT**2)-WRK2(1:nBasT**2)
    !! Save
    call CnstAB_SSDM(NBSQT,DPT2_AO,WRK1)
    call mma_deallocate(DPT2_AO)
  end if

  call mma_deallocate(RDMSA)
  call mma_deallocate(WRK1)
  call mma_deallocate(WRK2)

  call mma_deallocate(DPT2)

  ! ----- Calculate CI derivatives -----

  !! use quasi-canonical CSF rather than natural CSF
  CLag(:,:) = Zero

  !! 1) Explicit CI derivative
  !! a: Extract FIFA in the AS for explicit CI derivative
  !!    Note that this FIFA uses state-averaged density matrix
  !call sqprt(fifa_all,nbast)
  do jAsh=1,nAshT
    G1(:,jAsh) = FIFA_all(nFro(1)+nIsh(1)+nBas(1)*(nFro(1)+nIsh(1)+jAsh-1)+1:nFro(1)+nIsh(1)+nBas(1)*(nFro(1)+nIsh(1)+jAsh-1)+nAshT)
  end do
  !call sqprt(g1),nasht)
  !! Transform quasi-canonical to natural
  nCor = nFro(1)+nIsh(1)
  !call sqprt(trf,nbast)
  call DGemm_('N','N',nAshT,nAshT,nAshT,One,Trf(1+nBasT*nCor+nCor),nBasT,G1,nAshT,Zero,FIFA_all,nAshT)
  !call sqprt(FIFA_all,nasht)
  call DGemm_('N','T',nAshT,nAshT,nAshT,One,FIFA_all,nAshT,Trf(1+nBasT*nCor+nCor),nBasT,Zero,G1,nAshT)
  !call sqprt(g1,nasht)
  call DGemm_('N','N',nAshT,nAshT,nAshT,One,Trf(1+nBasT*nCor+nCor),nBasT,RDMEIG,nAshT,Zero,FIFA_all,nAshT)
  call DGemm_('N','T',nAshT,nAshT,nAshT,One,FIFA_all,nAshT,Trf(1+nBasT*nCor+nCor),nBasT,Zero,RDMEIG,nAshT)
  call mma_deallocate(Trf)

  !! b: FIFA in the inactive space
  TRC = Zero
  !do ISYM=1,NSYM
  do I=1,NISH(1) ! ISYM)
    !II = IOFF(ISYM)+nTri_Elem(I)
    TRC = TRC+FIFA(nTri_Elem(I))
  end do
  !end do
  ! Contribution from inactive orbitals:
  EINACT = Two*TRC
  !! This EINACT may be wrong. Perhaps, FIFA has to be
  !! back-transformed to natural orbital basis. However, this does
  !! not contribute to the final gradient as long as all the
  !! (internal) CI vectors are orthogonal.

  !! c: Finally, compute explicit CI derivative
  !! y_{I,T} = w_{ST} <I|f|S>
  !! Here, G1 is FIFA = ftu
  call CLagEigT(CLag,G1,SLag,EINACT)

  !! 2) Implicit CI derivative
  call CLagEig(.false.,.true.,nConf,nRoots,nState,nLev,CLag,RDMEIG)

# ifdef _MOLCAS_MPP_
  if (is_real_par()) call GADGOP(CLag,nConf*nState,'+')
# endif

  call mma_deallocate(RDMEIG)
  call mma_deallocate(G1)

end if

if (do_csf) then
  !! Eq (68)
  !! I think this contributes even in XMS-CASPT2.
  !! However, maybe I implemented other terms in a different way
  !! from BAGEL.
  call mma_allocate(CI1,nConf,Label='CI1')
  !! Do in the natural orbital in the XMS basis
  !! because CLag is like that
  !! CASSCF -> XMS
  call DGEMM_('T','N',nState,nState,nState,One,U0,nState,UEFF,nState,Zero,SLag,nState)
  UEFF(:,:) = SLag(:,:)

  EDIFF = ENERGY(iRoot2)-ENERGY(iRoot1)
  do iStat=1,nState
    call LoadCI_XMS('N',0,nConf,nState,CI1,iStat,U0)
    do jStat=1,nState
      CLag(:,jStat) = CLag(:,jStat)+(UEFF(iStat,iRoot1)*UEFF(jStat,iRoot2)-UEFF(iStat,iRoot2)*UEFF(jStat,iRoot1))*Half*EDIFF*CI1(:)
    end do
  end do

  !! XMS -> CASSCF
  call DGEMM_('N','N',nState,nState,nState,One,U0,nState,UEFF,nState,Zero,SLag,nState)
  UEFF(:,:) = SLag(:,:)

  call mma_deallocate(CI1)
end if

call mma_deallocate(SLag)

!! Finally, add to the total CI derivative array
CLagFull(:,:) = CLagFull(:,:)+CLag(:,:)

end subroutine XMS_Grad
