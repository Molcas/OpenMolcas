!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine NEWORB_RASSCF(CMOO,CMON,FP,FTR,VEC,WO,SQ,CMOX,D,OCCN)
! RASSCF program, version IBM-3090: SX section
!
! Purpose: To diagonalize the inactive,
!          and external parts of the Fock matrix FP=(FI+FA)
!          and order the eigenvalues and eigenvectors after energy.
!          The active orbitals are obtained by diagonalising the
!          corresponding RAS subspaces.
!          The new MO's are written onto JOBIPH
!          in address IADR15(2), and are used to produce the final
!          wave function. These are the orbitals printed in the
!          output section.
! Called from SXCTL if IFINAL=1 (after last MC iteration)
!
!      ********** IBM-3090 MOLCAS Release: 90 02 22 **********
!
! Calling arguments:
! D       : array of real, input
!           one body density matrix in MO-space
! FP      : array of real, input
!           MCSCF-Fock matrix
! FTR     : array of real
!           scratch array (local copy of FP, one sym. block only)
! SQ      : array of real
!           scratch array (keep in sym block of Fock)
! WO      : array of real
!           scratch array (keep in sym block of Fock)
! VEC     : array of real
!           scratch array (eigenvectors of FTR)
! CMOX    : array of real, input
!           scratch array (local copy of CMOO)
! CMOO    : array of real, input
!           old MO-coefficients
! CMON    : array of real, output
!           new MO-coefficients
! OCNN    : array of real, input/output
!           MO occupation numbers

use gas_data, only: NGAS, NGSSH
use PrintLevel, only: DEBUG
use output_ras, only: LF, IPRLOC
use general_data, only: NSYM, NTOT, JOBIPH, NASH, NBAS, NDEL, NFRO, NISH, NSSH, NTOT2
use rasscf_global, only: iFORDE, iOrbTyp, iOrdEM, iSupSM, FDIAG, ixSym, iTRI, iADR15
#ifdef _DMRG_
use qcmaquis_interface_cfg
use rasscf_global, only: DoDMRG
#else
use rasscf_global, only: DoBLOCKDMRG
#endif
#ifdef _HDF5_
use mh5, only: mh5_put_dset
use RASWfn, only: wfn_mocoef, wfn_occnum, wfn_orbene
#endif

implicit none
character(len=16), parameter :: ROUTINE = 'NEWORB  '
real*8 CMOO(*), CMON(*), FP(*), FTR(*), VEC(*), WO(*), SQ(*), D(*), OCCN(*), CMOX(*)
real*8 AVij, AVMx, Fact, FMin, Swap, VIJ
integer iPrLev, i, iAd15, iB, iBas, iGas, ii, iOff, iOrd, iST, iSTD, iSTFCK, iSTI, iSTM, iSTMO, iSTMO1, iSTMOA, iSYM, ixSymT, j, &
        jSel, k, Min, NA, NA1, NAB, NABT, NAO, NAO2, NAT, NB, NBF, NBT, ND, NDNB, NDO, NEO, NEO1, NEO2, NF, NFI_, NFNB, NFO, NI, &
        NI1, NIJ, NIO, NIO1, NIO2, NJ, NO1, NOO, NOT, NT, NTTR, NTU, NTUD, NU, NUT

! Local print level (if any)
IPRLEV = IPRLOC(4)
if (IPRLEV >= DEBUG) write(LF,*) ' Entering ',ROUTINE

! IORD = 1 Ordering of inactive and secondary orbitals

if ((ISUPSM == 1) .and. (IORDEM == 0)) IFORDE = 0
IORD = IFORDE

IB = 0
ISTMO1 = 1
ISTFCK = 0
ISTD = 0

! A long loop over symmetry:
do ISYM=1,NSYM
  NBF = NBAS(ISYM)
  NFO = NFRO(ISYM)
  NIO = NISH(ISYM)
  NAO = NASH(ISYM)
  NEO = NSSH(ISYM)
  NOO = NFO+NIO+NAO
  NOT = NIO+NAO+NEO
  ISTMO = ISTMO1+NFO*NBF
  !*********************************************************************
  !      Frozen orbitals (move MO's to CMON, and set zero to FDIAG)
  !*********************************************************************
  if (NFO /= 0) then
    NFNB = NBF*NFO
    call DCOPY_(NFNB,CMOO(ISTMO1),1,CMON(ISTMO1),1)
    do NF=1,NFO
      FDIAG(IB+NF) = 0.0d0
      OCCN(IB+NF) = 2.0d0
    end do
  end if
  !*********************************************************************
  !      Inactive part of the Fock matrix
  !*********************************************************************
  if (NIO /= 0) then
    ! MOVE FP TO TRIANGULAR FORM
    NIJ = 0
    do NI=1,NIO
      do NJ=1,NI
        NIJ = NIJ+1
        FTR(NIJ) = FP(NIJ+ISTFCK)
        if (IXSYM(IB+NFO+NI) /= IXSYM(IB+NFO+NJ)) FTR(NIJ) = 0.0d0
      end do
    end do
    ! DIAGONALIZE
    NIO2 = NIO**2
    call FZERO(VEC,NIO2)
    II = 1
    do NI=1,NIO
      VEC(II) = 1.0d0
      II = II+NIO+1
    end do
    call JACOB(FTR,VEC,NIO,NIO)

    ! Transform molecular orbitals

    call DGEMM_('N','N',NBF,NIO,NIO,1.0d0,CMOO(ISTMO),NBF,VEC,NIO,0.0d0,CMON(ISTMO),NBF)

    ! Sort eigenvalues and orbitals after energy

    ! Move eigenvalues to FDIAG and set occupation numbers to 2.

    II = 0
    NO1 = IB+NFO
    do NI=1,NIO
      II = II+NI
      FDIAG(NO1+NI) = FTR(II)
      OCCN(NO1+NI) = 2.0d0
    end do
    if ((NIO > 1) .and. (IORD /= 0)) then
      NIO1 = NIO-1
      do NI=1,NIO1
        NI1 = NI+1
        MIN = NI
        do NJ=NI1,NIO
          if (FDIAG(NO1+NJ) < FDIAG(NO1+MIN)) MIN = NJ
        end do
        if (MIN == NI) GO TO 20
        FMIN = FDIAG(NO1+MIN)
        FDIAG(NO1+MIN) = FDIAG(NO1+NI)
        FDIAG(NO1+NI) = FMIN
        ! (SVC) added: extra nodig voor verandering ordening met supsym
        IXSYMT = IXSYM(NO1+MIN)
        IXSYM(NO1+MIN) = IXSYM(NO1+NI)
        IXSYM(NO1+NI) = IXSYMT

        ISTI = ISTMO+NBF*(NI-1)
        ISTM = ISTMO+NBF*(MIN-1)
        call DSWAP_(NBF,CMON(ISTI),1,CMON(ISTM),1)
20      continue
      end do
    end if
  end if
  !*********************************************************************
  !      Active orbitals. Diagonalize corresponding diagonal block of D
  !*********************************************************************
  ! dongxia make a loop for all active spaces

  ioff = 0
  do igas=1,ngas
    if (igas > 1) ioff = ioff+ngssh(igas-1,isym)
    if (ngssh(igas,isym) /= 0) then
      ! MOVE D TO TRIANGULAR FORM
      NTU = 0
      ntud = istd+itri(ioff+1)
      do NT=1,ngssh(igas,isym)
        ntud = ntud+ioff
        nttr = nt+nio+ioff
        do NU=1,NT
          nut = nu+nio+ioff
          NTU = NTU+1
          ntud = ntud+1
          FTR(NTU) = D(ntud)
          if (IXSYM(IB+NFO+nttr) /= IXSYM(IB+NFO+nut)) ftr(ntu) = 0.0d0
        end do
      end do

#     ifdef _DMRG_
      ! YM:   If lRF, use Canonical instead of Natural, otherwise the energy expression will be
      !       mismatch due to the FI matrix still in Canonical form.
      ! YM:   The lRF flag can be removed if "TRACI" utility has MPS version
      !if ((.not. lRF) .or. (iOrbTyp == 2)) then
      ! YM:   Just use the same way as Block did, the NOs will be added later
      if ((.not. doDMRG) .or. (iOrbTyp == 2)) then
#     else
      ! NN.14 Skip making new orbitals for DMRG-CASSCF except for the case OutOrb = Canonical
      !       because the DMRG is orbital variant.
      if ((.not. DoBlockDMRG) .or. (iOrbTyp == 2)) then
#     endif

        ! DIAGONALIZE
        NAO2 = ngssh(igas,isym)**2
        call FZERO(VEC,NAO2)
        II = 1
        do NT=1,ngssh(igas,isym)
          VEC(II) = 1.0d0
          II = II+ngssh(igas,isym)+1
        end do
        call JACOB(FTR,VEC,ngssh(igas,isym),ngssh(igas,isym))
        !PAM01 Reorder to max-overlap agreement with input orbitals.
        !PAM01 This prevents numerical difficulty with the subsequent transformation
        !PAM01 of the CI vector, see TRAMAT routine.
        !PAM01 Similar later, with RAS2, RAS3.
        do I=1,ngssh(igas,isym)-1
          JSEL = I
          AVMX = abs(VEC(I+ngssh(igas,isym)*(I-1)))
          do J=I+1,ngssh(igas,isym)
            AVIJ = abs(VEC(I+ngssh(igas,isym)*(J-1)))
            if (AVIJ > AVMX) then
              JSEL = J
              AVMX = AVIJ
            end if
          end do
          !PAM01 Swap if necessary. Note phase control!!
          if (JSEL > I) then
            VIJ = VEC(I+ngssh(igas,isym)*(JSEL-1))
            if (VIJ > 0.0d0) then
              do K=1,ngssh(igas,isym)
                SWAP = VEC(K+ngssh(igas,isym)*(JSEL-1))
                VEC(K+ngssh(igas,isym)*(JSEL-1)) = -VEC(K+ngssh(igas,isym)*(I-1))
                VEC(K+ngssh(igas,isym)*(I-1)) = SWAP
              end do
            else
              do K=1,ngssh(igas,isym)
                SWAP = VEC(K+ngssh(igas,isym)*(JSEL-1))
                VEC(K+ngssh(igas,isym)*(JSEL-1)) = +VEC(K+ngssh(igas,isym)*(I-1))
                VEC(K+ngssh(igas,isym)*(I-1)) = -SWAP
              end do
            end if
            ! Also swap eigenvalues:
            SWAP = FTR((JSEL*(JSEL+1))/2)
            FTR((JSEL*(JSEL+1))/2) = FTR((I*(I+1))/2)
            FTR((I*(I+1))/2) = SWAP
          else
            !PAM01 If swap is not needed, still apply phase control!!
            if (VEC(I+ngssh(igas,isym)*(I-1)) < 0.0d0) then
              do K=1,ngssh(igas,isym)
                VEC(K+ngssh(igas,isym)*(I-1)) = -VEC(K+ngssh(igas,isym)*(I-1))
                VEC(K+ngssh(igas,isym)*I) = -VEC(K+ngssh(igas,isym)*I)
              end do
            end if
          end if
        end do
        !PAM01 We now have a MINIMAL pure rotation. This means that all
        !PAM01 upper-left subdeterminants are >0 and as big as possible.
        !PAM01 This eliminates numerical trouble in TRAMAT.

        ! Transform molecular orbitals
        ISTMOA = ISTMO+NBF*(NIO+ioff)
        call DGEMM_('N','N',NBF,ngssh(igas,isym),ngssh(igas,isym),1.0d0,CMOO(ISTMOA),NBF,VEC,ngssh(igas,isym),0.0d0,CMON(ISTMOA), &
                    NBF)
      else ! if ((.not. DoDMRG) .or. (iOrbTyp == 2))
        ! NN.14 Just copy CMO(Old) to CMO(new)
        ISTMOA = ISTMO+NBF*(NIO+ioff)
        call DCOPY_(NBF*ngssh(igas,isym),CMOO(ISTMOA),1,CMON(ISTMOA),1)
      end if

      ! Move eigenvalues to OCCN

      II = 0
      NO1 = IB+NFO+NIO+ioff
      do NT=1,ngssh(igas,isym)
        II = II+NT
        OCCN(NO1+NT) = FTR(II)
        FDIAG(NO1+NT) = 0.0d0
      end do

      ! FA:  no longer setting energies to 0 (though in principle ill-def).

      ! IFG: Pick the diagonal elements of the transformed Fock matrix
      NFI_ = (NFO+NIO+ioff)*(NFO+NIO+ioff+1)/2
      NO1 = IB+NFO+NIO+ioff
      do NT=1,ngssh(igas,isym)
        NFI_ = NFI_+NFO+NIO+ioff
        do NU=1,NT
          if (NU == NT) then
            FACT = FP(NU+NFI_+ISTFCK)
          else
            FACT = 2.0d0*FP(NU+NFI_+ISTFCK)
          end if
          do II=1,ngssh(igas,isym)
            FDIAG(NO1+II) = FDIAG(NO1+II)+FACT*VEC(NT+(II-1)*ngssh(igas,isym))*VEC(NU+(II-1)*ngssh(igas,isym))
          end do
        end do
        NFI_ = NFI_+NT
      end do

    end if  ! end of if (NGAS(1) /= 0)
  end do ! end loop over GAS spaces

  !*********************************************************************
  !      external part of the Fock matrix
  !*********************************************************************

  if (NEO /= 0) then
    ! Move fp to triangular form
    NAB = 0
    do NA=1,NEO
      do NB=1,NA
        NAB = NAB+1
        NAT = NA+NIO+NAO
        NBT = NB+NIO+NAO
        NABT = ISTFCK+(NAT**2-NAT)/2+NBT
        FTR(NAB) = FP(NABT)
        if (IXSYM(IB+NFO+NAT) /= IXSYM(IB+NFO+NBT)) FTR(NAB) = 0.0d0
      end do
    end do
    ! DIAGONALIZE
    NEO2 = NEO**2
    call FZERO(VEC,NEO2)
    II = 1
    do NA=1,NEO
      VEC(II) = 1.0d0
      II = II+NEO+1
    end do
    call JACOB(FTR,VEC,NEO,NEO)

    ! Move eigenvalues to FDIAG and set occupation numbers to zero.

    II = 0
    NO1 = IB+NFO+NIO+NAO
    do NA=1,NEO
      II = II+NA
      FDIAG(NO1+NA) = FTR(II)
      OCCN(NO1+NA) = 0.0d0
    end do

    ! Transform molecular orbitals

    ISTMOA = ISTMO+NBF*(NIO+NAO)
    call DGEMM_('N','N',NBF,NEO,NEO,1.0d0,CMOO(ISTMOA),NBF,VEC,NEO,0.0d0,CMON(ISTMOA),NBF)

    ! Sort eigenvalues and orbitals after energy

    if ((NEO > 1) .and. (IORD /= 0)) then
      NEO1 = NEO-1
      do NA=1,NEO1
        NA1 = NA+1
        MIN = NA
        do NB=NA1,NEO
          if (FDIAG(NO1+NB) < FDIAG(NO1+MIN)) MIN = NB
        end do
        if (MIN == NA) GO TO 99
        FMIN = FDIAG(NO1+MIN)
        FDIAG(NO1+MIN) = FDIAG(NO1+NA)
        FDIAG(NO1+NA) = FMIN
        ! (SVC) added: extra nodig voor verandering ordening met supsym
        IXSYMT = IXSYM(NO1+MIN)
        IXSYM(NO1+MIN) = IXSYM(NO1+NA)
        IXSYM(NO1+NA) = IXSYMT

        ISTI = ISTMOA+NBF*(NA-1)
        ISTM = ISTMOA+NBF*(MIN-1)
        call DSWAP_(NBF,CMON(ISTI),1,CMON(ISTM),1)
99      continue
      end do
    end if
  end if

  !*********************************************************************
  !      Deleted orbitals (move MO's and set zero to FDIAG and OCCN)
  !*********************************************************************
  NDO = NDEL(ISYM)
  if (NDO /= 0) then
    NDNB = NDO*NBF
    IST = ISTMO1+NBF*(NOO+NEO)
    call DCOPY_(NDNB,CMOO(IST),1,CMON(IST),1)
    do ND=1,NDO
      FDIAG(IB+NBF-NDO+ND) = 0.0d0
      OCCN(IB+NBF-NDO+ND) = 0.0d0
    end do
  end if

  IB = IB+NBF
  ISTFCK = ISTFCK+(NOT**2+NOT)/2
  ISTMO1 = ISTMO1+NBF**2
  ISTD = ISTD+(NAO**2+NAO)/2
  ! End of a long loop over symmetry.
end do

if (IPRLEV >= DEBUG) then
  write(LF,*) ' Diagonal elements of the FOCK matrix in NEWORB:'
  write(LF,'(1X,10F11.6)') (FDIAG(I),I=1,NTOT)
end if

!***********************************************************************
!     Orthogonalise new orbitals
!***********************************************************************
call ORTHO_rASSCF(WO,CMOX,CMON,SQ)
if (IPRLEV >= DEBUG) then
  write(LF,*)
  write(LF,*) ' CMO in NEWORB_RASSCF after diag and orthog'
  write(LF,*) ' ---------------------'
  write(LF,*)
  ioff = 0
  do iSym=1,nSym
    iBas = nBas(iSym)
    if (iBas /= 0) then
      write(6,*) 'Sym =',iSym
      do i=1,iBas
        write(6,*) (CMON(ioff+iBas*(i-1)+j),j=1,iBas)
      end do
      iOff = iOff+(iBas*iBas)
    end if
  end do
end if

! Save new CMOs on the JOBIPH

IAD15 = IADR15(2)
call dDAFILE(JOBIPH,1,CMON,NTOT2,IAD15)
call dDAFILE(JOBIPH,1,OCCN,NTOT,IAD15)

#ifdef _HDF5_
call mh5_put_dset(wfn_mocoef,CMON)
call mh5_put_dset(wfn_occnum,OCCN)
call mh5_put_dset(wfn_orbene,FDIAG)
#endif

return

end subroutine NEWORB_RASSCF
