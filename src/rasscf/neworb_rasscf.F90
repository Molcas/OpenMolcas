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
      SUBROUTINE NEWORB_RASSCF(CMOO,CMON,FP,FTR,VEC,WO,SQ,CMOX,D,OCCN)
!
!     RASSCF program, version IBM-3090: SX section
!
!     Purpose: To diagonalize the inactive,
!              and external parts of the Fock matrix FP=(FI+FA)
!              and order the eigenvalues and eigenvectors after energy.
!              The active orbitals are obtained by diagonalising the
!              corresponding RAS subspaces.
!              The new MO's are written onto JOBIPH
!              in address IADR15(2), and are used to produce the final
!              wave function. These are the orbitals printed in the
!              output section.
!     Called from SXCTL if IFINAL=1 (after last MC iteration)
!
!          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
!
!     Calling arguments:
!     D       : array of real, input
!               one body density matrix in MO-space
!     FP      : array of real, input
!               MCSCF-Fock matrix
!     FTR     : array of real
!               scratch array (local copy of FP, one sym. block only)
!     SQ      : array of real
!               scratch array (keep in sym block of Fock)
!     WO      : array of real
!               scratch array (keep in sym block of Fock)
!     VEC     : array of real
!               scratch array (eigenvectors of FTR)
!     CMOX    : array of real, input
!               scratch array (local copy of CMOO)
!     CMOO    : array of real, input
!               old MO-coefficients
!     CMON    : array of real, output
!               new MO-coefficients
!     OCNN    : array of real, input/output
!               MO occupation numbers

#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
#ifdef _HDF5_
      use mh5, only: mh5_put_dset
      use RASWfn, only: wfn_mocoef, wfn_occnum, wfn_orbene
#endif
      use gas_data, only: NGAS,NGSSH
      use rasscf_global, only: iFORDE, iOrbTyp, iOrdEM, iSupSM,         &
     &                         FDIAG, ixSym, iTRI, iADR15
#ifdef _DMRG_
      use rasscf_global, only: DoDMRG
#else
      use rasscf_global, only: DoBLOCKDMRG
#endif
      use PrintLevel, only: DEBUG
      use output_ras, only: LF,IPRLOC
      use general_data, only: NSYM,NTOT,JOBIPH,NASH,NBAS,NDEL,NFRO,NISH,&
     &                        NSSH,NTOT2

      IMPLICIT None

      Character(LEN=16), Parameter :: ROUTINE='NEWORB  '

      Real*8 CMOO(*),CMON(*),FP(*),FTR(*),VEC(*),                       &
     &          WO(*),SQ(*),D(*),OCCN(*),CMOX(*)

      Real*8 AVij, AVMx, Fact, FMin, Swap, VIJ
      Integer iPrLev, i, iAd15, iB, iBas, iGas, ii, iOff, iOrd, iST,    &
     &        iSTD, iSTFCK, iSTI, iSTM, iSTMO, iSTMO1, iSTMOA, iSYM,    &
     &        ixSymT, j, jSel, k, Min, NA, NA1, NAB, NABT, NAO, NAO2,   &
     &        NAT, NB, NBF, NBT, ND, NDNB, NDO, NEO, NEO1, NEO2, NF,    &
     &        NFI_, NFNB, NFO, NI, NI1, NIJ, NIO, NIO1, NIO2, NJ, NO1,  &
     &        NOO, NOT, NT, NTTR, NTU, NTUD, NU, NUT
! Local print level (if any)
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
!
!     IORD = 1 Ordering of inactive and secondary orbitals
!
      IF(ISUPSM.eq.1.AND.IORDEM.eq.0) IFORDE=0
      IORD=IFORDE

      IB=0
      ISTMO1=1
      ISTFCK=0
      ISTD=0

! A long loop over symmetry:
      DO ISYM=1,NSYM
       NBF=NBAS(ISYM)
       NFO=NFRO(ISYM)
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
       NEO=NSSH(ISYM)
       NOO=NFO+NIO+NAO
       NOT=NIO+NAO+NEO
       ISTMO=ISTMO1+NFO*NBF
!*********************************************************************
!      Frozen orbitals (move MO's to CMON, and set zero to FDIAG)
!*********************************************************************
       IF(NFO.NE.0) THEN
        NFNB=NBF*NFO
        CALL DCOPY_(NFNB,CMOO(ISTMO1),1,CMON(ISTMO1),1)
        DO  NF=1,NFO
         FDIAG(IB+NF)=0.0D0
         OCCN(IB+NF)=2.0D0
        END DO
       ENDIF
!*********************************************************************
!      Inactive part of the Fock matrix
!*********************************************************************
       IF(NIO.NE.0) THEN
!       MOVE FP TO TRIANGULAR FORM
        NIJ=0
        DO NI=1,NIO
         DO NJ=1,NI
          NIJ=NIJ+1
          FTR(NIJ)=FP(NIJ+ISTFCK)
          IF(IXSYM(IB+NFO+NI).NE.IXSYM(IB+NFO+NJ)) FTR(NIJ)=0.0D0
         END DO
        END DO
!       DIAGONALIZE
        NIO2=NIO**2
        CALL FZERO(VEC,NIO2)
        II=1
        DO NI=1,NIO
         VEC(II)=1.0D0
         II=II+NIO+1
        END DO
        CALL JACOB(FTR,VEC,NIO,NIO)
!
!       Transform molecular orbitals
!
        CALL DGEMM_('N','N',                                            &
     &              NBF,NIO,NIO,                                        &
     &              1.0d0,CMOO(ISTMO),NBF,                              &
     &              VEC,NIO,                                            &
     &              0.0d0,CMON(ISTMO),NBF)
!
!       Sort eigenvalues and orbitals after energy
!
!       Move eigenvalues to FDIAG and set occupation numbers to 2.
!
        II=0
        NO1=IB+NFO
        DO NI=1,NIO
         II=II+NI
         FDIAG(NO1+NI)=FTR(II)
         OCCN(NO1+NI)=2.0D0
        END DO
        IF(NIO.GT.1.AND.IORD.NE.0) THEN
         NIO1=NIO-1
         DO NI=1,NIO1
          NI1=NI+1
          MIN=NI
          DO NJ=NI1,NIO
           IF(FDIAG(NO1+NJ).LT.FDIAG(NO1+MIN)) MIN=NJ
          END DO
          IF(MIN.EQ.NI) GO TO 20
          FMIN=FDIAG(NO1+MIN)
          FDIAG(NO1+MIN)=FDIAG(NO1+NI)
          FDIAG(NO1+NI)=FMIN
! (SVC) added: extra nodig voor verandering ordening met supsym
          IXSYMT=IXSYM(NO1+MIN)
          IXSYM(NO1+MIN)=IXSYM(NO1+NI)
          IXSYM(NO1+NI)=IXSYMT

          ISTI=ISTMO+NBF*(NI-1)
          ISTM=ISTMO+NBF*(MIN-1)
          CALL DSWAP_(NBF,CMON(ISTI),1,CMON(ISTM),1)
20       CONTINUE
         END DO
        ENDIF
       ENDIF
!*********************************************************************
!      Active orbitals. Diagonalize corresponding diagonal block of D
!*********************************************************************
!.. dongxia make a loop for all active spaces
!
      ioff=0
      Do igas=1,ngas
       if(igas.gt.1)ioff=ioff+ngssh(igas-1,isym)
       IF(ngssh(igas,isym).NE.0) THEN
!       MOVE D TO TRIANGULAR FORM
        NTU=0
        ntud=istd+itri(ioff+1)
        DO NT=1,ngssh(igas,isym)
         ntud=ntud+ioff
! YM: change ntt --> nttr for the conflict if include rctfld.fh
         nttr=nt+nio+ioff
         DO NU=1,NT
          nut=nu+nio+ioff
          NTU=NTU+1
          ntud=ntud+1
          FTR(NTU)=D(ntud)
          IF(IXSYM(IB+NFO+nttr).NE.IXSYM(IB+NFO+nut)) ftr(ntu)=0.0d0
         END DO
        END DO

#ifdef _DMRG_
! YM:   If lRF, use Canonical instead of Natural, otherwise the energy expression will be
!       mismatch due to the FI matrix still in Canonical form.
!c YM:   The lRF flag can be removed if "TRACI" utility has MPS version
!c        If(.not.lRF.OR.iOrbTyp.EQ.2) Then
! YM:   Just use the same way as Block did, the NOs will be added later
        If(.not.doDMRG.OR.iOrbTyp.EQ.2) Then
#else
! NN.14 Skip making new orbitals for DMRG-CASSCF except for the case OutOrb = Canonical
!       because the DMRG is orbital variant.
        If(.NOT.DoBlockDMRG.OR.iOrbTyp.EQ.2) Then
#endif

!
!       DIAGONALIZE
        NAO2=ngssh(igas,isym)**2
        CALL FZERO(VEC,NAO2)
        II=1
        DO NT=1,ngssh(igas,isym)
         VEC(II)=1.0D0
         II=II+ngssh(igas,isym)+1
        END DO
        CALL JACOB(FTR,VEC,ngssh(igas,isym),ngssh(igas,isym))
!PAM01 Reorder to max-overlap agreement with input orbitals.
!PAM01 This prevents numerical difficulty with the subsequent transformation
!PAM01 of the CI vector, see TRAMAT routine.
!PAM01 Similar later, with RAS2, RAS3.
        DO I=1,ngssh(igas,isym)-1
         JSEL=I
         AVMX=ABS(VEC(I+ngssh(igas,isym)*(I-1)))
         DO J=I+1,ngssh(igas,isym)
          AVIJ=ABS(VEC(I+ngssh(igas,isym)*(J-1)))
          IF(AVIJ.GT.AVMX) THEN
            JSEL=J
            AVMX=AVIJ
          END IF
         END DO
!PAM01 Swap if necessary. Note phase control!!
         IF(JSEL.GT.I) THEN
           VIJ=VEC(I+ngssh(igas,isym)*(JSEL-1))
           IF(VIJ.GT.0.0D0) THEN
             DO K=1,ngssh(igas,isym)
               SWAP=VEC(K+ngssh(igas,isym)*(JSEL-1))
               VEC(K+ngssh(igas,isym)*(JSEL-1))=                        &
     &          -VEC(K+ngssh(igas,isym)*(I-1))
               VEC(K+ngssh(igas,isym)*(I-1))=SWAP
             END DO
           ELSE
             DO K=1,ngssh(igas,isym)
               SWAP=VEC(K+ngssh(igas,isym)*(JSEL-1))
               VEC(K+ngssh(igas,isym)*(JSEL-1))=                        &
     &          +VEC(K+ngssh(igas,isym)*(I-1))
               VEC(K+ngssh(igas,isym)*(I-1))=-SWAP
             END DO
           END IF
! Also swap eigenvalues:
           SWAP=FTR((JSEL*(JSEL+1))/2)
           FTR((JSEL*(JSEL+1))/2)=FTR((I*(I+1))/2)
           FTR((I*(I+1))/2)=SWAP
         ELSE
!PAM01 If swap is not needed, still apply phase control!!
           IF(VEC(I+ngssh(igas,isym)*(I-1)).LT.0.0D0) THEN
             DO K=1,ngssh(igas,isym)
              VEC(K+ngssh(igas,isym)*(I-1))=                            &
     &         -VEC(K+ngssh(igas,isym)*(I-1))
              VEC(K+ngssh(igas,isym)*I  )=                              &
     &         -VEC(K+ngssh(igas,isym)*I  )
             END DO
           END IF
         END IF
        END DO
!PAM01 We now have a MINIMAL pure rotation. This means that all
!PAM01 upper-left subdeterminants are >0 and as big as possible.
!PAM01 This eliminates numerical trouble in TRAMAT.

!
!       Transform molecular orbitals
!
        ISTMOA=ISTMO+NBF*(NIO+ioff)
        CALL DGEMM_('N','N',                                            &
     &              NBF,ngssh(igas,isym),ngssh(igas,isym),              &
     &              1.0d0,CMOO(ISTMOA),NBF,                             &
     &              VEC,ngssh(igas,isym),                               &
     &              0.0d0,CMON(ISTMOA),NBF)
! NN.14 Just copy CMO(Old) to CMO(new)
        Else ! If(.NOT.DoDMRG.OR.iOrbTyp.EQ.2)
          ISTMOA=ISTMO+NBF*(NIO+ioff)
         CALL DCOPY_(NBF*ngssh(igas,isym),CMOO(ISTMOA),1,CMON(ISTMOA),1)
        End If
!
!       Move eigenvalues to OCCN
!
        II=0
        NO1=IB+NFO+NIO+ioff
        DO NT=1,ngssh(igas,isym)
         II=II+NT
         OCCN(NO1+NT)=FTR(II)
         FDIAG(NO1+NT)=0.0D0
        END DO
!
!  FA:  no longer setting energies to 0 (though in principle ill-def).
!
!  IFG: Pick the diagonal elements of the transformed Fock matrix
        NFI_=(NFO+NIO+ioff)*(NFO+NIO+ioff+1)/2
        NO1=IB+NFO+NIO+ioff
        DO NT=1,ngssh(igas,isym)
         NFI_=NFI_+NFO+NIO+ioff
         DO NU=1,NT
          IF (NU.EQ.NT) THEN
           FACT=FP(NU+NFI_+ISTFCK)
          ELSE
           FACT=2.0d0*FP(NU+NFI_+ISTFCK)
          END IF
          DO II=1,ngssh(igas,isym)
           FDIAG(NO1+II)=FDIAG(NO1+II)+FACT                             &
     &                                *VEC(NT+(II-1)*ngssh(igas,isym))  &
     &                                *VEC(NU+(II-1)*ngssh(igas,isym))
          END DO
         END DO
         NFI_=NFI_+NT
        END DO

       ENDIF  ! end of if(NGAS(1).ne.0)
      END DO ! end loop over GAS spaces
!
!*********************************************************************
!      external part of the Fock matrix
!*********************************************************************
!
       IF(NEO.NE.0) THEN
!       Move fp to triangular form
        NAB=0
        DO NA=1,NEO
         DO NB=1,NA
          NAB=NAB+1
          NAT=NA+NIO+NAO
          NBT=NB+NIO+NAO
          NABT=ISTFCK+(NAT**2-NAT)/2+NBT
          FTR(NAB)=FP(NABT)
          IF(IXSYM(IB+NFO+NAT).NE.IXSYM(IB+NFO+NBT)) FTR(NAB)=0.0D0
         END DO
        END DO
!       DIAGONALIZE
        NEO2=NEO**2
        CALL FZERO(VEC,NEO2)
        II=1
        DO NA=1,NEO
         VEC(II)=1.0D0
         II=II+NEO+1
        END DO
        CALL JACOB(FTR,VEC,NEO,NEO)
!
!       Move eigenvalues to FDIAG and set occupation numbers to zero.
!
        II=0
        NO1=IB+NFO+NIO+NAO
        DO NA=1,NEO
         II=II+NA
         FDIAG(NO1+NA)=FTR(II)
         OCCN(NO1+NA)=0.0D0
        END DO
!
!       Transform molecular orbitals
!
        ISTMOA=ISTMO+NBF*(NIO+NAO)
        CALL DGEMM_('N','N',                                            &
     &              NBF,NEO,NEO,                                        &
     &              1.0d0,CMOO(ISTMOA),NBF,                             &
     &              VEC,NEO,                                            &
     &              0.0d0,CMON(ISTMOA),NBF)
!
!       Sort eigenvalues and orbitals after energy
!
        IF(NEO.GT.1.AND.IORD.NE.0) THEN
         NEO1=NEO-1
         DO NA=1,NEO1
          NA1=NA+1
          MIN=NA
          DO NB=NA1,NEO
           IF(FDIAG(NO1+NB).LT.FDIAG(NO1+MIN)) MIN=NB
          END DO
          IF(MIN.EQ.NA) GO TO 99
          FMIN=FDIAG(NO1+MIN)
          FDIAG(NO1+MIN)=FDIAG(NO1+NA)
          FDIAG(NO1+NA)=FMIN
! (SVC) added: extra nodig voor verandering ordening met supsym
          IXSYMT=IXSYM(NO1+MIN)
          IXSYM(NO1+MIN)=IXSYM(NO1+NA)
          IXSYM(NO1+NA)=IXSYMT

          ISTI=ISTMOA+NBF*(NA-1)
          ISTM=ISTMOA+NBF*(MIN-1)
          CALL DSWAP_(NBF,CMON(ISTI),1,CMON(ISTM),1)
99       CONTINUE
         END DO
        ENDIF
       ENDIF
!
!*********************************************************************
!      Deleted orbitals (move MO's and set zero to FDIAG and OCCN)
!*********************************************************************
       NDO=NDEL(ISYM)
       IF(NDO.NE.0) THEN
        NDNB=NDO*NBF
        IST=ISTMO1+NBF*(NOO+NEO)
        CALL DCOPY_(NDNB,CMOO(IST),1,CMON(IST),1)
        DO ND=1,NDO
         FDIAG(IB+NBF-NDO+ND)=0.0D0
         OCCN(IB+NBF-NDO+ND)=0.0D0
        END DO
       ENDIF
!
       IB=IB+NBF
       ISTFCK=ISTFCK+(NOT**2+NOT)/2
       ISTMO1=ISTMO1+NBF**2
       ISTD=ISTD+(NAO**2+NAO)/2
! End of a long loop over symmetry.
      END DO
!
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)' Diagonal elements of the FOCK matrix in NEWORB:'
        Write(LF,'(1X,10F11.6)') (FDIAG(I),I=1,NTOT)
      END IF
!
!*********************************************************************
!     Orthogonalise new orbitals
!*********************************************************************
      CALL ORTHO_rASSCF(WO,CMOX,CMON,SQ)
        If ( IPRLEV.ge.DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' CMO in NEWORB_RASSCF after diag and orthog'
         Write(LF,*) ' ---------------------'
         Write(LF,*)
         ioff=0
         Do iSym = 1,nSym
          iBas = nBas(iSym)
          if(iBas.ne.0) then
            write(6,*) 'Sym =', iSym
            do i= 1,iBas
              write(6,*) (CMON(ioff+iBas*(i-1)+j),j=1,iBas)
            end do
            iOff = iOff + (iBas*iBas)
          end if
         End Do
        End If
!
!     Save new CMOs on the JOBIPH
!
      IAD15=IADR15(2)
      CALL dDAFILE(JOBIPH,1,CMON,NTOT2,IAD15)
      CALL dDAFILE(JOBIPH,1,OCCN,NTOT,IAD15)

#ifdef _HDF5_
      call mh5_put_dset(wfn_mocoef,CMON)
      call mh5_put_dset(wfn_occnum,OCCN)
      call mh5_put_dset(wfn_orbene,FDIAG)
#endif

      RETURN
      END
