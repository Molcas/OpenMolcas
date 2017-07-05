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
      SUBROUTINE NEWORB_RASSCF(CMOO,CMON,FP,FTR,VEC,WO,SQ,CMOX,D,OCCN)
C
C     RASSCF program, version IBM-3090: SX section
C
C     Purpose: To diagonalize the inactive,
C              and external parts of the Fock matrix FP=(FI+FA)
C              and order the eigenvalues and eigenvectors after energy.
C              The active orbitals are obtained by diagonalising the
C              corresponding RAS subspaces.
C              The new MO's are written onto JOBIPH
C              in address IADR15(2), and are used to produce the final
C              wave function. These are the orbitals printed in the
C              output section.
C     Called from SXCTL if IFINAL=1 (after last MC iteration)
C
C          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
C
C     Calling arguments:
*     D       : array of real, input
*               one body density matrix in MO-space
*     FP      : array of real, input
*               MCSCF-Fock matrix
*     FTR     : array of real
*               scratch array (local copy of FP, one sym. block only)
*     SQ      : array of real
*               scratch array (keep in sym block of Fock)
*     WO      : array of real
*               scratch array (keep in sym block of Fock)
*     VEC     : array of real
*               scratch array (eigenvectors of FTR)
*     CMOX    : array of real, input
*               scratch array (local copy of CMOO)
*     CMOO    : array of real, input
*               old MO-coefficients
*     CMON    : array of real, output
*               new MO-coefficients
*     OCNN    : array of real, input/output
*               MO occupation numbers

#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='NEWORB  ')
#include "SysDef.fh"
#include "raswfn.fh"
#include "rctfld.fh"

      DIMENSION CMOO(*),CMON(*),FP(*),FTR(*),VEC(*),
     *          WO(*),SQ(*),D(*),OCCN(*),CMOX(*)

      Call qEnter(routine)
C Local print level (if any)
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
C
C     IORD = 1 Ordering of inactive and secondary orbitals
C
      IF(ISUPSM.eq.1.AND.IORDEM.eq.0) IFORDE=0
      IORD=IFORDE

      IB=0
      ISTMO1=1
      ISTFCK=0
      ISTD=0

* A long loop over symmetry:
      DO ISYM=1,NSYM
       NBF=NBAS(ISYM)
       NFO=NFRO(ISYM)
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
       NEO=NSSH(ISYM)
       NOO=NFO+NIO+NAO
       NOT=NIO+NAO+NEO
       ISTMO=ISTMO1+NFO*NBF
C*********************************************************************
C      Frozen orbitals (move MO's to CMON, and set zero to FDIAG)
C*********************************************************************
       IF(NFO.NE.0) THEN
        NFNB=NBF*NFO
        CALL DCOPY_(NFNB,CMOO(ISTMO1),1,CMON(ISTMO1),1)
        DO  NF=1,NFO
         FDIAG(IB+NF)=0.0D0
         OCCN(IB+NF)=2.0D0
        END DO
       ENDIF
C*********************************************************************
C      Inactive part of the Fock matrix
C*********************************************************************
       IF(NIO.NE.0) THEN
C       MOVE FP TO TRIANGULAR FORM
        NIJ=0
        DO NI=1,NIO
         DO NJ=1,NI
          NIJ=NIJ+1
          FTR(NIJ)=FP(NIJ+ISTFCK)
          IF(IXSYM(IB+NFO+NI).NE.IXSYM(IB+NFO+NJ)) FTR(NIJ)=0.0D0
         END DO
        END DO
C       DIAGONALIZE
        NIO2=NIO**2
        CALL VCLR(VEC,1,NIO2)
        II=1
        DO NI=1,NIO
         VEC(II)=1.0D0
         II=II+NIO+1
        END DO
        CALL JACOB(FTR,VEC,NIO,NIO)
C
C       Transform molecular orbitals
C
        CALL DGEMM_('N','N',
     &              NBF,NIO,NIO,
     &              1.0d0,CMOO(ISTMO),NBF,
     &              VEC,NIO,
     &              0.0d0,CMON(ISTMO),NBF)
C
C       Sort eigenvalues and orbitals after energy
C
C       Move eigenvalues to FDIAG and set occupation numbers to 2.
C
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
* (SVC) added: extra nodig voor verandering ordening met supsym
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
C*********************************************************************
C      Active orbitals. Diagonalize corresponding diagonal block of D
C*********************************************************************
C.. dongxia make a loop for all active spaces
C
      ioff=0
      Do igas=1,ngas
       if(igas.gt.1)ioff=ioff+ngssh(igas-1,isym)
       IF(ngssh(igas,isym).NE.0) THEN
C       MOVE D TO TRIANGULAR FORM
        NTU=0
        ntud=istd+itri(ioff+1)
        DO NT=1,ngssh(igas,isym)
         ntud=ntud+ioff
C YM: change ntt --> nttr for the confict if include rctfld.fh
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
C YM:   If lRF, use Canonical instead of Natural, otherwise the energy expression will be
c       mismatch due to the FI matrix still in Canonical form.
cc YM:   The lRF flag can be removed if "TRACI" utility has MPS version
cc        If(.not.lRF.OR.iOrbTyp.EQ.2) Then
c YM:   Just use the same way as Block did, the NOs will be added later
        If(.not.doDMRG.OR.iOrbTyp.EQ.2) Then
#else
C NN.14 Skip making new orbitals for DMRG-CASSCF except for the case OutOrb = Canonical
C       because the DMRG is orbital variant.
        If(.NOT.DoBlockDMRG.OR.iOrbTyp.EQ.2) Then
#endif

C
C       DIAGONALIZE
        NAO2=ngssh(igas,isym)**2
        CALL VCLR(VEC,1,NAO2)
        II=1
        DO NT=1,ngssh(igas,isym)
         VEC(II)=1.0D0
         II=II+ngssh(igas,isym)+1
        END DO
        CALL JACOB(FTR,VEC,ngssh(igas,isym),ngssh(igas,isym))
CPAM01 Reorder to max-overlap agreement with input orbitals.
CPAM01 This prevents numerical difficulty with the subsequent transformation
CPAM01 of the CI vector, see TRAMAT routine.
CPAM01 Similar later, with RAS2, RAS3.
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
CPAM01 Swap if necessary. Note phase control!!
         IF(JSEL.GT.I) THEN
           VIJ=VEC(I+ngssh(igas,isym)*(JSEL-1))
           IF(VIJ.GT.0.0D0) THEN
             DO K=1,ngssh(igas,isym)
               SWAP=VEC(K+ngssh(igas,isym)*(JSEL-1))
               VEC(K+ngssh(igas,isym)*(JSEL-1))=
     &          -VEC(K+ngssh(igas,isym)*(I-1))
               VEC(K+ngssh(igas,isym)*(I-1))=SWAP
             END DO
           ELSE
             DO K=1,ngssh(igas,isym)
               SWAP=VEC(K+ngssh(igas,isym)*(JSEL-1))
               VEC(K+ngssh(igas,isym)*(JSEL-1))=
     &          +VEC(K+ngssh(igas,isym)*(I-1))
               VEC(K+ngssh(igas,isym)*(I-1))=-SWAP
             END DO
           END IF
C Also swap eigenvalues:
           SWAP=FTR((JSEL*(JSEL+1))/2)
           FTR((JSEL*(JSEL+1))/2)=FTR((I*(I+1))/2)
           FTR((I*(I+1))/2)=SWAP
         ELSE
CPAM01 If swap is not needed, still apply phase control!!
           IF(VEC(I+ngssh(igas,isym)*(I-1)).LT.0.0D0) THEN
             DO K=1,ngssh(igas,isym)
              VEC(K+ngssh(igas,isym)*(I-1))=
     &         -VEC(K+ngssh(igas,isym)*(I-1))
              VEC(K+ngssh(igas,isym)*I  )=
     &         -VEC(K+ngssh(igas,isym)*I  )
             END DO
           END IF
         END IF
        END DO
CPAM01 We now have a MINIMAL pure rotation. This means that all
CPAM01 upper-left subdeterminants are >0 and as big as possible.
CPAM01 This eliminates numerical trouble in TRAMAT.

C
C       Transform molecular orbitals
C
        ISTMOA=ISTMO+NBF*(NIO+ioff)
        CALL DGEMM_('N','N',
     &              NBF,ngssh(igas,isym),ngssh(igas,isym),
     &              1.0d0,CMOO(ISTMOA),NBF,
     &              VEC,ngssh(igas,isym),
     &              0.0d0,CMON(ISTMOA),NBF)
C NN.14 Just copy CMO(Old) to CMO(new)
        Else ! If(.NOT.DoDMRG.OR.iOrbTyp.EQ.2)
          ISTMOA=ISTMO+NBF*(NIO+ioff)
         CALL DCOPY_(NBF*ngssh(igas,isym),CMOO(ISTMOA),1,CMON(ISTMOA),1)
        End If
C
C       Move eigenvalues to OCCN and set energies to zero (not defined).
C
        II=0
        NO1=IB+NFO+NIO+ioff
        DO NT=1,ngssh(igas,isym)
         II=II+NT
         OCCN(NO1+NT)=FTR(II)
         FDIAG(NO1+NT)=0.0D0
        END DO

       ENDIF  ! end of if(NGAS(1).ne.0)
      END DO ! end loop over GAS spaces
C
C*********************************************************************
C      external part of the Fock matrix
C*********************************************************************
C
       IF(NEO.NE.0) THEN
C       Move fp to triangular form
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
C       DIAGONALIZE
        NEO2=NEO**2
        CALL VCLR(VEC,1,NEO2)
        II=1
        DO NA=1,NEO
         VEC(II)=1.0D0
         II=II+NEO+1
        END DO
        CALL JACOB(FTR,VEC,NEO,NEO)
C
C       Move eigenvalues to FDIAG and set occupation numbers to zero.
C
        II=0
        NO1=IB+NFO+NIO+NAO
        DO NA=1,NEO
         II=II+NA
         FDIAG(NO1+NA)=FTR(II)
         OCCN(NO1+NA)=0.0D0
        END DO
C
C       Transform molecular orbitals
C
        ISTMOA=ISTMO+NBF*(NIO+NAO)
        CALL DGEMM_('N','N',
     &              NBF,NEO,NEO,
     &              1.0d0,CMOO(ISTMOA),NBF,
     &              VEC,NEO,
     &              0.0d0,CMON(ISTMOA),NBF)
C
C       Sort eigenvalues and orbitals after energy
C
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
* (SVC) added: extra nodig voor verandering ordening met supsym
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
C
C*********************************************************************
C      Deleted orbitals (move MO's and set zero to FDIAG and OCCN)
C*********************************************************************
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
C
       IB=IB+NBF
       ISTFCK=ISTFCK+(NOT**2+NOT)/2
       ISTMO1=ISTMO1+NBF**2
       ISTD=ISTD+(NAO**2+NAO)/2
* End of a long loop over symmetry.
      END DO
C
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)' Diagonal elements of the FOCK matrix in NEWORB:'
        Write(LF,'(1X,10F11.6)') (FDIAG(I),I=1,NTOT)
      END IF
C
C*********************************************************************
C     Orthogonalise new orbitals
C*********************************************************************
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
C
C     Save new CMOs on the JOBIPH
C
      IAD15=IADR15(2)
      CALL dDAFILE(JOBIPH,1,CMON,NTOT2,IAD15)
      CALL dDAFILE(JOBIPH,1,OCCN,NTOT,IAD15)

#ifdef _HDF5_
      call mh5_put_dset(wfn_mocoef,CMON)
      call mh5_put_dset(wfn_occnum,OCCN)
      call mh5_put_dset(wfn_orbene,FDIAG)
#endif

      CALL QEXIT('NEWORB')
      RETURN
      END
