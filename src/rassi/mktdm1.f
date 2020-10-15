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
      SUBROUTINE MKTDM1(LSYM1,MPLET1,MSPROJ1,IFSBTAB1,
     &    LSYM2,MPLET2,MSPROJ2,IFSBTAB2,ISSTAB,MAPORB,
     &    DET1,DET2,SIJ,NASHT,TDM1,TSDM1,WTDM1,ISTATE,JSTATE,
     &    job1,job2,ist,jst)

      !> module dependencies
#ifdef _DMRG_
      use rasscf_data, only: doDMRG
      use qcmaquis_interface_cfg
      use qcmaquis_interface_utility_routines, only:
     &    pretty_print_util
      use qcmaquis_info
      use qcmaquis_interface_mpssi
#endif

      IMPLICIT NONE
      INTEGER LSYM1,MPLET1,MSPROJ1,LSYM2,MPLET2,MSPROJ2
#ifdef _DMRG_
      Real*8, Allocatable:: TDMAA(:), TDMBB(:)
#endif
      INTEGER IFSBTAB1(*),IFSBTAB2(*),ISSTAB(*),MAPORB(*)
      INTEGER IORB,ISORB,ISYOP,ITABS,IUABS,JORB,JSORB,LORBTB
      INTEGER MS2OP,NASHT,NASORB,NSPD1
      INTEGER, INTENT(IN) :: ISTATE, JSTATE,job1,job2,ist,jst
      REAL*8 DET1(*),DET2(*)
      REAL*8 SIJ,TDM1(NASHT,NASHT),TSDM1(NASHT,NASHT),WTDM1(NASHT,NASHT)
      REAL*8 S1,S2,SM,SM1,SM2,GAA,GAB,GBA,GBB
      REAL*8 OVERLAP_RASSI,TMATEL,RED,FACT,CGCOEF,DCLEBS

#include "symmul.fh"
#include "stdalloc.fh"
#include "WrkSpc.fh"
      Real*8, Allocatable:: SPD1(:)

C Given two CI expansions, using a biorthonormal set of SD''s,
C calculate the following quantities:
C (1) The overlap
C (2) The spin-summed 1-particle transition density matrix
C (3) The WE-reduced transition spin density matrix
C in the biorthonormal active orbital basis.

      LORBTB=ISSTAB(3)
C Pick out nr of active orbitals from orbital table:
      NASORB=IWORK(LORBTB+3)

C Overlap:

      SIJ=0.0D00

      IF(MPLET1.EQ.MPLET2.AND.MSPROJ1.EQ.MSPROJ2) THEN
#ifdef _DMRG_
        if(.not.doDMRG)then
#endif
          SIJ=OVERLAP_RASSI(IFSBTAB1,IFSBTAB2,DET1,DET2)
#ifdef _DMRG_
        else
            sij = qcmaquis_mpssi_overlap(
     &         qcm_prefixes(job1),
     &         ist,
     &         qcm_prefixes(job2),
     &         jst,
     &         .true.)
        end if ! DMRG or not
#endif
      END IF ! mmplet and msproj check

C General 1-particle transition density matrix:
      NSPD1=NASORB**2
#ifdef _DMRG_
      if (.not.doDMRG) then
#endif
        Call mma_allocate(SPD1,nSPD1,Label='SPD1')
        SPD1(:)=0.0D0
#ifdef _DMRG_
      else
        ! For DMRG, we only need the AA and BB spin components
        ! Let's allocate two different arrays for them because that's
        ! easier for the new QCMaquis interface
        NSPD1=NASHT**2
        Call mma_allocate(TDMAA,nSPD1,Label='TDMAA')
        Call mma_allocate(TDMBB,nSPD1,Label='TDMBB')
        TDMAA(:)=0.0D0
        TDMBB(:)=0.0D0
      end if
#endif
      ISYOP = MUL(LSYM1,LSYM2)
      MS2OP = MSPROJ1-MSPROJ2

      IF(ABS(MS2OP).LE.2) THEN
#ifdef _DMRG_
        if(.not.doDMRG)then
#endif
          !> spind constructs the 1-particle transition density matrix
          !> output in SPD1
          !> main input: DET1 and DET2
          CALL SPIND(ISYOP,MS2OP,IWORK(LORBTB),ISSTAB,
     &               IFSBTAB1,IFSBTAB2,DET1,DET2,SPD1)

#ifdef _DMRG_
        else
          if(isyop /= 1)
     & stop 'MPS property density with spatial symm irrep > 1: FIXME!'
          ! calculate 1-TDMs: Must always be calculated with the higher multiplicity as <T|o|S>
          ! where T always has a higher multiplicity than S

          if (MPLET1.lt.MPLET2) then
            call qcmaquis_mpssi_get_onetdm_spin(
     &             qcm_prefixes(job2),
     &             jst,
     &             qcm_prefixes(job1),
     &             ist,
     &             TDMAA,
     &             TDMBB,
     &             NSPD1)
          else
            call qcmaquis_mpssi_get_onetdm_spin(
     &             qcm_prefixes(job1),
     &             ist,
     &             qcm_prefixes(job2),
     &             jst,
     &             TDMAA,
     &             TDMBB,
     &             NSPD1)
          end if

c         if(debug_dmrg_rassi_code)then
c           write(6,*) 'density for i, j',istate,jstate
c           write(6,*) 'dimension: ',nasorb**2, '--> #nact', nasorb
c           call pretty_print_util(SPD1,1,nasorb,1,nasorb,
c      &                           nasorb,nasorb,1,6)
c         end if
      end if
#endif
      END IF
C Create a scalar, and an WE-reduced spin, transition density matrix.
C The scalar matrix is simply the usual spin-summed density matrix.
C The WE-reduced matrix is the one defined through Wigner-Eckarts thm:
C  <A S1 M1 | T[K]_Q |B S2 M2>
C                   = FACT*CG(S2 M2 K Q;S1 M1)*<A S1 || T[K] ||B S2>
C i.e. with the usual Clebsch-Gordan factor, and a prefactor
C   FACT=(-1)**(MAX(S1,S2)-S1)/SQRT(2*S1+1)
      S1=DBLE(MPLET1-1)*0.5D0
      S2=DBLE(MPLET2-1)*0.5D0
      SM1=DBLE(MSPROJ1)*0.5D0
      SM2=DBLE(MSPROJ2)*0.5D0
      DO IORB=1,NASHT
       ISORB=2*IORB-1
       DO JORB=1,NASHT
        JSORB=2*JORB-1
#ifdef _DMRG_
        if (.not.doDMRG) then
#endif
          GAA=SPD1(0+ISORB+NASORB*(JSORB-1))
          GAB=SPD1(0+ISORB+NASORB*(JSORB  ))
          GBA=SPD1(1+ISORB+NASORB*(JSORB-1))
          GBB=SPD1(1+ISORB+NASORB*(JSORB  ))
#ifdef _DMRG_
        else
          GAA=TDMAA(JORB+NASHT*(IORB-1))
          GBB=TDMBB(JORB+NASHT*(IORB-1))
          ! transpose from row-major order,
          ! as it comes from C++ this way
          GAB=0.0D0
          GBA=0.0D0
        end if
#endif

C Position determined by active orbital index in external order:
        ITABS=MAPORB(ISORB)
        IUABS=MAPORB(JSORB)

c #ifdef _DMRG_
c         if(debug_dmrg_rassi_code)then
c           write(6,'(a,2i3,4f12.8)') ' i,j: GAA,GBB,GAB,GBA ==> ',
c      &                itabs,iuabs,gaa,gbb,gab,gba
c         end if
c #endif

        !> scalar TDM
        TDM1(ITABS,IUABS)=GAA+GBB
        !> spin TDM
        TSDM1(ITABS,IUABS)=GAA-GBB

C Clebsch-Gordan coefficient:
        SM=SM1-SM2
        FACT=1.0D0/SQRT(DBLE(MPLET1))
        IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
        CGCOEF=FACT*DCLEBS(S2,1.0D0,S1,SM2,SM,SM1)
C Spin tensor component matrix element
        IF(MSPROJ2.EQ.MSPROJ1+2) THEN
          TMATEL=SQRT(2.0D0)*GBA
        ELSE IF (MSPROJ2.EQ.MSPROJ1-2) THEN
          TMATEL=-SQRT(2.0D0)*GAB
        ELSE IF (MSPROJ2.EQ.MSPROJ1) THEN
          TMATEL=0.5D0*(GBB-GAA)
        END IF
C Thus obtain reduced matrix element from Wigner-Eckart theorem:
        RED=0.0D0
        IF(CGCOEF.NE.0.0D0) THEN
          RED=TMATEL/CGCOEF
        ELSE
          IF(ABS(TMATEL).GT.1.0D-12)THEN
            Call WarningMessage(1,'A possible bug was detected.')
            WRITE(6,*)' WARNING: Non-zero matrix element computed'
            WRITE(6,*)' which should be zero by spin symmetry!'
            WRITE(6,*)'              Spins S1, S2:',S1,S2
            WRITE(6,*)' Spin projections SM1, SM2:',SM1,SM2
            WRITE(6,*)'    Operator has S=1.0, SM:',SM
            WRITE(6,*)' Clebsch-Gordan:',CGCOEF
            WRITE(6,*)' Size is TMATEL=',TMATEL
          END IF
        END IF

        !> W-reduced TDM
        WTDM1(ITABS,IUABS)=RED

       END DO
      END DO

c #ifdef _DMRG_
c       if(debug_dmrg_rassi_code)then
c         !> debug print
c         write(6,*) '1-tdm density for i, j',istate,jstate
c         call pretty_print_util(tdm1,1,nasht,1,nasht,
c      &                         nasht,nasht,1,6)
c         write(6,*) '1-tdm sp-density for i, j',istate,jstate
c         call pretty_print_util(tsdm1,1,nasht,1,nasht,
c      &                         nasht,nasht,1,6)
c         write(6,*) 'w-reduced tdm for i, j',istate,jstate
c         call pretty_print_util(wtdm1,1,nasht,1,nasht,
c      &                         nasht,nasht,1,6)
c       end if
c #endif
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(ISTATE)
        CALL Unused_integer(JSTATE)
      END IF

#ifdef _DMRG_
      if (.not.doDMRG) then
#endif
      CALL mma_deallocate(SPD1)

#ifdef _DMRG_
      else
        call mma_deallocate(TDMAA)
        call mma_deallocate(TDMBB)
      end if
#endif
      RETURN
      END
