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
     &    lLROOT,job1,job2,ist,jst)

      !> module dependencies
#ifdef _DMRG_
      use qcmaquis_interface_cfg
      use qcmaquis_interface_wrapper
      use qcmaquis_interface_utility_routines, only:
     &    pretty_print_util
      use qcmaquis_info
#endif

      IMPLICIT NONE
      INTEGER LSYM1,MPLET1,MSPROJ1,LSYM2,MPLET2,MSPROJ2,lLROOT
      INTEGER IFSBTAB1(*),IFSBTAB2(*),ISSTAB(*),MAPORB(*)
      INTEGER IORB,ISORB,ISYOP,ITABS,IUABS,JORB,JSORB,LORBTB
      INTEGER LSPD1,MS2OP,NASHT,NASORB,NSPD1
      INTEGER, INTENT(IN) :: ISTATE, JSTATE,job1,job2,ist,jst
      REAL*8 DET1(*),DET2(*)
      REAL*8 SIJ,TDM1(NASHT,NASHT),TSDM1(NASHT,NASHT),WTDM1(NASHT,NASHT)
      REAL*8 S1,S2,SM,SM1,SM2,GAA,GAB,GBA,GBB
      REAL*8 OVERLAP_RASSI,TMATEL,RED,FACT,CGCOEF,DCLEBS
#ifdef _DMRG_
      LOGICAL :: debug_dmrg_rassi_code = .false.
#endif

#ifdef _DMRG_
      ! strings for conversion of the qcmaquis h5 checkpoint names from 2u1 to su2u1
      character(len=3) :: mplet1s, msproj1s
      ! new checkpoint names
      character(len=2300) :: checkpoint1_2u1,checkpoint2_2u1
#endif

#include "symmul.fh"
#include "WrkSpc.fh"

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
          if (doMPSSICheckpoints) then
            if (dmrg_external%MPSrotated) then
              write(mplet1s,'(I3)')  MPLET1-1
              write(msproj1s,'(I3)')  MSPROJ1

              checkpoint1_2u1 = qcm_group_names(job1)%states(ist)
     &(1:len_trim(qcm_group_names(job1)%states(ist))-3)
     &//"."//trim(adjustl(mplet1s))//"."//trim(adjustl(msproj1s))//".h5"
              checkpoint2_2u1 = qcm_group_names(job2)%states(jst)
     &(1:len_trim(qcm_group_names(job2)%states(jst))-3)
     &//"."//trim(adjustl(mplet1s))//"."//trim(adjustl(msproj1s))//".h5"

              call dmrg_interface_ctl(
     &                            task        = 'overlapU',
     &                            energy      = sij,
     &                            checkpoint1 = checkpoint1_2u1,
     &                            checkpoint2 = checkpoint2_2u1
     &                           )
            else
              call dmrg_interface_ctl(
     &                            task   = 'overlap ',
     &                            energy = sij,
     &                            checkpoint1 =
     &                            qcm_group_names(job1)%states(ist),
     &                            checkpoint2 =
     &                            qcm_group_names(job2)%states(jst)
     &                           )
            end if
          else
          ! Leon: TODO: Add possibility to calculate overlap of rotated MPS without using checkpoint names
            call dmrg_interface_ctl(
     &                            task   = 'overlap ',
     &                            energy = sij,
     &                            state  = iWork(lLROOT+istate-1),
     &                            stateL = iWork(lLROOT+jstate-1)
     &                           )
          end if

        end if ! DMRG or not
#endif

        if(debug_dmrg_rassi_code)then
          write(6,*) 'overlap sij is for i,j',istate,jstate,sij
          call flush(6)
        end if

      END IF ! mmplet and msproj check

C General 1-particle transition density matrix:
      NSPD1=NASORB**2
      CALL GETMEM('SPD1','Allo','Real',LSPD1,NSPD1)
      CALL DCOPY_(NSPD1,0.0D0,0,WORK(LSPD1),1)
      ISYOP = MUL(LSYM1,LSYM2)
      MS2OP = MSPROJ1-MSPROJ2

      IF(ABS(MS2OP).LE.2) THEN
#ifdef _DMRG_
        if(.not.doDMRG)then
#endif
          !> spind constructs the 1-particle transition density matrix
          !> output in WORK(LSPD1)
          !> main input: DET1 and DET2
          CALL SPIND(ISYOP,MS2OP,IWORK(LORBTB),ISSTAB,
     &               IFSBTAB1,IFSBTAB2,DET1,DET2,WORK(LSPD1))

#ifdef _DMRG_
        else
          if(isyop /= 1)
     & stop 'MPS property density with spatial symm irrep > 1: FIXME!'
          if (doMPSSICheckpoints) then
            call dmrg_interface_ctl(
     &                            task       = 'imp rdmY',
     &                            x1         = work(lspd1),
     &                            ndim       = nasorb,
     &                            checkpoint1=
     &                            qcm_group_names(job1)%states(ist),
     &                            checkpoint2=
     &                            qcm_group_names(job2)%states(jst),
     &                            msproj     = msproj1,
     &                            msprojL    = msproj2,
     &                            multiplet  = MPLET1-1, ! (MPLET1 == 2*S+1) and we need 2*S
     &                            multipletL = MPLET2-1, ! (MPLET2 == 2*S+1) and we need 2*S
     &                            rdm1       = .true.,
     &                            rdm2       = .false.
     &                           )
          else
            call dmrg_interface_ctl(
     &                            task       = 'imp rdmY',
     &                            x1         = work(lspd1),
     &                            ndim       = nasorb,
     &                            state      = iWork(lLROOT+istate-1),
     &                            stateL     = iWork(lLROOT+jstate-1),
     &                            msproj     = msproj1,
     &                            msprojL    = msproj2,
     &                            multiplet  = MPLET1-1, ! (MPLET1 == 2*S+1) and we need 2*S
     &                            multipletL = MPLET2-1, ! (MPLET2 == 2*S+1) and we need 2*S
     &                            rdm1       = .true.,
     &                            rdm2       = .false.
     &                           )
          end if
        end if
#endif

#ifdef _DMRG_
        if(debug_dmrg_rassi_code)then
          write(6,*) 'density for i, j',istate,jstate
          write(6,*) 'dimension: ',nasorb**2, '--> #nact', nasorb
          call pretty_print_util(WORK(LSPD1),1,nasorb,1,nasorb,
     &                           nasorb,nasorb,1,6)
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
        GAA=WORK(LSPD1-1+ISORB+NASORB*(JSORB-1))
        GAB=WORK(LSPD1-1+ISORB+NASORB*(JSORB  ))
        GBA=WORK(LSPD1  +ISORB+NASORB*(JSORB-1))
        GBB=WORK(LSPD1  +ISORB+NASORB*(JSORB  ))

C Position determined by active orbital index in external order:
        ITABS=MAPORB(ISORB)
        IUABS=MAPORB(JSORB)

#ifdef _DMRG_
        if(debug_dmrg_rassi_code)then
          write(6,'(a,2i3,4f12.8)') ' i,j: GAA,GBB,GAB,GBA ==> ',
     &                itabs,iuabs,gaa,gbb,gab,gba
        end if
#endif

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

#ifdef _DMRG_
      if(debug_dmrg_rassi_code)then
        !> debug print
        write(6,*) '1-tdm density for i, j',istate,jstate
        call pretty_print_util(tdm1,1,nasht,1,nasht,
     &                         nasht,nasht,1,6)
        write(6,*) '1-tdm sp-density for i, j',istate,jstate
        call pretty_print_util(tsdm1,1,nasht,1,nasht,
     &                         nasht,nasht,1,6)
        write(6,*) 'w-reduced tdm for i, j',istate,jstate
        call pretty_print_util(wtdm1,1,nasht,1,nasht,
     &                         nasht,nasht,1,6)
      end if
#else
c Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_integer(ISTATE)
        CALL Unused_integer(JSTATE)
      END IF
#endif

      CALL GETMEM('SPD1','Free','Real',LSPD1,NSPD1)

      RETURN
      END
