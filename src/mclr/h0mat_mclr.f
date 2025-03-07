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
      SUBROUTINE H0MAT_MCLR(H0,ISBDET, ISBCNF,
     &                 MXP1,MXP2,MXQ,NOCOB,
     &                 NPRCIV,NOCSF,IREFSM,
     &                 IDC,PSSIGN,ECORE,
     &                 VEC1,VEC2,H0SCR,iH0SCR,ieaw)
      Use Str_Info, only: CNSM,DFTP,DTOC
      use MCLR_Data, only: NAELCI, NBELCI
      use MCLR_Data, only: NCNASM
* Obtain preconditioner space corresponding to internal space INTSPC
* Obtain Hamiltonian matrices corresponding to this subspace
*
* Construct Preconditioner blocks of Hamiltonian matrix
*
* ======
*.Output
* ======
*
* CSF : NP1CSF,NP2CSF,NQCSF : Number of CSFs in the 3 primary subspaces
*
* NPRCIV : Number of parameters in preconditioner space
*
      IMPLICIT None
*. Offsets for CSF information
*
      Real*8 H0(*)
      Integer ISBDET(*)
      Integer ISBCNF(*)
      Integer MXP1,MXP2,MXQ,NOCOB,NPRCIV,NOCSF,IREFSM,IDC
      Real*8 PSSIGN,ECORE
      Real*8 vec1(*),vec2(*),h0scr(*)
      Integer ih0scr(*)
      Integer ieaw
*     Integer, Allocatable:: IOCOC(:)
      Integer iprt,intspc,NAEL,NBEL,ICOMBI,IPWAY,NINOB,NP1CNF,NP1CSF,
     &        NP2CNF,NP2CSF,NPCNF,NQCNF,NQCSF
      Real*8  PSIGN
*
* Info on actual internal subspace
      iprt=100
      intspc=1
*     luhdia=0
*     IATP = IASTFI(INTSPC)
*     IBTP = IBSTFI(INTSPC)
*     MNR1 = MNR1IC(INTSPC)
*     MXR1 = MXR1IC(INTSPC)
*     MNR3 = MNR3IC(INTSPC)
*     MXR3 = MXR3IC(INTSPC)
      NAEL = NAELCI(INTSPC)
      NBEL = NBELCI(INTSPC)
*
*     NOCTPA = NOCTYP(IATP)
*     NOCTPB = NOCTYP(IBTP)
*. Allowed combination of alpha and beta strings
*     Call mma_allocate(IOCOC,NOCTPA*NOCTPB,Label='IOCOC')
*     CALL IAIBCM_MCLR(MNR1,MXR3,NOCTPA,NOCTPB,
*    &            Str(IATP)%EL1,Str(IATP)%EL3,
*    &            Str(IBTP)%EL1,Str(IBTP)%EL3,
*    &            IOCOC,NTEST)
*
      IF(IDC.EQ.1) THEN
        ICOMBI = 0
        PSIGN = 0.0D0
      ELSE
        PSIGN = PSSIGN
        ICOMBI = 1
      END IF
*
*     IF( NOCSF .NE. 0) THEN
*.Combinations expansion, PQ preconditioner
*
*       IHAMSM = 1
*       IWAY = 1
* strings are unsigned
*       ISTSGN = 0
*       CALL H0SD(LUHDIA,LBLK,VEC1,IWAY,NSBDET,NAEL,NBEL,
*    &            ISMOST(1,IREFSM),IOCOC,
*    &            IHAMSM,H0,NOCOB,0,
*    &            ECORE,ICOMBI,PSIGN,NPRCIV,SBEVC,
*    &            SBEVL,1,NCIVAR,ISBDET,ISBIA,ISBIB,
*    &            MXP1,MXP2,MXQ,
*    &            MP1CSF,MP2CSF,MQCSF,
*    &            Str(IATP)%OCSTR, Str(IBTP)%OCSTR,
*    &            ISTSGN,IDUMMY,IDUMMY,
*    &            INTSPC,IPRT)
*     ELSE IF (NOCSF.EQ.0) THEN
*.CSF basis,PQ preconditioner
        IPWAY = 2
        CALL H0CSF(H0,ISBDET,ISBCNF,
     &       MXP1,MXP2,MXQ,
     &       DTOC,
     &       DFTP,CNSM(ieaw)%ICONF,
     &       IREFSM,ECORE,NINOB,NOCOB,
     &       H0SCR,iH0SCR,NCNASM(IREFSM),
     &       NAEL+NBEL,NAEL,NBEL,IPWAY,
     &       NP1CSF,NP1CNF,NP2CSF,NP2CNF,NQCSF,NQCNF,
     &       NPRCIV,NPCNF,
     &       VEC1,VEC2,IPRT,INTSPC,ICOMBI,PSIGN)
*     END IF
*
*     Call mma_deallocate(IOCOC)

      IF (.FALSE.) CALL Unused_integer(NOCSF)
      END SUBROUTINE H0MAT_MCLR
