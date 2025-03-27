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

subroutine H0MAT_MCLR(H0,ISBDET,ISBCNF,MXP1,MXP2,MXQ,NOCOB,NPRCIV,NOCSF,IREFSM,IDC,PSSIGN,ECORE,VEC1,VEC2,H0SCR,iH0SCR,ieaw)
! Obtain preconditioner space corresponding to internal space INTSPC
! Obtain Hamiltonian matrices corresponding to this subspace
!
! Construct Preconditioner blocks of Hamiltonian matrix
!
! ======
!.Output
! ======
!
! CSF : NP1CSF,NP2CSF,NQCSF : Number of CSFs in the 3 primary subspaces
!
! NPRCIV : Number of parameters in preconditioner space

use Str_Info, only: CNSM, DFTP, DTOC
use MCLR_Data, only: NAELCI, NBELCI
use MCLR_Data, only: NCNASM
use Constants, only: Zero

implicit none
! Offsets for CSF information
real*8 H0(*)
integer ISBDET(*)
integer ISBCNF(*)
integer MXP1, MXP2, MXQ, NOCOB, NPRCIV, NOCSF, IREFSM, IDC
real*8 PSSIGN, ECORE
real*8 vec1(*), vec2(*), h0scr(*)
integer ih0scr(*)
integer ieaw
!integer, allocatable :: IOCOC(:)
integer intspc, NAEL, NBEL, ICOMBI, IPWAY, NINOB, NP1CNF, NP1CSF, NP2CNF, NP2CSF, NPCNF, NQCNF, NQCSF
real*8 PSIGN

! Info on actual internal subspace
intspc = 1
!luhdia = 0
!IATP = IASTFI(INTSPC)
!IBTP = IBSTFI(INTSPC)
!MNR1 = MNR1IC(INTSPC)
!MXR1 = MXR1IC(INTSPC)
!MNR3 = MNR3IC(INTSPC)
!MXR3 = MXR3IC(INTSPC)
NAEL = NAELCI(INTSPC)
NBEL = NBELCI(INTSPC)

!NOCTPA = NOCTYP(IATP)
!NOCTPB = NOCTYP(IBTP)
! Allowed combination of alpha and beta strings
!call mma_allocate(IOCOC,NOCTPA*NOCTPB,Label='IOCOC')
!call IAIBCM_MCLR(MNR1,MXR3,NOCTPA,NOCTPB,Str(IATP)%EL1,Str(IATP)%EL3,Str(IBTP)%EL1,Str(IBTP)%EL3,IOCOC)

if (IDC == 1) then
  ICOMBI = 0
  PSIGN = Zero
else
  PSIGN = PSSIGN
  ICOMBI = 1
end if

!if (NOCSF /= 0) then
!  ! Combinations expansion, PQ preconditioner
!
!  IHAMSM = 1
!  IWAY = 1
!  ! strings are unsigned
!  ISTSGN = 0
!  call H0SD(LUHDIA,LBLK,VEC1,IWAY,NSBDET,NAEL,NBEL,ISMOST(1,IREFSM),IOCOC,IHAMSM,H0,NOCOB,0,ECORE,ICOMBI,PSIGN,NPRCIV,SBEVC, &
!            SBEVL,1,NCIVAR,ISBDET,ISBIA,ISBIB,MXP1,MXP2,MXQ,MP1CSF,MP2CSF,MQCSF,Str(IATP)%OCSTR, Str(IBTP)%OCSTR,ISTSGN,IDUMMY, &
!            IDUMMY,INTSPC)
!else if (NOCSF == 0) then
!  CSF basis,PQ preconditioner
IPWAY = 2
call H0CSF(H0,ISBDET,ISBCNF,MXP1,MXP2,MXQ,DTOC,DFTP,CNSM(ieaw)%ICONF,IREFSM,ECORE,NINOB,NOCOB,H0SCR,iH0SCR,NCNASM(IREFSM), &
           NAEL+NBEL,NAEL,NBEL,IPWAY,NP1CSF,NP1CNF,NP2CSF,NP2CNF,NQCSF,NQCNF,NPRCIV,NPCNF,VEC1,VEC2,INTSPC,ICOMBI,PSIGN)
!end if

!call mma_deallocate(IOCOC)

if (.false.) call Unused_integer(NOCSF)

end subroutine H0MAT_MCLR
