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

#ifdef _DMRG_
subroutine SG_Setup_RASSCF(DBG,SkipGUGA,initial_occ)
#else
subroutine SG_Setup_RASSCF(DBG,SkipGUGA)
#endif

use Molcas, only: MxLev
use fciqmc, only: DoNECI
use fcidump, only: DumpOnly
use CC_CI_mod, only: Do_CC_CI
use gas_data, only: iDoGAS
use rasscf_global, only: DoBlockDMRG
use general_data, only: nSym, nActel, iSpin, nHole1, nElec3, nRs1, nRs2, nRs3, STSYM, nConf
use general_data, only: CIS, EXS, SGS
#ifdef _DMRG_
use rasscf_global, only: DoDMRG
use input_ras, only: Key
use stdalloc, only: mma_deallocate
#endif
use gas_data, only: NGAS, NGSSH
use rasscf_global, only: NSM
use sguga, only: SG_Init_Simple, MKCOT, MKCLIST, MKSGNUM
use definitions, only: u6, iwp, wp

#ifdef _DMRG_
integer(kind=iwp), allocatable, intent(inout) :: initial_occ(:,:)
#endif

logical(kind=iwp), intent(inout):: DBG,SkipGUGA
real(kind=wp) Eterna_1, Eterna_2, dum1, dum2, dum3
integer(kind=iwp) :: IGAS, ISYM, NLEV, NSTA, iq, Level(MxLev)

NLEV = 0
do IGAS=1,NGAS
  do ISYM=1,NSYM
    NSTA = NLEV+1
    NLEV = NLEV+NGSSH(IGAS,ISYM)
    NSM(NSTA:NLEV) = ISYM
  end do
end do
Level(1:MxLev)=[(iq,iq=1,MxLev)]

! Construct the Guga tables

if (.not. (DoNECI .or. Do_CC_CI .or. DumpOnly .or. SkipGUGA)) then
  ! right now skip most part of gugactl for GAS, but only call mkism.
  if (.not. iDoGas) then
    ! DMRG calculation no need the SG_Init_RASSCF subroutine
#   ifdef _DMRG_
    if (Key('DMRG') .or. doDMRG) then
      call mma_deallocate(initial_occ)
      SkipGUGA = .true.
    else
#   endif
      call Timing(Eterna_1,dum1,dum2,dum3)
      if (DBG) write(u6,*) ' Call SG_Init_Simple'
      Call SG_Init_Simple(nSym,nActEl,iSpin,SGS,CIS,                 &
                          EXS,nHole1,nElec3,nRs1,nRs2,nRs3,          &
                          xLevel=Level,xL2Act=Level,                 &
                          xNLEV=NLEV,xNSM=NSM)

      if (SGS%NVERT0 == 0) then
         CIS%NCSF(STSYM) = 0
      Else
         if (doBlockDMRG) then
            CIS%NCSF(STSYM) = 1
         Else

            ! FORM VARIOUS OFFSET TABLES:
            ! NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
            !       TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.

            call MKCOT(SGS,CIS)

            ! CONSTRUCT THE CASE LIST

            call MKCLIST(SGS,CIS)

            ! SET UP ENUMERATION TABLES

            call MKSGNUM(STSYM,SGS,CIS,EXS)

            if (NActEl == 0) CIS%NCSF(STSYM) = 1

            !     (SGS%IFRAS-1) IS THE NUMBER OF SYMMETRIES CONTAINING ACTIVE ORBITALS
            !     IF THIS IS GREATER THAN 1 ORBITAL REORDERING INTEGRALS IS REQUIRED
            !     SET UP THE REINDEXING TABLE
         End If
      End If
      call SETSXCI()
      NCONF = CIS%NCSF(STSYM)

      call Timing(Eterna_2,dum1,dum2,dum3)
#   ifdef _DMRG_
    end if
#   endif
  end if
end if

end subroutine sg_setup_rasscf
