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

subroutine SG_Setup_RASSCF(DBG,SkipGUGA,initial_occ)

use Molcas, only: MxLev
use fciqmc, only: DoNECI
use fcidump, only: DumpOnly
use CC_CI_mod, only: Do_CC_CI
use gas_data, only: iDoGAS, NGAS, NGSSH
use rasscf_global, only: DoBlockDMRG, NSM
use general_data, only: iSpin, nActel, nConf, nElec3, nHole1, nRs1, nRs2, nRs3, nSym, STSYM
use general_data, only: CIS, EXS, SGS
use sguga, only: MKCOT, MKSGNUM, SG_Init_Simple
#ifdef _DMRG_
use rasscf_global, only: DoDMRG
use input_ras, only: Key
use stdalloc, only: mma_deallocate
#endif
use rasdef, only: nRas,nRasEl,nRsPrt, IFRAS
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(inout):: DBG,SkipGUGA
integer(kind=iwp), allocatable, optional, intent(inout) :: initial_occ(:,:)
integer(kind=iwp) :: IGAS, iq, ISYM, Level(MxLev), NLEV, NSTA, nRs1T
real(kind=wp) :: dum1, dum2, dum3, Eterna_1, Eterna_2

NLEV = 0
do IGAS=1,NGAS
  do ISYM=1,NSYM
    NSTA = NLEV+1
    NLEV = NLEV+NGSSH(IGAS,ISYM)
    NSM(NSTA:NLEV) = ISYM
  end do
end do
Level(1:MxLev)=[(iq,iq=1,MxLev)]

If (nHole1+nElec3/=0) Then
   IFRAS=1
   nRsPrt=3
   nRas(:,1)=nRs1(:)
   nRas(:,2)=nRs2(:)
   nRas(:,3)=nRs3(:)
   nRs1T=Sum(nRs1(1:nSym))
   nRasEl(1)=2*nRs1T-nHole1
   nRasEl(2)=nActel-nElec3
   nRasEl(3)=nActel
Else
   IFRAS=0
   nRsPrt=1
   nRas(:,1)=nRs2(:)
   nRasEl(1)=nActel
End If
Do ISYM=1,NSYM
   IF (SUM(nRAS(ISYM,1:nRsPrt))/=0) IFRAS=IFRAS+1
END DO
SGS%IFRAS=IFRAS

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
#   else
#   include "macros.fh"
    unused_opt(initial_occ)
#   endif
      call Timing(Eterna_1,dum1,dum2,dum3)
      if (DBG) write(u6,*) ' Call SG_Init_Simple'
      call SG_Init_Simple(nSym,nActEl,iSpin,SGS,CIS,EXS,                &
                          xLevel=Level,xL2Act=Level,xNLEV=NLEV,xNSM=NSM)

      if (SGS%NVERT0 == 0) then
         CIS%NCSF(STSYM) = 0
      else
         if (doBlockDMRG) then
            CIS%NCSF(STSYM) = 1
        else

            ! FORM VARIOUS OFFSET TABLES:
            ! NOTE: NIPWLK AND DOWNWLK ARE THE NUMER OF INTEGER WORDS USED
            !       TO STORE THE UPPER AND LOWER WALKS IN PACKED FORM.

            ! CONSTRUCT THE CASE LIST

          call MKCOT(SGS,CIS)

            ! SET UP ENUMERATION TABLES

            call MKSGNUM(STSYM,SGS,CIS,EXS)

            if (NActEl == 0) CIS%NCSF(STSYM) = 1

            !     (SGS%IFRAS-1) IS THE NUMBER OF SYMMETRIES CONTAINING ACTIVE ORBITALS
            !     IF THIS IS GREATER THAN 1 ORBITAL REORDERING INTEGRALS IS REQUIRED
            !     SET UP THE REINDEXING TABLE
        end if
      end if
      call SETSXCI()
      NCONF = CIS%NCSF(STSYM)

      call Timing(Eterna_2,dum1,dum2,dum3)
#   ifdef _DMRG_
    end if
#   endif
  end if
end if

end subroutine sg_setup_rasscf
