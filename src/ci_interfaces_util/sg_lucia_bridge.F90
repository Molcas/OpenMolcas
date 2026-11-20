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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> The bridge between the CSF basis of Lucia and the split-graph (SGUGA)
!> basis: reordering of CI vectors and the walks of the Lucia CSFs. It uses
!> Lucia's configuration tables and the SGUGA graph, so it belongs to the
!> CI interface (the gateway), not to either environment.
module SG_Lucia_Bridge

use sguga, only: CIS, EXS, MkCOT, MkSgNum, SG_NUM, SG_Num_Setup, SG_PHASE, SGS, StepVec
use Definitions, only: wp, iwp, u6

implicit none
private

public :: SG_Lucia_Walk, SG_ReOrd, SG_ReOrd_Orig, SG_Verify_Lucia_Walk

contains

! This is temporarily for the MCLR code
subroutine SG_ReOrd_orig(iState,IREFSM,IMODE,nConf,CIOLD,CINEW)
!***********************************************************************
!                                                                      *
!     Rearrange CI-vectors                                             *
!     iMode=0 --> from SGA to split graph GUGA order                   *
!     iMode=1 --> from split graph GUGA to SGA order                   *
!                                                                      *
!     calling arguments:                                               *
!     iRefSm  : integer                                                *
!               state symmetry                                         *
!     iMode   : integer                                                *
!               switch selecting reordering mode (see above)           *
!     nSm     : array of integer                                       *
!               symmetry per active orbital                            *
!     CIold   : array of real                                          *
!               incoming CI vector                                     *
!     CInew   : array of real                                          *
!               outgoing CI vector                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher and J. Olsen                                      *
!     University of Lund, Sweden, 1990                                 *
!***********************************************************************

use lucia_data, only: MINOP, NCNFTP, NCSFTP, NTYP
use Lucia_data, only: CFTP, CONF_Occ
use Molcas, only: MxAct
use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iState, IREFSM, IMODE, nConf
real(kind=wp), intent(in) :: CIOLD(nConf)
real(kind=wp), intent(out) :: CINEW(nConf)
integer(kind=iwp) :: IC, ICL, ICNBS, ICNBS0, ICSBAS, ICSFJP, IIBCL, IIBOP, IICSF, IOPEN, IP, IPBAS, ISG, ITYP, &
                     IWALK(mxAct), JOCC, KCNF(MxAct), KOCC, KORB
real(kind=wp) :: Fact
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif

call SG_Num_Setup(iState,IREFSM)


ICSFJP = 0
ICNBS0 = 0 ! dummy initialize
IPBAS = 0 ! dummy initialize
! LOOP OVER CONFIGURATIONS TYPES
do ITYP=1,NTYP
  IOPEN = ITYP+MINOP-1
  ICL = (SGS(iState)%nActEl-IOPEN)/2
  ! BASE ADDRESS FOR CONFIGURATION OF THIS TYPE
  if (ITYP == 1) then
    ICNBS0 = 1
  else
    ICNBS0 = ICNBS0+NCNFTP(ITYP-1,IREFSM)*(SGS(iState)%nActEl+IOPEN-1)/2
  end if
  ! BASE ADDRESS FOR PROTOTYPE SPIN COUPLINGS
  if (ITYP == 1) then
    IPBAS = 1
  else
    IPBAS = IPBAS+NCSFTP(ITYP-1)*(IOPEN-1)
  end if

  ! LOOP OVER NUMBER OF CONFIGURATIONS OF TYPE ITYP AND PROTOTYPE
  ! SPIN COUPLINGS

  do IC=1,NCNFTP(ITYP,IREFSM)
    ICNBS = ICNBS0+(IC-1)*(IOPEN+ICL)
    do IICSF=1,NCSFTP(ITYP)
      ICSFJP = ICSFJP+1
      ICSBAS = IPBAS+(IICSF-1)*IOPEN
      KCNF(:)=0
      ! Obtain configuration in standard RASSCF form
      IIBOP = 1
      IIBCL = 1
      JOCC = ICL+IOPEN
      do KOCC=0,JOCC-1
        KORB = Conf_Occ(IREFSM)%A(ICNBS+KOCC)
        if (KORB < 0) then
          ! Doubly occupied orbital
          KCNF(IIBCL) = abs(KORB)
          IIBCL = IIBCL+1
        else
          ! Singly occupied orbital
          KCNF(ICL+IIBOP) = KORB
          IIBOP = IIBOP+1
        end if
      end do

      ! COMPUTE STEP VECTOR
      call STEPVEC(KCNF(1:ICL),KCNF(ICL+1),ICL,IOPEN,CFTP(ICSBAS),SGS(iState)%nLev,IWALK)

      ! GET SPLIT GRAPH ORDERING NUMBER
      ISG = SG_NUM(SGS(iState),EXS(istate),IWALK)
      ! GET PHASE PHASE FACTOR
      IP = SG_PHASE(SGS(istate),IWALK)
      Fact = merge(-One,One,IP < 0)
      if (IMODE == 0) then
        CINEW(ISG) = Fact * CIOLD(ICSFJP)
      else
        CINEW(ICSFJP) = Fact * CIOLD(ISG)
      end if
    end do
  end do
end do

#ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) ' OLD CI-VECTOR IN SUBROUTINE REORD (MAX. 200 ELEMENTS)'
write(u6,'(10F12.8)') (CIOLD(I),I=1,min(200,ICSFJP))
  write(u6,*) ' NEW CI-VECTOR IN SUBROUTINE REORD (MAX. 200 ELEMENTS)'
write(u6,'(10F12.8)') (CINEW(I),I=1,min(200,ICSFJP))
  write(u6,*)
#endif

end subroutine SG_Reord_Orig

subroutine SG_ReOrd(iState,IREFSM,IMODE,nConf,CIOLD,CINEW)

use Molcas, only: MxAct
use Constants, only: One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none

integer(kind=iwp), intent(in)  :: iState
integer(kind=iwp), intent(in)  :: IREFSM
integer(kind=iwp), intent(in)  :: IMODE
integer(kind=iwp), intent(in)  :: nConf

real(kind=wp), intent(in)  :: CIOLD(nConf)
real(kind=wp), intent(out) :: CINEW(nConf)

integer(kind=iwp) :: IC
integer(kind=iwp) :: ICL
integer(kind=iwp) :: ICSFJP
integer(kind=iwp) :: IICSF
integer(kind=iwp) :: IOPEN
integer(kind=iwp) :: IP
integer(kind=iwp) :: ISG
integer(kind=iwp) :: MinOpen, nOC
integer(kind=iwp) :: MaxOpen

integer(kind=iwp) :: IWALK(MxAct)
integer(kind=iwp) :: KCNF(MxAct)

real(kind=wp) :: Fact

#ifdef _DEBUGPRINT_
integer(kind=iwp) :: I
#endif
If (SGS(iState)%nRsPrt /= 1) Then
   Call SG_ReOrd_orig(iState,IREFSM,IMODE,nConf,CIOLD,CINEW)
   Return
endif

if (.not. allocated(CIS(iState)%ICASE)) call MkCOT(SGS(iState),CIS(iState))
call MkSgNum(IREFSM,SGS(iState),CIS(iState),EXS(iState))

ICSFJP = 0

MinOpen = SGS(iState)%iSpin - 1
MaxOpen = min(SGS(iState)%nActEl, &
              2*SGS(iState)%nLev-SGS(iState)%nActEl)
if (mod(MaxOpen-MinOpen,2) /= 0) MaxOpen = MaxOpen - 1

! loop over open-shell sectors. That is, run over all orbital configurations
! that have the same number of open-shells
do iOpen = MinOpen, MaxOpen, 2

   ICL = (SGS(iState)%nActEl-iOpen)/2

   nOC = SGS(iState)%CDCInfo(iOpen,SGS(iState)%iSpin)%nOCSym(IREFSM)

!  Run over all prototype configurations of the sector
   do IC = 1, nOC

      call GetOccupation_SGUGA(iState,iOpen,SGS(iState)%iSpin,IREFSM,IC,KCNF)

!     Loop over all orbital combinations for this particular sector.
      do IICSF = 1, SGS(iState)%CDCInfo(iOpen,SGS(iState)%iSpin)%nCSF

         ICSFJP = ICSFJP + 1

         call BuildWalk(iState,ICL,iOpen,IICSF,KCNF,IWALK)

#ifdef _SANTITY_CHECK_
         if (any(IWALK(1:SGS(iState)%nLev) < 0) .or. &
             any(IWALK(1:SGS(iState)%nLev) > 3)) then

            write(u6,*) 'INVALID IWALK'
            write(u6,'(20I4)') IWALK(1:SGS(iState)%nLev)
            call Abend()
         endif
#endif

         ISG = SG_NUM(SGS(iState),EXS(iState),IWALK)
         IP = SG_PHASE(SGS(iState),IWALK)
         Fact = One
         if (IP < 0) Fact = -One
         if (IMODE == 0) then
            CINEW(ISG) = Fact*CIOLD(ICSFJP)
         else
            CINEW(ICSFJP) = Fact*CIOLD(ISG)
         endif
      enddo
   enddo
enddo

#ifdef _DEBUGPRINT_

write(u6,*)
write(u6,*) ' OLD CI-VECTOR IN SUBROUTINE REORD (MAX. 200 ELEMENTS)'
write(u6,'(10F12.8)') (CIOLD(I),I=1,min(200,ICSFJP))

write(u6,*)
write(u6,*) ' NEW CI-VECTOR IN SUBROUTINE REORD (MAX. 200 ELEMENTS)'
write(u6,'(10F12.8)') (CINEW(I),I=1,min(200,ICSFJP))
write(u6,*)

#endif

end subroutine SG_ReOrd

!> Step vector IWALK of the CSF ICSF of the Lucia (SGA) CSF basis of
!> symmetry IREFSM and the phase IP relating it to the GUGA CSF, as in
!> SG_ReOrd: CINEW(SG_NUM(IWALK)) = IP*CIOLD(ICSF) for IMODE = 0.
!> Only for a single RAS partition (nRsPrt = 1, the case of SG_ReOrd proper).
subroutine SG_Lucia_Walk(iState,IREFSM,ICSF,IWALK,IP)

use Molcas, only: MxAct

integer(kind=iwp), intent(in) :: iState, IREFSM, ICSF
integer(kind=iwp), intent(out) :: IWALK(MxAct), IP
integer(kind=iwp) :: iOpen, ICL, MinOpen, MaxOpen, nOC, nCSF, iOff, IC, IICSF
integer(kind=iwp) :: KCNF(MxAct)

if (SGS(iState)%nRsPrt /= 1) then
  call SG_Lucia_Walk_CFTP(iState,IREFSM,ICSF,IWALK,IP)
  return
end if

MinOpen = SGS(iState)%iSpin-1
MaxOpen = min(SGS(iState)%nActEl,2*SGS(iState)%nLev-SGS(iState)%nActEl)
if (mod(MaxOpen-MinOpen,2) /= 0) MaxOpen = MaxOpen-1

iOff = 0
do iOpen=MinOpen,MaxOpen,2
  nOC = SGS(iState)%CDCInfo(iOpen,SGS(iState)%iSpin)%nOCSym(IREFSM)
  nCSF = SGS(iState)%CDCInfo(iOpen,SGS(iState)%iSpin)%nCSF
  if (ICSF <= iOff+nOC*nCSF) then
    ICL = (SGS(iState)%nActEl-iOpen)/2
    IC = (ICSF-iOff-1)/nCSF+1
    IICSF = ICSF-iOff-(IC-1)*nCSF
    call GetOccupation_SGUGA(iState,iOpen,SGS(iState)%iSpin,IREFSM,IC,KCNF)
    call BuildWalk(iState,ICL,iOpen,IICSF,KCNF,IWALK)
    IP = SG_PHASE(SGS(iState),IWALK)
    return
  end if
  iOff = iOff+nOC*nCSF
end do

write(u6,*) 'SG_Lucia_Walk: CSF index out of range:',ICSF
call Abend()

end subroutine SG_Lucia_Walk

!> SG_Lucia_Walk from Lucia's own configuration list (Conf_Occ) and prototype
!> spin couplings (CFTP), as in SG_ReOrd_orig; valid for any RAS partitioning
!> (used for nRsPrt /= 1, where the CDCInfo tables are not built).
subroutine SG_Lucia_Walk_CFTP(iState,IREFSM,ICSF,IWALK,IP)

use lucia_data, only: MINOP, NCNFTP, NCSFTP, NTYP
use Lucia_data, only: CFTP, CONF_Occ
use Molcas, only: MxAct

integer(kind=iwp), intent(in) :: iState, IREFSM, ICSF
integer(kind=iwp), intent(out) :: IWALK(MxAct), IP
integer(kind=iwp) :: IC, ICL, ICNBS, ICNBS0, ICSBAS, IIBCL, IIBOP, IICSF, IOFF, IOPEN, IPBAS, ITYP, JOCC, KCNF(MxAct), KOCC, KORB, &
                     nBlk

ICNBS0 = 1
IPBAS = 1
IOFF = 0
do ITYP=1,NTYP
  IOPEN = ITYP+MINOP-1
  ICL = (SGS(iState)%nActEl-IOPEN)/2
  if (ITYP > 1) then
    ICNBS0 = ICNBS0+NCNFTP(ITYP-1,IREFSM)*(SGS(iState)%nActEl+IOPEN-1)/2
    IPBAS = IPBAS+NCSFTP(ITYP-1)*(IOPEN-1)
  end if
  nBlk = NCNFTP(ITYP,IREFSM)*NCSFTP(ITYP)
  if (ICSF <= IOFF+nBlk) then
    IC = (ICSF-IOFF-1)/NCSFTP(ITYP)+1
    IICSF = ICSF-IOFF-(IC-1)*NCSFTP(ITYP)
    ICNBS = ICNBS0+(IC-1)*(IOPEN+ICL)
    ICSBAS = IPBAS+(IICSF-1)*IOPEN
    ! configuration in standard RASSCF form: closed orbitals, then open ones
    KCNF(:) = 0
    IIBOP = 1
    IIBCL = 1
    JOCC = ICL+IOPEN
    do KOCC=0,JOCC-1
      KORB = Conf_Occ(IREFSM)%A(ICNBS+KOCC)
      if (KORB < 0) then
        KCNF(IIBCL) = abs(KORB)
        IIBCL = IIBCL+1
      else
        KCNF(ICL+IIBOP) = KORB
        IIBOP = IIBOP+1
      end if
    end do
    call STEPVEC(KCNF(1:ICL),KCNF(ICL+1),ICL,IOPEN,CFTP(ICSBAS),SGS(iState)%nLev,IWALK)
    IP = SG_PHASE(SGS(iState),IWALK)
    return
  end if
  IOFF = IOFF+nBlk
end do

write(u6,*) 'SG_Lucia_Walk_CFTP: CSF index out of range:',ICSF
call Abend()

end subroutine SG_Lucia_Walk_CFTP

!> Check SG_Lucia_Walk for all nConf CSFs of symmetry IREFSM: the split-graph
!> index (SG_NUM) and phase of each walk against the reordering SG_ReOrd
!> (mode 0) of the index vector; for a single RAS partition also the walks
!> of the CDCInfo route against those of the CFTP route.
subroutine SG_Verify_Lucia_Walk(iState,IREFSM,nConf)

use Molcas, only: MxAct
use stdalloc, only: mma_allocate, mma_deallocate

integer(kind=iwp), intent(in) :: iState, IREFSM, nConf
integer(kind=iwp) :: I, IP, IP2, ISG, nBad
integer(kind=iwp) :: IWALK(MxAct), IWALK2(MxAct)
real(kind=wp), allocatable :: Idx(:), Map(:)

call mma_allocate(Idx,nConf,Label='Idx')
call mma_allocate(Map,nConf,Label='Map')
do I=1,nConf
  Idx(I) = real(I,kind=wp)
end do
call SG_ReOrd(iState,IREFSM,0,nConf,Idx,Map)
nBad = 0
do I=1,nConf
  call SG_Lucia_Walk(iState,IREFSM,I,IWALK,IP)
  ISG = SG_NUM(SGS(iState),EXS(iState),IWALK)
  if ((ISG < 1) .or. (ISG > nConf)) then
    nBad = nBad+1
  else if (nint(Map(ISG)) /= merge(-I,I,IP < 0)) then
    nBad = nBad+1
  end if
  if (SGS(iState)%nRsPrt == 1) then
    call SG_Lucia_Walk_CFTP(iState,IREFSM,I,IWALK2,IP2)
    if (any(IWALK(1:SGS(iState)%nLev) /= IWALK2(1:SGS(iState)%nLev)) .or. (IP /= IP2)) nBad = nBad+1
  end if
end do
call mma_deallocate(Map)
call mma_deallocate(Idx)
if (nBad /= 0) then
  write(u6,*) 'SG_Verify_Lucia_Walk: walks of the Lucia CSFs differ from SG_ReOrd, errors:',nBad
  call Abend()
end if
write(u6,'(1X,A,I0,A)') 'SG_Verify_Lucia_Walk: walks of ',nConf,' Lucia CSFs verified against SG_ReOrd'

end subroutine SG_Verify_Lucia_Walk

subroutine GetOccupation_SGUGA(iState,iOpen,mult,iSym,iOC,KCNF)
   use Molcas, only: MxAct
   integer(kind=iwp), intent(in)  :: iState
   integer(kind=iwp), intent(in)  :: iOpen
   integer(kind=iwp), intent(in)  :: mult
   integer(kind=iwp), intent(in)  :: iSym
   integer(kind=iwp), intent(in)  :: iOC
   integer(kind=iwp), intent(out) :: KCNF(MxAct)

   integer(kind=iwp) :: iDOC
   integer(kind=iwp) :: iSOC
   integer(kind=iwp) :: iCL
   integer(kind=iwp) :: i

   KCNF(:) = 0

   iCL = SGS(iState)%CDCInfo(iOpen,mult)%nClosed

   iDOC = SGS(iState)%CDCInfo(iOpen,mult)%OCDOC(iOC,iSym)
   iSOC = SGS(iState)%CDCInfo(iOpen,mult)%OCSOC(iOC,iSym)

   if (iCL > 0) then
      do i=1,iCL
         KCNF(i) = SGS(iState)%CDCInfo(iOpen,mult)%DOC(i,iDOC)
      end do
   end if

   if (iOpen > 0) then
      do i=1,iOpen
         KCNF(iCL+i) = SGS(iState)%CDCInfo(iOpen,mult)%SOC(i,iSOC)
      end do
   end if

end subroutine GetOccupation_SGUGA

subroutine BuildWalk(iState,ICL,IOPEN,IICSF,KCNF,IWALK)
   use Molcas, only: MxAct
   integer(kind=iwp), intent(in)  :: istate
   integer(kind=iwp), intent(in)  :: ICL
   integer(kind=iwp), intent(in)  :: IOPEN
   integer(kind=iwp), intent(in)  :: IICSF
   integer(kind=iwp), intent(in)  :: KCNF(MxAct)
   integer(kind=iwp), intent(out) :: IWALK(MxAct)

   integer(kind=iwp) :: i
   integer(kind=iwp) :: Twos(0:IOPEN)

   if (iOpen == 0) then
      IWALK(:) = 0
      do i=1,ICL
         IWALK(KCNF(i)) = 3
      end do
      return
   endif

   Twos(0:IOPEN) = SGS(iState)%CDCInfo(IOPEN,SGS(iState)%iSpin)%Twos(0:IOPEN,IICSF)

   IWALK(:) = 0

   ! Doubly occupied orbitals
   do i=1,ICL
      IWALK(KCNF(i)) = 3
   end do

   ! Singly occupied orbitals
   do i=1,IOPEN
      if (Twos(i) > Twos(i-1)) then
         IWALK(KCNF(ICL+i)) = 1
      else
         IWALK(KCNF(ICL+i)) = 2
      end if
   end do

End subroutine BuildWalk

!===========================================================
! Build walk from Lucia occupation class and CFTP data
!===========================================================


end module SG_Lucia_Bridge
