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
! Copyright (C) 1998, Roland Lindh                                     *
!               2004,2012, Thomas Bondo Pedersen                       *
!               2007, Francesco Aquilante                              *
!***********************************************************************

subroutine MP2_Driver(ireturn)
!***********************************************************************
!     New MP2 driver for Molcas                                        *
!                                                                      *
!    Author: R. Lindh                                                  *
!            Dept. of Chemical Physics                                 *
!            University of Lund, Sweden                                *
!            November 7, 1998                                          *
!                                                                      *
!    Modified:                                                         *
!                                                                      *
!       - code using Cholesky vectors directly                         *
!         October 2004, T. B. Pedersen                                 *
!         Dept. of Theoretical Chemistry                               *
!         University of Lund, Sweden                                   *
!                                                                      *
!       - code for the "Scaled Opposite-Spin" (SOS) MP2                *
!         May 2007, F. Aquilante                                       *
!         Dept. of Theoretical Chemistry                               *
!         University of Lund, Sweden                                   *
!                                                                      *
!       - code for Laplace-SOS-MP2 for Cholesky/DF and LDF             *
!         November-December 2012, T. B. Pedersen                       *
!         Centre for Theoretical and Computational Chemistry           *
!         Dept. of Chemistry                                           *
!         University of Oslo, Norway                                   *
!***********************************************************************

use MBPT2_Global, only: CMO, DoCholesky, DoDF, DoLDF, EOcc, EOrb, EVir, FnIntA, FnIntM, iPL, LuHLF1, LuHLF2, LuHLF3, LuIntA, &
                        LuIntM, MBPT2_Clean, NamAct, nBas
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: i, iOpt, iPrc, irc, iSym, iTol, iTst, iType, lthCMO, lthEOr, nAsh(8), nDel_tra(8), nFro_tra(8), nIsh(8), nOccT
real(kind=wp) :: E0, E2BJAI, ESCF, ESSMP2, Etot, REFC, Shanks1_E, t1dg, t1nrm, TCPE(4), TCPT, TIOE(4), TIOT
logical(kind=iwp) :: Conventional, IsDirect, Exists, Ready
character(len=8) :: Method, Method1
real(kind=wp), allocatable :: T1amp(:)
logical(kind=iwp), parameter :: Debug = .false.
integer(kind=iwp), external :: Cho_X_GetTol, ip_of_Work
real(kind=wp), external :: ddot_, Seconds
#include "Molcas.fh"
#include "trafo.fh"
#include "corbinf.fh"
#include "chomp2_cfg.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
TCPT = seconds()

!***********************************************************************
! Check so it is a RHF-SCF reference that is being used.
! TBP, November 2012: do not quit, just issue a warning!
call Get_cArray('Relax Method',Method1,8)
if ((Method1(1:7) /= 'RHF-SCF') .and. (Method1(1:5) /= 'MBPT2')) then
  write(u6,*)
  call WarningMessage(1,'MP2 implementation intended for RHF references only')
  write(u6,'(A,A)') 'MBPT2 WARNING: Reference function according to RunFile:',Method1
  write(u6,'(A)') 'I''ll assume you know what you''re doing and continue'
  write(u6,*)
  call xFlush(u6)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Figure out if it's a cholesky run, DF run, LDF run
DoCholesky = .false.
DoDF = .false.
DoLDF = .false.
call DecideOnCholesky(DoCholesky)
call DecideOnDF(DoDF)
call DecideOnLocalDF(DoLDF)
!                                                                      *
!***********************************************************************
!                                                                      *
! read COMFILE Interface file from SCF and allocate memory...

call TIMING(TCPE(1),TCPT,TIOE(1),TIOT)
call RdMBPT()
lthCMO = size(CMO)
lthEOr = size(EOrb)

call mma_allocate(EOcc,lthEOr,label='EOcc')
call mma_allocate(EVir,lthEOr,label='EVir')
EOcc(:) = Zero
EVir(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! read Program Input Cards...

call RdInp(CMO,EOrb,EOcc,EVir,iTst,ESCF)
!                                                                      *
!***********************************************************************
!                                                                      *
call DelGHOST_MBPT()
!                                                                      *
!***********************************************************************
!                                                                      *
E2BJAI = Zero
REFC = Zero

Wref = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out input parameters

call PrInp_MBPT2(EOcc,EVir,iTst)
!                                                                      *
!***********************************************************************
!                                                                      *
! Copy pointers to orbital energies to chomp2_dec.fh
! Needed for amplitude Cholesky decomposition.

if (DoCholesky) then
  call ChoMP2_SetPtsOen(EOcc,EVir)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Check, if there is an ORDINT file...

call f_Inquire(FNINTA,Exists)
call DecideOnDirect(.true.,Exists,IsDirect,DoCholesky)
call TIMING(TCPE(2),TCPT,TIOE(2),TIOT)
call TIMING(TCPE(3),TCPT,TIOE(3),TIOT)

if (DoT1amp) then
  l_T1 = nOcc(1)*nExt(1)
  nOccT = nOcc(1)
  do iSym=2,nSym
    l_T1 = l_T1+nOcc(iSym)*nExt(iSym)
    nOccT = nOccT+nOcc(iSym)
  end do
  call mma_allocate(T1amp,l_T1,label='T1amp')
  call Thouless_T1(CMO,nSym,nBas,nFro,nOcc,nExt,T1amp)
  t1nrm = ddot_(l_T1,T1amp,1,T1amp,1)
  t1dg = sqrt(t1nrm/nOccT)
  write(u6,'(A,F8.4)') '       T1 diagnostic : ',t1dg
  write(u6,*)
  iOffT1(1) = ip_of_Work(T1amp(1))-1
  do i=2,nSym
    iOffT1(i) = iOffT1(i-1)+nOcc(i-1)*nExt(i-1)
  end do
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
if (DoLDF) then ! LDF
  Conventional = .false.
  Ready = .false.
  if (Laplace .and. SOS_MP2) then
    call WarningMessage(2,'LDF-Laplace-SOS-MP2 not implemented yet!')
    call SysHalt('mp2_driver')
  else
    call WarningMessage(2,'Only LDF-Laplace-SOS-MP2 implemented!')
    call SysHalt('mp2_driver')
  end if
else if (DoCholesky .and. (ChoAlg > 0) .and. (.not. SOS_mp2) .and. (.not. FNOMP2) .and. (.not. LovMP2)) then
  Conventional = .false.
  Ready = .false.
  call ChoMP2_Drv(irc,E2BJAI,CMO,EOcc,EVir)
  if (irc /= 0) then
    write(u6,*) 'MP2_Driver: ChoMP2_Drv returned ',irc
    call SysAbendMsg('MP2_Driver','Non-zero return code from ChoMP2_Drv',' ')
  else
    Ready = .true.
  end if
else if (DoCholesky .and. SOS_mp2) then ! CD/DF SOS-MP2
  Conventional = .false.
  Ready = .false.
  if (Laplace) then
    call ChoMP2_Drv(irc,E2BJAI,CMO,EOcc,EVir)
    if (irc /= 0) then
      write(u6,*) 'MP2_Driver: ChoMP2_Drv returned ',irc
      call SysAbendMsg('MP2_Driver','Non-zero return code from ChoMP2_Drv',' ')
    else
      Ready = .true.
    end if
  else
    call Cho_SOSmp2_Drv(irc,E2BJAI,CMO,EOcc,EVir)
    if (irc /= 0) then
      write(u6,*) 'SOS-MP2_Driver: Cho_SOSmp2_Drv returned ',irc
      call SysAbendMsg('SOS-MP2_Driver','Non-zero return code from Cho_SOSmp2_Drv',' ')
    else
      Ready = .true.
    end if
  end if
! CD/DF Frozen Natural Orbital MP2
else if (DoCholesky .and. FNOMP2) then
  Conventional = .false.
  Ready = .false.
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A)') ' Start FNO-MP2 section '
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A,8I4)')
  write(u6,'(A,I3,A)') ' NOs specified: ',int(vkept*100),'% of the total virtual space'
  call FNOMP2_Drv(irc,E2BJAI,CMO,EOcc,EVir)
  if (irc /= 0) then
    write(u6,*) 'MP2 driver: FNOMP2_Drv returned ',irc
    call SysAbendMsg('MP2 driver','Non-zero return code from FNOMP2_Drv',' ')
  else
    Ready = .true.
  end if
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A)') ' End FNO-MP2 section '
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A,8I4)')
  write(u6,'(A,8I4)')
else if (DoCholesky .and. LovMP2) then ! CD/DF Localized O-V MP2
  Conventional = .false.
  Ready = .false.
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A)') ' Start LovMP2 section '
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A,8I4)')
  call LovMP2_Drv(irc,E2BJAI,CMO,EOcc,EVir,NamAct,nActa,ThrLov,DoMP2,all_Vir)
  if (irc /= 0) then
    write(u6,*) 'MP2 driver: LovMP2_Drv returned ',irc
    call SysAbendMsg('MP2 driver','Non-zero return code from LovMP2_Drv',' ')
  else
    Ready = .true.
  end if
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A)') ' End LovMP2 section '
  write(u6,'(A)') '-------------------------------------------------------'
  write(u6,'(A,8I4)')
  write(u6,'(A,8I4)')
else ! conventional (possibly with Cholesky)
  Conventional = .true.
  if (DoCholesky) then
    if ((ChoAlg == 0) .and. (iPL >= 2)) then
      write(u6,*) 'Conventional algorithm used.'
      write(u6,*)
      write(u6,*) 'Integrals generated from Cholesky vectors (Algorithm 0):'
    else
      call SysHalt('mp2_driver')
    end if
  else
    if (iPL >= 2) write(u6,*) 'Conventional algorithm used...'
  end if
  if (iTst /= 0) then
    call finalize()
    return
  end if

  ! Use the driver from the CASPT2 code.

  call DaName_MF_wa(LuIntM,FnIntM)

  if (DoDens) then

    ! Same as just above but keeping frozen and deleted orbitals in the
    ! transformation
    do i=1,8
      nIsh(i) = nOrb(i)+nDel(i)
    end do
    nFro_tra(:) = 0
    nDel_tra(:) = 0
  else
    nIsh(:) = nOcc
  end if

  nAsh(:) = 0
  if (.not. DoDens) then
    ! PAM Jan 2013: Set correct nOrb:
    do i=1,8
      norb(i) = norb(i)-nfro(i)
    end do
    call SetUp_CASPT2_Tra(nSym,nBas,nOrb,nIsh,nAsh,nFro,nDel,CMO,lthCMO,LuIntM,LuHlf1,LuHlf2,LuHlf3)
    ! End of patch
  else
    call SetUp_CASPT2_Tra(nSym,nBas,nIsh,nIsh,nAsh,nFro_tra,nDel_tra,CMO,lthCMO,LuIntM,LuHlf1,LuHlf2,LuHlf3)
  end if
  if (.not. DoCholesky) then
    iRC = -1
    iOpt = 0
    call OpnOrd(iRC,iOpt,FnIntA,LuIntA)
    if (iRC /= 0) then
      write(u6,*) 'mp2_driver: error opening MOLINT'
      call Abend()
    end if
  end if

  iType = 1 ! Means that TraCtl_Drv is called by MP2
  call TraCtl_Drv(iType,.true.,8)

  call DaClos(LuHlf1)
  call DaClos(LuHlf2)
  call DaClos(LuHlf3)

  call TIMING(TCPE(3),TCPT,TIOE(3),TIOT)
  ! PRINT TRANSFORMED INTEGRALS (USED ONLY FOR DEBUGGING PURPOSES)
  if (Debug) then
    iPrc = 0
    call RDInt2_MP2(iPrc)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (DoDens) then

    ! Obtain variational density if requested by user

    call Mp2Dens_drv(E2BJAI,REFC)
  else

    ! COMPUTE CORRELATION ENERGY CONTRIBUTION

    call BJAI(IADOUT,EOcc,EVir,E2BJAI,REFC)
  end if
  call TIMING(TCPE(4),TCPT,TIOE(4),TIOT)

  Ready = .true.

  ! Close files (for Cholesky, OrdInt was never opened)
  if (.not. DoCholesky) then
    iRc = -1
    call ClsOrd(iRc)
    if (iRc /= 0) then
      write(u6,*) 'MP2_Driver: Error closing ORDINT'
      call Abend()
    end if
  end if
  call DaClos(LuIntM)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out resulting energy, etc.

if (Ready) then

  Wref = One/(One+Wref)   ! Note: this is Cref**2

  if (DoLDF) then
    call WarningMessage(2,'LDF should not be implemented....')
    call SysHalt('mp2_driver')
  else if (DoCholesky .and. (ChoAlg > 0) .and. (.not. SOS_mp2) .and. (.not. LovMP2) .and. (.not. FNOMP2)) then
    if (iPL >= 2) write(u6,'(3(/6X,A,F20.10,A)//6X,A,F20.10,A//6X,A,F15.5)') &
      ' SCF energy                           =',ESCF,' a.u.', &
      ' Second-order correlation energy      =',E2BJAI,' a.u.', &
      ' ( Opposite-Spin contribution         =',-EOSMP2,' )', &
      ' Total energy                         =',E2BJAI+ESCF,' a.u.', &
      ' Reference weight ( Cref**2 )         =',Wref
  else if (DoCholesky .and. LovMP2) then
    ESSMP2 = E2BJAI+EOSMP2
    if (iPL >= 2) write(u6,'(4(/6X,A,F20.10,A)//6X,A,F20.10,A//6X,A,F15.5)') &
      ' SCF energy                           =',ESCF,' a.u.', &
      ' Second-order correlation energy      =',E2BJAI,' a.u.', &
      ' ( Opposite-Spin contribution         =',-EOSMP2,' )', &
      ' ( Same-Spin contribution             =',ESSMP2,' )', &
      ' Total energy                         =',E2BJAI+ESCF,' a.u.', &
      ' Reference weight ( Cref**2 )         =',Wref
  else if (DoCholesky .and. FNOMP2) then
    ESSMP2 = E2BJAI+EOSMP2-XEMP2
    if (iPL >= 2) write(u6,'(5(/6X,A,F20.10,A)//6X,A,F20.10,A//6X,A,F15.5)') &
      ' SCF energy                           =',ESCF,' a.u.', &
      ' Second-order correlation energy      =',E2BJAI,' a.u.', &
      ' ( Opposite-Spin contribution         =',-EOSMP2,' )', &
      ' ( Same-Spin contribution             =',ESSMP2,' )', &
      ' ( Truncation error estimate          =',XEMP2,' )', &
      ' Total energy                         =',E2BJAI+ESCF,' a.u.', &
      ' Reference weight ( Cref**2 )         =',Wref
  else if (DoCholesky .and. SOS_mp2) then ! CD/DF SOS-MP2 energy
    E2BJAI = C_os*E2BJAI
    if (iPL >= 2) then
      if (Laplace) then
        write(u6,'(/,6X,A,I4)') ' Number of Laplace grid points:',Laplace_nGridPoints
        write(u6,'(3(/6X,A,F20.10,A)//6X,A,F20.10,A)') ' Opposite-Spin (OS) scaling factor    =',C_os,'     ', &
                                                       ' SCF energy                           =',ESCF,' a.u.', &
                                                       ' L-SOS 2nd-order correlation energy   =',E2BJAI,' a.u.', &
                                                       ' Total L-SOS-MP2 energy               =',E2BJAI+ESCF,' a.u.'
      else
        write(u6,'(3(/6X,A,F20.10,A)//6X,A,F20.10,A)') ' Opposite-Spin (OS) scaling factor    =',C_os,'     ', &
                                                       ' SCF energy                           =',ESCF,' a.u.', &
                                                       ' SOS 2nd-order correlation energy     =',E2BJAI,' a.u.', &
                                                       ' Total SOS-MP2 energy                 =',E2BJAI+ESCF,' a.u.'
      end if
    end if
  else
    WRef = REFC**2
    if (iPL >= 2) write(u6,'(2(/6X,A,F20.10,A)//6X,A,F20.10,A/6X,A,F15.5)') &
      ' SCF energy                           =',ESCF,' a.u.', &
      ' Second-order correlation energy      =',E2BJAI,' a.u.', &
      ' Total energy                         =',E2BJAI+ESCF,' a.u.', &
      ' Reference weight ( Cref**2 )         =',Wref
  end if
  if (iPL >= 2) then
    ETot = E2BJAI+ESCF
    write(u6,*)
    call PrintResult(u6,'(6X,A,T50,F19.10)','Total MBPT2 energy',0,'',[ETot],1)
  end if
  call xFlush(u6)
  if (iPL >= 2) write(u6,*)
  call Store_Energies(1,E2BJAI+ESCF,1)
  Method = 'MBPT2   '
  call Put_cArray('Relax Method',Method,8)
  if (DoDens) call Prpt()
  if (iPL >= 2) write(u6,*)
else
  write(u6,*) ' Energy evaluation not completed.'
end if

! Shanks-type series convergence acceleration

call Compute_Shanks(ESCF,E2BJAI+ESCF,EOrb,lthEOr,nBas,nFro,nOcc,nSym,E0,Shanks1_E)

if (iPL >= 2) then
  write(u6,'(6X,A,A,(F20.10),A)') ' Zeroth-order energy (E0)   ','          =',E0,' a.u.'
  write(u6,*)

  write(u6,'(6X,A,A,(F20.10),A)') ' Shanks-type energy S1(E)   ','          =',Shanks1_E,' a.u.'
  write(u6,*)
  write(u6,*)
end if

! PRINT PROCESSING AND TIMING INFORMATION

if (Conventional .and. (iPL >= 2)) then
  write(u6,'(//6X,A//6X,A/6X,A/)') ' Data processing and timing information:', &
                                   ' Section                                              time(sec)', &
                                   '                                                    CPU  Elapsed'
  write(u6,'(6X,A,2(2X,F8.2))') 'Input data processing                        ',TCPE(2)-TCPE(1),TIOE(2)-TIOE(1)
  write(u6,'(6X,A,2(2X,F8.2))') 'Transformation of integrals                  ',TCPE(3)-TCPE(2),TIOE(3)-TIOE(2)
  write(u6,'(6X,A,2(2X,F8.2))') 'MBPT2 calculations (BJAI)                    ',TCPE(4)-TCPE(3),TIOE(4)-TIOE(3)
  write(u6,'(6X,A,2(2X,F8.2))') 'Total MBPT2 calculations                     ',TCPE(4)-TCPE(1),TIOE(4)-TIOE(1)
  write(u6,*)
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
iTol = Cho_X_GetTol(8)
call Add_Info('E_MP2',[E2BJAI+ESCF],1,iTol)
call Add_Info('HF ref weight',[WREF],1,iTol)
!                                                                      *
!***********************************************************************
!                                                                      *

! Close up calculation
call finalize()

return

contains

subroutine finalize()
  call MBPT2_Clean()
  if (DoT1amp) call mma_deallocate(T1amp)
  ireturn = 0
end subroutine finalize

end subroutine MP2_Driver
