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
! Copyright (C) 1989-1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine Seward(ireturn)
!***********************************************************************
! In 1867, William Seward, for 2 cents per acre, purchased             *
! Alaska, a valueless wasteland of ice and snow.                       *
!                                                                      *
! In 1990, Roland Lindh and Ungsik Ryu worked on molecular             *
! integral evaluation, an exhausted scientific area with no            *
! room for innovation.                                                 *
!                                                                      *
! Bowen Liu                                                            *
! April, 1990                                                          *
!***********************************************************************
!***********************************************************************
!                                                                      *
!  Object: Driver for the one and two electron integral program        *
!          SEWARD. SEWARD computes integrals for cartesian and         *
!          spherical harmonic gaussian basis functions.                *
!                                                                      *
!  Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA     *
!          July '89 - May '90                                          *
!                                                                      *
!          Roland Lindh, Dept. of Theoretical Chemistry, University of *
!          Lund, SWEDEN. Modified to use Schwartz inequality for pre-  *
!          screening, July 1991.                                       *
!***********************************************************************

use Period, only: AdCell
use GeoList, only: Centr, Chrg, Mass
use MpmC, only: Coor_MPM
use Basis_Info, only: Basis_Info_Dmp, Basis_Info_Free, Basis_Info_Get, Basis_Info_Init, nBas
use Center_Info, only: Center_Info_Dmp, Center_Info_Free, Center_Info_Get, Center_Info_Init
use Symmetry_Info, only: nIrrep, lIrrep
use LundIO, only: Buf, iDisk, lBuf, Lu_28
use DKH_Info, only: DKroll
use OneDat, only: sNew
use Gateway_Info, only: NEMO, Do_GuessOrb, Do_FckInt, lRP_Post, PkAcc
use RICD_Info, only: Do_RI, Cholesky, DiagCheck, LocalDF
#ifdef _FDE_
use Embedding_Global, only: embPot, embPotInBasis
#endif
use Gateway_global, only: Fake_ERIs, G_Mode, GS_Mode, iPack, iWROpt, Onenly, Primitive_Pass, PrPrt, Run_Mode, S_Mode, Test
use stdalloc, only: mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "warnings.h"
#include "print.fh"
integer(kind=iwp) :: i, iOpt, iRC, iRout, iWrOpt_Save, Lu_One, LuSpool, MaxDax, nChoV(8), nDiff, nDNA
real(kind=wp) :: ChFracMem, DiagErr(4), Dummy(2), TCpu1, TCpu2, TWall1, Twall2
logical(kind=iwp) :: PrPrt_Save, Exists, DoRys, lOPTO, IsBorn, Do_OneEl
!-SVC: identify runfile with a fingerprint
character(len=256) :: cDNA
integer(kind=iwp), external :: ip_of_Work, isFreeUnit
logical(kind=iwp), external :: Reduce_Prt
external :: Integral_WrOut, Integral_WrOut2, Integral_RI_3
interface
  subroutine get_genome(cDNA,nDNA) bind(C,name='get_genome_')
    use, intrinsic :: iso_c_binding, only: c_char
    use Definitions, only: MOLCAS_C_INT
    character(kind=c_char) :: cDNA(*)
    integer(kind=MOLCAS_C_INT) :: nDNA
  end subroutine get_genome
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
!call Seward_Banner()
lOPTO = .false.
call CWTime(TCpu1,TWall1)

! Prologue

iRout = 1
PrPrt_Save = .false. ! dummy initialize
!                                                                      *
!***********************************************************************
!                                                                      *
! Figure out the run_mode
!
! Seward can be run in two different modes!
! GS_Mode: does the work of both Gateway and Seward
! S_Mode:  only the work of Seward
!
! Check if the run file is there

call f_Inquire('RUNFILE',Exists)
if (Exists) then
  call Qpg_iScalar('Run_Mode',Exists)
  if (Exists) then

    ! The Run_mode of the runfile is either GS_Mode or G_Mode

    call Get_iScalar('Run_Mode',Run_Mode)

    ! If the Run_mode is that Gateway is in action then Seward
    ! should be run in S_mode.

    if (Run_Mode == G_Mode) Run_Mode = S_Mode
  else
    Run_Mode = GS_Mode
  end if
else

  ! Seward runs without Gateway

  Run_Mode = GS_Mode
  call MkRun(iRC,0)
  call Put_iScalar('Run_Mode',Run_Mode)

  ! Determine and save the fingerprint of the runfile in a field with
  ! label 'BirthCertificate' if it is empty.  This allows us to
  ! uniquely identify the runfile and any later associated files.

  call qpg_cArray('BirthCertificate',IsBorn,nDNA)
  if (.not. IsBorn) then
    call Get_Genome(cDNA,nDNA)
    call Put_cArray('BirthCertificate',cDNA,nDNA)
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the memory size available

call SetMem('Clear=Off')
!                                                                      *
!***********************************************************************
!                                                                      *
! If Seward is run in S_mode most of the input is already on the
! runfile. If Seward is run in GS_Mode it will handle the input and
! runfile in the conventional way.

call Seward_Init()
if (Run_Mode == S_Mode) then

  ! S_Mode

  DoRys = .true.
  nDiff = 0
  call GetInf(DoRys,nDiff)
  Primitive_Pass = .true.
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! GS_Mode
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call Funi_Init()
  call Basis_Info_Init()
  call Center_Info_Init()

end if ! Run_Mode == S_Mode
!                                                                      *
!***********************************************************************
!                                                                      *
! columbus support: initialize additional items in Runfile
! default: no mixed operation
call Put_iScalar('Columbus',0)
call Put_iScalar('colgradmode',0)
dummy(1) = Zero
dummy(2) = Zero
call Put_dArray('MR-CISD energy',dummy,2)
call NQGrid_Init()
!                                                                      *
!***********************************************************************
!                                                                      *
! Spool the input

LuSpool = 21
call SpoolInp(LuSpool)
! Read the input from input file

call RdCtl_Seward(LuSpool,lOPTO,Do_OneEl)
if (Run_Mode /= S_Mode) then
  call Basis_Info_Dmp()
  call Basis_Info_Free()
  call Basis_Info_Get()
  call Center_Info_Dmp()
  call Center_Info_Free()
  call Center_Info_Get()
end if

call Close_LuSpool(LuSpool)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
do

  ! Process the input.

  call Input_Seward(lOPTO)

  if (Primitive_Pass) then
    PrPrt_Save = PrPrt
    PrPrt = .false.
  else
    PrPrt = PrPrt_Save
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Compute the Nuclear potential energy

  if (.not. Primitive_Pass) call DrvN0()
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (Show) then

    ! Print out basis set information

    write(u6,*)
    write(u6,'(6X,A)') 'Basis set specifications :'
    write(u6,'(6X,A,T30,8(2X,A))') 'Symmetry species',(lIrrep(i),i=0,nIrrep-1)
    write(u6,'(6X,A,T30,8I5)') 'Basis functions',(nBas(i),i=0,nIrrep-1)
    write(u6,*)

  end if

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! If only test case then clean up!

  if (Test) exit
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Write/update information on the run file.

  if (.not. Primitive_Pass) then
    call DmpInf()
    call basis2run()
  end if
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! ONE-ELECTRON INTEGRAL SECTION
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
# ifdef _FDE_
  ! Embedding
  if (embPot .and. (.not. embPotInBasis)) then
    call embPotInit(.false.)
  end if
# endif

  Lu_One = 2
  iOpt = ibset(0,sNew)
  iRC = -1

  ! Generate primimitive integrals only if needed.

  if (Primitive_Pass .and. (DKroll .or. Nemo)) then
    call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
    call OneBas('PRIM')
  else
    call OpnOne(iRC,iOpt,'ONEINT',Lu_One)
  end if
  if (iRC /= 0) then
    call WarningMessage(2,' *** Error in subroutine INPUT ***;     Abend in subroutine OpnOne')
    call Abend()
  end if

  if (Do_OneEl .and. ((.not. Primitive_Pass) .or. DKroll .or. NEMO)) call Drv1El()

  iOpt = 0
  iRC = -1
  call ClsOne(iRC,iOpt)
  if (iRC /= 0) then
    call WarningMessage(2,' *** Error in SEWARD main ***;  Abend in subroutine ClsOne')
    call Abend()
  end if

# ifdef _FDE_
  ! Embedding
  if (embPot .and. (.not. embPotInBasis)) call embPotFreeMem()
# endif
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! If a pass in which primitive integrals where computed do a second
  ! pass.

  if (.not. Primitive_Pass) exit
  Primitive_Pass = .false.
  call Free_iSD()
end do

if (.not. Test) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Branch out if only one-electron integrals are to be computed!

  if (.not. Onenly) then

    ! If ERIs/CD/RI already available, one may want not to redo it!

    if (Fake_ERIs) then
      call set_fake_ERIs()
    else
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
      ! TWO-ELECTRON INTEGRAL SECTION
      !                                                                *
      !*****************************************************************
      !*****************************************************************
      !                                                                *
      if (iWRopt == 0) then

        !----- Molcas format

        if (Cholesky) then ! Cholesky decomposition
          call Cho_MCA_Drv()
          call Get_iArray('NumCho',nChoV,nIrrep)
          write(u6,'(6X,A,T30,8I5)') 'Cholesky vectors',(nChoV(i),i=1,nIrrep)
          write(u6,*)
          write(u6,*)
        else if (Do_RI) then
          if (LocalDF) then
            call Drv2El_LocalDF()
          else
            if (nPrint(iRout) >= 6) then
              write(u6,*)
              write(u6,'(A)') 'Seward processing 2-center and 3-center ERIs'
              write(u6,*)
            end if
            call Drv2El_3Center_RI(Integral_RI_3,Zero)
            call Get_iArray('NumCho',nChoV,nIrrep)
            if (nPrint(iRout) >= 6) then
              write(u6,'(6X,A,T30,8I5)') 'RI vectors',(nChoV(i),i=1,nIrrep)
              write(u6,*)
              write(u6,*)
            end if
          end if
        else
          iWrOpt_Save = iWrOpt
          iWrOpt = 0
          call Sort0()

          call Drv2El(Integral_WrOut2,Zero)

          call Sort1B()
          call Sort2()
          call Sort3(MaxDax)

          if ((.not. Reduce_Prt()) .and. (nPrint(iRout) >= 6)) then
            write(u6,*)
            write(u6,'(A)') ' Integrals are written in MOLCAS2 format'
            if (iPack /= 0) then
              write(u6,'(A)') ' No packing of integrals has been applied'
            else
              write(u6,'(A,G11.4)') ' Packing accuracy =',PkAcc
              write(u6,'(A,I10)') ' Highest disk address written',MaxDax
            end if
            write(u6,'(A)') ' Diagonal and subdiagonal, symmetry allowed 2-el integral blocks are stored on Disk'
          end if
          iWrOpt = iWrOpt_Save
        end if
        !                                                              *
        !***************************************************************
        !                                                              *
      else if (iWRopt == 1) then

        !----- Molecule format (Molcas 1.0)

        Lu_28 = 28
        Lu_28 = isfreeunit(Lu_28)
        call DaName_MF(Lu_28,'BASINT')
        iDisk = 0
        lBuf = ip_of_Work(Buf%r_End)-ip_of_Work(Buf%Buf(1))

        call Drv2El(Integral_WrOut,Zero)

        call dDafile(Lu_28,1,Buf%Buf,lBuf,iDisk)
        Buf%nUt = -1
        call dDafile(Lu_28,1,Buf%Buf,lBuf,iDisk)
        write(u6,*)
        write(u6,'(A)') ' Integrals are written in MOLCAS1 format'
        !write(u6,'(I10,A)') IntTot,' Integrals written on Disk'
        !                                                              *
        !***************************************************************
        !                                                              *
      else

        call WarningMessage(2,'Seward: Invalid value of iWRopt!')
        call Abend()

      end if

    end if ! Fake_ERIs
  end if   ! Onenly
end if     ! Test
!                                                                      *
!***********************************************************************
!                                                                      *
! At the end of the calculation free all memory to check for
! corruption of the memory.

call ClsSew()
if (allocated(AdCell)) call mma_deallocate(AdCell)
call mma_deallocate(Coor_MPM)
call mma_deallocate(Chrg)
call mma_deallocate(Mass)
call mma_deallocate(Centr)
!                                                                      *
!***********************************************************************
!                                                                      *
call CWTime(TCpu2,TWall2)
!                                                                      *
!***********************************************************************
!                                                                      *
! Diagonal ERI check

if (Cholesky .or. Do_RI) then
  if (DiagCheck) then
    write(u6,*) ' ==== Start Diagonal ERI check  ===='
    call Cho_X_init(irc,ChFracMem)
    if (irc /= 0) then
      call WarningMessage(2,' Seward: Non-zero rc in Cho_X_init.')
      call Abend()
    end if
    call Cho_X_CheckDiag(irc,DiagErr)
    if (irc /= 0) then
      call WarningMessage(2,' Seward: Non-zero rc in Cho_X_CheckDiag.')
      call Abend()
    end if
    call Cho_X_Final(irc)
    if (irc /= 0) then
      call WarningMessage(2,' Seward: Non-zero rc in Cho_X_Final.')
      call Abend()
    end if
    write(u6,*)
    write(u6,*) ' ====  End  Diagonal ERI check  ===='
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Automatic run of GuessOrb

if (Do_GuessOrb .and. Do_FckInt) call GuessOrb(iReturn,.false.)
if ((.not. Prprt) .and. Do_OneEl) call Put_NucAttr()
!                                                                      *
!***********************************************************************
!                                                                      *
! Epilogue

if (nPrint(iRout) >= 6) call FastIO('STATUS')

if (Test) then
  ireturn = _RC_EXIT_EXPECTED_
else if (lRP_Post) then
  ireturn = _RC_INVOKED_OTHER_MODULE_
else
  ireturn = _RC_ALL_IS_WELL_
end if

return

end subroutine Seward
