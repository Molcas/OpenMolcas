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
! Copyright (C) 2006, Roland Lindh                                     *
!***********************************************************************

subroutine Gateway(iReturn)
!***********************************************************************
!                                                                      *
!     In the fall of 2006 Seward, Alaska and McKinley were in need     *
!     of a sibling. The new code would be the gateway into the         *
!     MOLCAS world for our future GUI. After some thoughts we          *
!     decided to call the gateway Gateway after the small village      *
!     Gateway, Alaska.                                                 *
!                                                                      *
!     Roland Lindh, 26th September 2006.                               *
!                                                                      *
!***********************************************************************

use Period, only: AdCell
use GeoList, only: Centr, Chrg, Mass
use MpmC, only: Coor_MPM
use Basis_Info, only: dbsc, nBas, nCnttp, basis_info_dmp, basis_info_init, basis_info_get, basis_info_free
use Center_Info, only: dc, center_info_dmp, center_info_init, center_info_get, center_info_free
use external_centers, only: iXPolType, XF
use Gateway_global, only: DirInt, Expert, G_Mode, Primitive_Pass, Run_Mode
use Sizes_of_Seward, only: S
use RICD_Info, only: Cho_OneCenter, Chol => Cholesky, Do_DCCD, Do_RI
use Cholesky, only: Cho_1Center
use Symmetry_Info, only: nIrrep, VarR, VarT
use rctfld_module, only: lLangevin, lRF, nPCM_Info, PCM
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: iReturn
#include "Molcas.fh"
#include "print.fh"
integer(kind=iwp) :: iCnt, iCnttp, iNuc, iOption, iRc, iter_S, LuSpool, mdc, nDNA, nNuc
integer(kind=iwp), parameter :: nMamn = MaxBfn+MaxBfn_Aux
character(len=LenIn) :: xLblCnt(MxAtom)
logical(kind=iwp) :: lOPTO, Pseudo, Do_OneEl, IsBorn, Found
!-SVC: identify runfile with a fingerprint
character(len=256) :: cDNA
character(len=LenIn8), allocatable :: Mamn(:)
real(kind=wp), allocatable :: DCo(:,:), DCh(:), DCh_Eff(:)
integer(kind=iwp), allocatable :: nStab(:)
integer(kind=iwp), external :: AixRm
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
!call Gateway_banner()
iReturn = 0

! If Gateway is running the Run_Mode on the runfile should always be G_Mode.

Run_Mode = G_Mode
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
!                                                                      *
!***********************************************************************
!                                                                      *
! Get the memory size available

call SetMem('Clear=Off')
!                                                                      *
!***********************************************************************
!                                                                      *
! Initialize common blocks and start from scratch

call Seward_Init()
call Funi_Init()
call NQGrid_Init()
call Basis_Info_Init()
call Center_Info_Init()
!                                                                      *
!***********************************************************************
!                                                                      *
! Spool the input

call SpoolInp(LuSpool)
!                                                                      *
!***********************************************************************
!                                                                      *
! Remove possible leftover files

call f_Inquire('UDC.Gateway',Found)
if (Found) iRC = AixRm('UDC.Gateway')
call f_Inquire('UDC.NG',Found)
if (Found) iRC = AixRm('UDC.NG')
!                                                                      *
!***********************************************************************
!                                                                      *
! Read the input.

lOPTO = .false.
call RdCtl_Seward(LuSpool,lOPTO,Do_OneEl)

! Write the Basis_Info data to file. Release the arrays and read
! them back from the runfile now allocating them to the proper size.

call Basis_Info_Dmp()
call Basis_Info_Free()
call Basis_Info_Get()
call Center_Info_Dmp()
call Center_Info_Free()
call Center_Info_Get()
!                                                                      *
!***********************************************************************
!                                                                      *
! Close Spool file

call Close_LuSpool(LuSpool)
!                                                                      *
!***********************************************************************
!                                                                      *
! Print out section

call Print_Symmetry()
call Flip_Flop(.false.)
call Print_Basis(lOPTO)
call Print_Geometry(0)
call Print_Isotopes()
if (nPrint(2) > 0) nPrint(117) = 6
call RigRot(Centr,Mass,S%kCentr)
call Print_Basis2()
call Print_OpInfo()
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate the SO/AO basis set here.

Primitive_Pass = .false.
call Flip_Flop(Primitive_Pass)
call mma_allocate(Mamn,nMamn,label='Mamn')
call SOCtl_Seward(Mamn,nMamn)
!                                                                      *
!***********************************************************************
!                                                                      *
if (lRF .and. (.not. PCM)) VarT = .true.
Pseudo = .false.
do iCnttp=1,nCnttp
  Pseudo = Pseudo .or. (dbsc(iCnttp)%pChrg .and. dbsc(iCnttp)%Fixed)
end do
if (allocated(XF) .or. Pseudo) then
  VarR = .true.
  VarT = .true.
end if
call DmpInf()
!                                                                      *
!***********************************************************************
!                                                                      *
! Produce minimal set of entries on the runfile to facilitate
! Grid_It's and ExpBas's needs.

call Drvn0()

call Put_cArray('Unique Basis Names',Mamn(1),(LenIn8)*S%nDim)
call Put_iArray('NBAS',nBas,nIrrep)
call basis2run()
call mma_deallocate(Mamn)

! Generate list of unique atoms

nNuc = 0
do iCnttp=1,nCnttp
  if ((.not. dbsc(iCnttp)%pChrg) .and. (.not. dbsc(iCnttp)%Frag) .and. (.not. dbsc(iCnttp)%Aux)) nNuc = nNuc+dbsc(iCnttp)%nCntr
end do

call mma_allocate(DCo,3,nNuc)
call mma_allocate(DCh,nNuc)
call mma_allocate(DCh_Eff,nNuc)
call mma_allocate(nStab,nNuc)
mdc = 0
iNuc = 0
do iCnttp=1,nCnttp
  if ((.not. dbsc(iCnttp)%pChrg) .and. (.not. dbsc(iCnttp)%Frag) .and. (.not. dbsc(iCnttp)%Aux)) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      iNuc = iNuc+1
      DCo(1:3,iNuc) = dbsc(iCnttp)%Coor(1:3,iCnt)
      DCh_Eff(iNuc) = dbsc(iCnttp)%Charge
      DCh(iNuc) = real(dbsc(iCnttp)%AtmNr,kind=wp)
      xLblCnt(iNuc) = dc(mdc)%LblCnt(1:LenIn)
      nStab(iNuc) = dc(mdc)%nStab
    end do
  else
    mdc = mdc+dbsc(iCnttp)%nCntr
  end if
end do
call Put_iScalar('Unique atoms',nNuc)
call Put_dArray('Unique Coordinates',DCo,3*nNuc)
call Put_dArray('Nuclear charge',DCh,nNuc)
call Put_dArray('Effective nuclear Charge',DCh_Eff,nNuc)
call Put_cArray('Unique Atom Names',xLblCnt(1),LenIn*nNuc)
call Put_iArray('nStab',nStab,nNuc)

call mma_deallocate(nStab)
call mma_deallocate(DCo)
call mma_deallocate(DCh)
call mma_deallocate(DCh_Eff)

! Manipulate the option flag

iOption = 0
if (DirInt) iOption = ibset(iOption,0)
if (Expert) iOption = ibset(iOption,1)
if (lRF) iOption = ibset(iOption,2)
if (lLangevin .or. (iXPolType > 0)) iOption = ibset(iOption,3)
if (PCM) then
  iOption = ibset(iOption,4)
  nPCM_Info = 0
  call Put_iScalar('PCM info length',nPCM_Info)
end if
iOption = ibset(iOption,5)
! 2el-integrals from the Cholesky vectors
if (Chol .or. Do_RI) iOption = ibset(iOption,9)
! RI-Option
if (Do_RI) iOption = ibset(iOption,10)
! 1C-CD
if (Chol .and. Cho_1Center) iOption = ibset(iOption,12)
Cho_OneCenter = Cho_1Center
if (Do_DCCD) iOption = ibset(iOption,13)
call Put_iScalar('System BitSwitch',iOption)
iter_S = 0
call Put_iScalar('Saddle Iter',iter_S)
!                                                                      *
!***********************************************************************
!                                                                      *
call ClsSew()
if (allocated(AdCell)) call mma_deallocate(AdCell)
call mma_deallocate(Coor_MPM)
call mma_deallocate(Chrg)
call mma_deallocate(Mass)
call mma_deallocate(Centr)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Gateway
