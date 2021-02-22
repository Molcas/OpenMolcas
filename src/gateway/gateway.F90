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

use Period
use GeoList
use MpmC
use Basis_Info
use Center_Info
use external_centers, only: iXPolType, XF
use Temporary_parameters, only: Primitive_Pass, Expert, VarR, VarT, DirInt
use Sizes_of_Seward, only: S
use RICD_Info, only: Do_RI, Cholesky, Cho_OneCenter
use Symmetry_Info, only: nIrrep
use Definitions, only: iwp

implicit real*8(A-H,O-Z)
integer AixRm
external AixRm
#include "Molcas.fh"
#include "status.fh"
#include "gateway.fh"
#include "rctfld.fh"
#include "stdalloc.fh"
#include "print.fh"
character xLblCnt(MxAtom)*(LENIN)
parameter(nMamn=MaxBfn+MaxBfn_Aux)
character*(LENIN8), allocatable :: Mamn(:)
logical lOPTO, Pseudo, Do_OneEl
logical Cho_1Center
!VV      LOGICAL GA_USES_MA,GA_MEMORY_LIMITED
!-SVC: identify runfile with a fingerprint
character cDNA*256
logical IsBorn, Found
real*8, allocatable :: DCo(:,:), DCh(:), DCh_Eff(:)
integer, allocatable :: nStab(:)
!                                                                      *
!***********************************************************************
!                                                                      *
! Call Gateway_banner()
iReturn = 0

! If Gateway is running the Run_Mode on the runfile should always
! be G_Mode.

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
!     Spool the input

LuSpool = 21
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
! them back from the runfile now allocating them to the proper
! size.

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
call DmpInf()
!                                                                      *
!***********************************************************************
!                                                                      *
! Produce minimal set of entries on the runfile to facilitate
! Grid_It's and ExpBas's needs.

call Drvn0()

call Put_cArray('Unique Basis Names',Mamn(1),(LENIN8)*S%nDim)
call Put_iArray('NBAS',nBas,nIrrep)
call basis2run()
call mma_deallocate(Mamn)

! Generate list of unique atoms

nNuc = 0
do iCnttp=1,nCnttp
  if (.not. dbsc(iCnttp)%pChrg .and. .not. dbsc(iCnttp)%Frag .and. .not. dbsc(iCnttp)%Aux) nNuc = nNuc+dbsc(iCnttp)%nCntr
end do

call mma_allocate(DCo,3,nNuc)
call mma_allocate(DCh,nNuc)
call mma_allocate(DCh_Eff,nNuc)
call mma_allocate(nStab,nNuc)
mdc = 0
iNuc = 0
do iCnttp=1,nCnttp
  if (.not. dbsc(iCnttp)%pChrg .and. .not. dbsc(iCnttp)%Frag .and. .not. dbsc(iCnttp)%Aux) then
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdc = mdc+1
      iNuc = iNuc+1
      DCo(1:3,iNuc) = dbsc(iCnttp)%Coor(1:3,iCnt)
      DCh_Eff(iNuc) = dbsc(iCnttp)%Charge
      DCh(iNuc) = dble(dbsc(iCnttp)%AtmNr)
      xLblCnt(iNuc) = dc(mdc)%LblCnt(1:LENIN)
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
call Put_cArray('Unique Atom Names',xLblCnt(1),LENIN*nNuc)
call Put_iArray('nStab',nStab,nNuc)

call mma_deallocate(nStab)
call mma_deallocate(DCo)
call mma_deallocate(DCh)
call mma_deallocate(DCh_Eff)

! Manipulate the option flag

iOption = 0
if (DirInt) iOption = ior(iOption,1)
if (Expert) iOption = ior(iOption,2)
if (lRF) iOption = ior(iOption,4)
if (lLangevin .or. iXPolType > 0) iOption = ior(iOption,8)
if (PCM) then
  iOption = ior(iOption,16)
  nPCM_Info = 0
  call Put_iScalar('PCM info length',nPCM_Info)
end if
iOption = ior(iOption,32)
if (lRF .and. .not. PCM) iOption = ior(iOption,2**7)
Pseudo = .false.
do iCnttp=1,nCnttp
  Pseudo = Pseudo .or. (dbsc(iCnttp)%pChrg .and. dbsc(iCnttp)%Fixed)
end do
if (allocated(XF) .or. Pseudo) then
  iOption = ior(iOption,2**7)
  iOption = ior(iOption,2**8)
end if
if (VarT) iOption = ior(iOption,2**7)
if (VarR) iOption = ior(iOption,2**8)
! 2el-integrals from the Cholesky vectors
if (Cholesky .or. Do_RI) iOption = ior(iOption,2**9)
! RI-Option
if (Do_RI) iOption = ior(iOption,2**10)
! 1C-CD
Cho_1Center = Get_Cho_1Center()
if (Cholesky .and. Cho_1Center) iOption = ior(iOption,2**12)
Cho_OneCenter = Cho_1Center
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

contains

function Get_Cho_1Center() result(res)
  logical(kind=iwp) :: res
# include "cholesky.fh"
  res = Cho_1Center
end function Get_Cho_1Center

end subroutine Gateway
