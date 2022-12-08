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
! Copyright (C) 2010, Mickael G. Delcey                                *
!               2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine processRP(KeepGroup,SymThr)

use External_Centers, only: nRP, RP_Centers
#ifndef _HAVE_EXTRA_
use XYZ, only: Clear_XYZ, Out_Raw, Parse_Group, Read_XYZ, Symmetry
#endif
use Definitions, only: wp, iwp

implicit none
character(len=180), intent(in) :: KeepGroup
real(kind=wp), intent(inout) :: SymThr
integer(kind=iwp) :: i, LuRP
character(len=180) :: KG
#ifdef _HAVE_EXTRA_
character(len=180) :: Key, KWord, minGroup
character(len=180), external :: Get_Ln
#else
real(kind=wp), allocatable :: DumRot(:,:,:), DumTrans(:,:)
#endif
integer(kind=iwp), external :: IsFreeUnit

!***********************************************************************
!                                                                      *
!    A silly routine to try to handle symmetry in RP-Coord section     *
!    verify that R, P and the standard structures have one common      *
!    symmetry. The user should preferably use a defined symmetry if    *
!    the TS symmetry is lower than R and P                             *
!                                                                      *
!    M.G. Delcey     June 2010                                         *
!    Lund University                                                   *
!                                                                      *
!    Adaptation to OpenMolcas                                          *
!    I. Fdez. Galvan   July 2017                                       *
!                                                                      *
!***********************************************************************

! If C1, all is already done

KG = KeepGroup
call UpCase(KG)
if ((KG(1:1) == 'E') .or. (KG(1:2) == 'C1')) then
  KG = 'NOSYM'
end if
if (KG(1:5) == 'NOSYM') return

#ifdef _HAVE_EXTRA_
LuRP = IsFreeUnit(10)
call molcas_open(LuRP,'findsym.out')
KWord = Get_Ln(LuRP)
minGroup = Get_Ln(LuRP)
close(LuRP)
if (index(KWord,'C1') /= 0) return

! Else findsym will find the group

call findsymf('findsym.RP1',KG,SymThr,ierr)
if (ierr /= 0) call Error(1)
call molcas_open(LuRP,'findsym.out')
Key = Get_Ln(LuRP)
if (index(Key,'C1') /= 0) then
  close(LuRP)
  return
end if
Key = Get_Ln(LuRP)
KWord = Get_Ln(LuRP)
call Get_I(1,nRP,1)
! Only read the second structure in findsym.out
do i=1,nRP+3
  read(LuRP,*)
end do
do i=1,nRP
  read(LuRP,*) KWord,j,(RP_Centers(j,i,2),j=1,3)
end do
nRP = nRP*3
close(LuRP)

! Compare with the one of the current structure
! no need if user defined

if (KG(1:4) == 'FULL') then
  if (Key /= minGroup) call Error(2)
end if

call findsymf('findsym.RP2',KG,SymThr,ierr)
if (ierr /= 0) call Error(1)
call molcas_open(LuRP,'findsym.out')
read(LuRP,*)
KWord = Get_Ln(LuRP)

if (KG(1:4) == 'FULL') then
  if (KWord /= minGroup) then
    close(LuRP)
    call Error(2)
  end if
end if

KWord = Get_Ln(LuRP)
call Get_I(1,i,1)
if (3*i /= nRP) call Error(1)
do i=1,nRP/3+3
  read(LuRP,*)
end do
do i=1,nRP/3
  read(LuRP,*) KWord,j,(RP_Centers(j,i,2),j=1,3)
end do
#else
! If manual symmetry, trust it
! Only check if automatic symmetry requested
if (KG(1:4) /= 'FULL') return
! The detected symmetry for the two structures must match
LuRP = 10
LuRP = IsFreeUnit(LuRP)
call Molcas_Open(LuRP,'findsym.RP1')
call Read_XYZ(LuRP,DumRot,DumTrans)
close(LuRP)
call Parse_Group(KeepGroup,SymThr)
nRP = Out_Raw(RP_Centers(:,:,1))
call Clear_XYZ()

KG = trim(Symmetry)
LuRP = IsFreeUnit(LuRP)
call Molcas_Open(LuRP,'findsym.RP2')
call Read_XYZ(LuRP,DumRot,DumTrans)
close(LuRP)
call Parse_Group(KeepGroup,SymThr)
i = Out_Raw(RP_Centers(:,:,2))
if (i /= nRP) call Error(1)
call Clear_XYZ()
if (Symmetry /= KG) call Error(2)
#endif

return

contains

subroutine Error(code)

  integer(kind=iwp), intent(in) :: code

  select case (code)
    case (1)
      call WarningMessage(2,'Error in RP-Coord section, check symmetry')
    case (2)
      call WarningMessage(2,'Error in RP-Coord section, structures do not have the same symmetry. '// &
                          'Please define manually the symmetry group.')
  end select
  call Quit_OnUserError()

end subroutine Error

end subroutine processRP
