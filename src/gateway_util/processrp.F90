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

use External_Centers
use iso_c_binding
#ifndef _HAVE_EXTRA_
use XYZ
#endif

implicit real*8(a-h,o-z)
character KeepGroup*180, KWord*180, KG*180
#ifdef _HAVE_EXTRA_
character*180 minGroup, Key
character*180, external :: Get_Ln
#else
real*8, dimension(:,:,:), allocatable :: DumRot
real*8, dimension(:,:), allocatable :: DumTrans
real*8, pointer :: p1Dim(:)
#endif
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
if (KG(1:5) == 'NOSYM') Go To 30

#ifdef _HAVE_EXTRA_
LuRP = IsFreeUnit(10)
call molcas_open(LuRP,'findsym.out')
KWord = Get_Ln(LuRP)
minGroup = Get_Ln(LuRP)
close(LuRP)
if (index(KWord,'C1') /= 0) Go To 30

! Else findsym will find the group

call findsymf('findsym.RP1',KG,SymThr,ierr)
if (ierr /= 0) Go To 20
call molcas_open(LuRP,'findsym.out')
Key = Get_Ln(LuRP)
if (index(Key,'C1') /= 0) then
  close(LuRP)
  Go To 30
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

write(6,*) Key,minGroup
if (KG(1:4) == 'FULL') then
  if (Key /= minGroup) Go To 21
end if

call findsymf('findsym.RP2',KG,SymThr,ierr)
if (ierr /= 0) Go To 20
call molcas_open(LuRP,'findsym.out')
read(LuRP,*)
KWord = Get_Ln(LuRP)

if (KG(1:4) == 'FULL') then
  if (KWord /= minGroup) then
    close(LuRP)
    Go To 21
  end if
end if

KWord = Get_Ln(LuRP)
call Get_I(1,i,1)
if (3*i /= nRP) Go To 20
do i=1,nRP/3+3
  read(LuRP,*)
end do
do i=1,nRP/3
  read(LuRP,*) KWord,j,(RP_Centers(j,i,2),j=1,3)
end do
#else
! If manual symmetry, trust it
! Only check if automatic symmetry requested
if (KG(1:4) /= 'FULL') Go To 30
! The detected symmetry for the two structures must match
LuRP = 10
LuRP = IsFreeUnit(LuRP)
call Molcas_Open(LuRP,'findsym.RP1')
call Read_XYZ(LuRP,DumRot,DumTrans)
close(LuRP)
call Parse_Group(KeepGroup,SymThr)
call c_f_pointer(c_loc(RP_Centers(1,1,1)),p1Dim,[nRP])
nRP = Out_Raw(p1Dim)
nullify(p1Dim)
call Clear_XYZ()

KG = trim(Symmetry)
LuRP = IsFreeUnit(LuRP)
call Molcas_Open(LuRP,'findsym.RP2')
call Read_XYZ(LuRP,DumRot,DumTrans)
close(LuRP)
call Parse_Group(KeepGroup,SymThr)
call c_f_pointer(c_loc(RP_Centers(1,1,2)),p1Dim,[nRP])
i = Out_Raw(p1Dim)
nullify(p1Dim)
if (i /= nRP) Go To 20
call Clear_XYZ()
if (Symmetry /= KG) Go To 21
#endif
30 continue

return

20 continue
call WarningMessage(2,'Error in RP-Coord section, check symmetry')
call Quit_OnUserError()
21 continue
KWord = 'Error in RP-Coord section, structures do not have the same symmetry. Please define manually the symmetry group.'
call WarningMessage(2,KWord)
call Quit_OnUserError()

end subroutine processRP
