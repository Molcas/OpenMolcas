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

subroutine Chk_Numerical(LuSpool,Numerical)

use Alaska_Info, only: Auto, DefRoot, ForceNAC, iRlxRoot
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuSpool
logical(kind=iwp), intent(out) :: Numerical
#include "nac.fh"
integer(kind=iwp) :: iDNG, iGO, iRoot, iRoot0, istatus, LuWr
real(kind=wp) :: rDelta
logical(kind=iwp) :: Is_Root_Set, DNG, KeepOld, Found
character(len=180) :: KWord, Key
character(len=180), external :: Get_Ln

call Qpg_iScalar('DNG',DNG)
if (DNG) then
  call Get_iScalar('DNG',iDNG)
  Numerical = iDNG == 1
else
  Numerical = .false.
end if
LuWr = u6

! Setting the defaults
!    iRoot     : Which root to optimize the geometry for
!    rDelta    : Displacements are chosen as r_nearest_neighbor * rDelta
iRoot = 1
rDelta = 0.01_wp
NACstates(1) = 0
NACstates(2) = 0
isNAC = .false.
KeepOld = .false.
DefRoot = .true.
ForceNAC = .false.
iRlxRoot = 1
Auto = .false.

! RASSCF will set the default root to the root that is relaxed.

Is_Root_Set = .false.
call Qpg_iScalar('NumGradRoot',Is_Root_Set)
if (Is_Root_Set) then
  call Get_iScalar('NumGradRoot',iRoot)
end if

rewind(LuSpool)
call RdNLst(LuSpool,'ALASKA')
KWord = ' &ALASKA'
do
  read(LuSpool,'(A72)',iostat=istatus) Key
  if (istatus > 0) then
    call WarningMessage(2,'Chk_Numerical: Error reading the input')
    write(LuWr,'(A,A)') 'Last read line=',KWord
    call Quit_OnUserError()
  else if (istatus < 0) then
    exit
  end if
  KWord = Key
  call UpCase(KWord)
  select case (KWord(1:4))
    case ('NUME')
      Numerical = .true.
    case ('ROOT')
      Key = Get_Ln(LuSpool)
      call Get_I1(1,iRoot)
      DefRoot = .false.
    case ('DELT')
      Key = Get_Ln(LuSpool)
      call Get_F1(1,rDelta)
    case ('NAC ')
      Key = Get_Ln(LuSpool)
      call Get_I(1,NACstates,2)
      isNAC = .true.
      DefRoot = .false.
    case ('KEEP')
      KeepOld = .true.
    case ('AUTO')
      Auto = .true.
    case ('END ')
      exit
    case default
  end select
end do

call Get_iScalar('Grad ready',iGO)
iGO = ibclr(iGo,0)
call Put_iScalar('Grad ready',iGO)

! Put on the runfile which root and delta to use

call qpg_iScalar('Relax CASSCF root',Found)
if (Found) then
  call Get_iScalar('Relax CASSCF root',iRoot0)
  call Put_iScalar('NumGradRoot',iRoot)
  call Put_iScalar('Relax CASSCF root',iRoot)
else
  iRoot0 = 0
end if
call Put_dScalar('Numerical Gradient rDelta',rDelta)

! These are really input options for numerical_gradient

call Put_lScalar('Keep old gradient',KeepOld)

return

end subroutine Chk_Numerical
