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

subroutine Get_Molecule()

use Slapaf_Info, only: Q_nuclear, Coor, Grd, Weights, AtomLbl
use Symmetry_Info, only: VarR, VarT

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "Molcas.fh"
logical Found
integer Columbus

!                                                                      *
!***********************************************************************
!                                                                      *
! Read initial data

call Get_iScalar('Unique atoms',nsAtom)

call mma_allocate(Coor,3,nsAtom,Label='Coor')
call Get_dArray('Unique Coordinates',Coor,3*nsAtom)

call mma_allocate(Q_nuclear,nsAtom)
call Get_dArray('Nuclear charge',Q_nuclear,nsAtom)

call Get_iScalar('Grad ready',iGO)
iJustGrad = iand(iGO,2**0)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate gradient (it will be read later)
! (This should eventually be removed, as Grd is unused...)

call Get_iScalar('Columbus',columbus)
if ((iJustGrad == 1) .and. (columbus == 1)) then

  ! C&M mode

  call Get_iScalar('ColGradMode',iMode)
  if (iMode == 0) then
    call mma_allocate(Grd,3,nsAtom,Label='Grd')
    call Get_dArray_chk('GRAD',Grd,3*nsAtom)
  else if (iMode <= 3) then
    call qpg_dArray('Grad State1',Found,Length)
    if (.not. Found .or. (Length == 0)) then
      call SysAbendmsg('Get_Molecule','Did not find:','Grad State1')
    end if
    if (length /= 3*nsAtom) then
      call WarningMessage(2,'Init: length /= 3*nsAtom')
      write(6,*) 'Grad'
      write(6,*) 'length,nsAtom=',length,nsAtom
      call Abend()
    end if
    call mma_allocate(Grd,3,nsAtom,Label='Grd')
    call Get_dArray('Grad State1',Grd,3*nsAtom)

  end if
  iJustGrad = 0
  iGO = ior(iGO,iJustGrad)
  call Put_iScalar('Grad ready',iGO)
else

  ! M mode

  call mma_allocate(Grd,3,nsAtom,Label='Grd')
  Grd(:,:) = Zero
end if

call mma_allocate(AtomLbl,nsAtom,Label='AtomLbl')
call Get_cArray('Unique Atom Names',AtomLbl,LenIn*nsAtom)
!                                                                      *
!***********************************************************************
!                                                                      *
! Check if method is translational or rotational variant.

iPL = iPrintLevel(-1)
if ((VarT .or. VarR) .and. (iPL > 0)) then
  write(6,*)
  if (VarT) write(6,*) '    Gradient is translational variant!'
  if (VarR) write(6,*) '    Gradient is rotational variant!'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Read weights

call Qpg_dArray('Weights',Found,nData)
if (Found .and. (nData >= nsAtom)) then
  ! The weights array length is actually the total number of atoms,
  ! not just symmetry-unique, but the symmetry-unique ones are first
  call mma_allocate(Weights,nData,Label='Weights')
  call Get_dArray('Weights',Weights,nData)
else
  call SysAbendMsg('Get_Molecule','No or wrong weights were found in the RUNFILE.','')
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Get_Molecule
