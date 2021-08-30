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

subroutine InfoToMp(nSym,nBas,nOcc,Energy_Without_FFPT,Ene_Occ,nOcOb,UserDen,Restart)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(8), nOcc
real(kind=wp), intent(out) :: Energy_Without_FFPT, Ene_Occ(nOcc)
integer(kind=iwp), intent(out) :: nOcOb
logical(kind=iwp), intent(in) :: UserDen, Restart
integer(kind=iwp) :: i, iDum(1), iErr, iSym, iWarn, Lu_, nVec
character(len=40) :: VTitle
character(len=8) :: Method
character(len=6) :: FName
real(kind=wp), allocatable :: Occ(:), Vec(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Set some variables needed for creating the MpProp file.              *
!                                                                      *
!***********************************************************************
!                                                                      *
nOcOb = 0
if (.not. UserDen) then
  nVec = 0
  do iSym=1,nSym
    nVec = nVec+nBas(iSym)**2
  end do
  if (Restart) then
    call Get_dScalar('MpProp Energy',Energy_Without_FFPT)
    call Get_dArray('MpProp Orb Ener',Ene_occ,nOcc)
    call Get_iScalar('MpProp nOcOb',nOcOb)
  else
    call Get_dScalar('Last energy',Energy_Without_FFPT)
    call Put_dScalar('MpProp Energy',Energy_Without_FFPT)

    call mma_allocate(Vec,nVec,label='Vec')
    call mma_allocate(Occ,nOcc,label='Occ')

    Lu_ = 11
    FName = 'INPORB'
    iDum = 0
    iWarn = 2
    call RdVec(FName,Lu_,'COE',nSym,nBas,nBas,Vec,Occ,Ene_Occ,iDum,VTitle,iWarn,iErr)
    close(Lu_)

    do i=1,nOcc
      if (Occ(i) /= Zero) nOcOb = nOcOb+1
    end do
    call Put_dArray('MpProp Orb Ener',Ene_occ,nOcc)
    call Put_iScalar('MpProp nOcOb',nOcOb)

    call mma_deallocate(Vec)
    call mma_deallocate(Occ)
  end if
else
  !Here we go if user give density as input. Need to put some
  !dummy values for later; also we give a value to the relax
  !method, which for now is called 'external'.
  Energy_Without_FFPT = Zero
  Ene_Occ(:) = Zero
  write(Method,'(A)') 'External'
  call Put_cArray('Relax Method',Method,8)
end if

return

end subroutine InfoToMp
