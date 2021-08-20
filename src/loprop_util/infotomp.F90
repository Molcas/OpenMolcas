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

subroutine InfoToMp(nSym,nBas,Energy_Without_FFPT,ip_Ene_Occ,nOcOb,UserDen,Restart)

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
integer nBas(8)
character*8 Method
character*40 VTitle
character*6 FName
logical UserDen, Restart
dimension iDum(1)

!                                                                      *
!***********************************************************************
!                                                                      *
! Set some variables needed for creating the MpProp file.              *
!                                                                      *
!***********************************************************************
!                                                                      *
NOCOB = 0
if (.not. UserDen) then
  nVec = 0
  nOcc = 0
  do iSym=1,nSym
    nVec = nVec+nBas(iSym)**2
    nOcc = nOcc+nBas(iSym)
  end do
  call Allocate_Work(ip_Ene_Occ,nOcc)
  if (Restart) then
    call Get_dScalar('MpProp Energy',Energy_Without_FFPT)
    call Get_dArray('MpProp Orb Ener',Work(ip_Ene_occ),nOcc)
    call Get_iScalar('MpProp nOcOb',nOcOb)
  else
    call Get_dScalar('Last energy',Energy_Without_FFPT)
    call Put_dScalar('MpProp Energy',Energy_Without_FFPT)

    call Allocate_Work(ip_Vec,nVec)
    call Allocate_Work(ip_Occ,nOcc)

    Lu_ = 11
    FName = 'INPORB'
    iDum = 0
    iWarn = 2
    call RdVec(FName,Lu_,'COE',nSym,nBas,nBas,Work(ip_Vec),Work(ip_Occ),Work(ip_Ene_Occ),iDum,VTitle,iWarn,iErr)
    close(Lu_)

    do i=0,nOcc-1
      if (Work(ip_Occ+i) /= 0.0d0) then
        nOcOb = nOcOb+1
      end if
    end do
    call Put_dArray('MpProp Orb Ener',Work(ip_Ene_occ),nOcc)
    call Put_iScalar('MpProp nOcOb',nOcOb)

    call Free_Work(ip_Vec)
    call Free_Work(ip_Occ)
  end if
else
  !Here we go if user give density as input. Need to put some
  !dummy values for later; also we give a value to the relax
  !method, which for now is called 'external'.
  Energy_Without_FFPT = 0.0d0
  nOcc = 0
  do iSym=1,nSym
    nOcc = nOcc+nBas(iSym)
  end do
  call Allocate_Work(ip_Ene_Occ,nOcc)
  do i=0,nOcc-1
    Work(ip_Ene_Occ+i) = 0.0d0
  end do
  write(Method,'(A)') 'External'
  call Put_cArray('Relax Method',Method,8)
end if

return

end subroutine InfoToMp
