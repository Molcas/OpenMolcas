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

subroutine AllenGinsberg(QMMethod,Eint,Poli,dNuc,Cha,Dip,Qua,MxBaux,iVEC,nDim,iExtr_Atm,lEig,iEig,iQ_Atoms,ip_ExpCento,E_Nuc_Part, &
                         lSlater,Eint_Nuc)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
#include "warnings.h"
dimension Eint(MxQCen,10), Poli(MxQCen,10)
dimension Cha(MxBaux,MxQCen), Dip(MxBaux,3,MxQCen)
dimension Qua(MxBaux,6,MxQCen), dNuc(MxAt)
dimension iExtr_Atm(MxAt), Eint_Nuc(MxAt)
dimension iCenSet(MxAt**2)
character QMMethod*5
logical lEig, Check1, Check2, lSlater

! Set up centre index set. The order of centres are decided by
! the MpProp-program and are hence collected in the get_center
! routine.

! Atom centres

do iAt=1,MxAt
  if (iExtr_Atm(iAt) == -1) then
    goto 902
  else
    iCenSet(iAt) = iExtr_Atm(iAt)
  end if
end do
NExtrAt = iAt
902 continue
NExtrAt = iAt-1

! Bond centres

kaunter = iQ_Atoms
kaunt = NExtrAt
do iAt=2,iQ_Atoms
  do jAt=1,iAt-1
    kaunter = kaunter+1
    Check1 = .false.
    Check2 = .false.
    do i1=1,NExtrAt
      if (iAt == iCenSet(i1)) Check1 = .true.
      if (jAt == iCenSet(i1)) Check2 = .true.
    end do
    if (Check1 .and. Check2) then
      kaunt = kaunt+1
      iCenSet(kaunt) = kaunter
    end if
  end do
end do

! A minor check.

NExpect = NExtrAt*(nExtrAt+1)/2
NTotal = kaunt
if (NTotal /= NExpect) then
  write(6,*)
  write(6,*) ' Error in atom specification for partial perturbation extraction.'
  call Quit(_RC_GENERAL_ERROR_)
end if

! Compute partial nuclear contribution.

E_Nuc_Part = 0.0d0
do iAt=1,NExtrAt
  iCx = iCenSet(iAt)
  if (lSlater) then
    E_Nuc_Part = E_Nuc_Part-(Eint_Nuc(iCx)+Poli(iCx,1))*dNuc(iCx)
  else
    E_Nuc_Part = E_Nuc_Part-(Eint(iCx,1)+Poli(iCx,1))*dNuc(iCx)
  end if
end do

! Set up matrix elements for the partial perturbations.
! Compare with hel, helstate, polink and polins.

call GetMem('VelPart','Allo','Real',iVelP,nDim*(nDim+1)/2)
call GetMem('VpoPart','Allo','Real',iVpoP,nDim*(nDim+1)/2)
call dcopy_(nDim*(nDim+1)/2,[ZERO],iZERO,Work(iVelP),iONE)
call dcopy_(nDim*(nDim+1)/2,[ZERO],iZERO,Work(iVpoP),iONE)
kk = 0
do i=1,nDim
  do j=1,i
    kk = kk+1
    do k=1,NTotal
      iCx = iCenSet(k)
      dMp = Cha(kk,iCx)
      Work(iVelP+kk-1) = Work(iVelP+kk-1)+Eint(iCx,1)*dMp
      Work(iVpoP+kk-1) = Work(iVpoP+kk-1)+Poli(iCx,1)*dMp
      dMp = Dip(kk,1,iCx)
      Work(iVelP+kk-1) = Work(iVelP+kk-1)+Eint(iCx,2)*dMp
      Work(iVpoP+kk-1) = Work(iVpoP+kk-1)+Poli(iCx,2)*dMp
      dMp = Dip(kk,2,iCx)
      Work(iVelP+kk-1) = Work(iVelP+kk-1)+Eint(iCx,3)*dMp
      Work(iVpoP+kk-1) = Work(iVpoP+kk-1)+Poli(iCx,3)*dMp
      dMp = Dip(kk,3,iCx)
      Work(iVelP+kk-1) = Work(iVelP+kk-1)+Eint(iCx,4)*dMp
      Work(iVpoP+kk-1) = Work(iVpoP+kk-1)+Poli(iCx,4)*dMp
      dMp = Qua(kk,1,iCx)
      Work(iVelP+kk-1) = Work(iVelP+kk-1)+Eint(iCx,5)*dMp
      Work(iVpoP+kk-1) = Work(iVpoP+kk-1)+Poli(iCx,5)*dMp
      dMp = Qua(kk,3,iCx)
      Work(iVelP+kk-1) = Work(iVelP+kk-1)+Eint(iCx,7)*dMp
      Work(iVpoP+kk-1) = Work(iVpoP+kk-1)+Poli(iCx,7)*dMp
      dMp = Qua(kk,6,iCx)
      Work(iVelP+kk-1) = Work(iVelP+kk-1)+Eint(iCx,10)*dMp
      Work(iVpoP+kk-1) = Work(iVpoP+kk-1)+Poli(iCx,10)*dMp
      dMp = Qua(kk,2,iCx)
      Work(iVelP+kk-1) = Work(iVelP+kk-1)+Eint(iCx,6)*dMp*2.0d0
      Work(iVpoP+kk-1) = Work(iVpoP+kk-1)+Poli(iCx,6)*dMp*2.0d0
      dMp = Qua(kk,4,iCx)
      Work(iVelP+kk-1) = Work(iVelP+kk-1)+Eint(iCx,8)*dMp*2.0d0
      Work(iVpoP+kk-1) = Work(iVpoP+kk-1)+Poli(iCx,8)*dMp*2.0d0
      dMp = Qua(kk,5,iCx)
      Work(iVelP+kk-1) = Work(iVelP+kk-1)+Eint(iCx,9)*dMp*2.0d0
      Work(iVpoP+kk-1) = Work(iVpoP+kk-1)+Poli(iCx,9)*dMp*2.0d0
    end do
  end do
end do

! Collect expectation value for the partial perturbation.

call Expectus(QMMethod,Work(iVelP),Work(iVelP),Work(iVpoP),Work(iVpoP),MxBaux,iVEC,nDim,lEig,iEig,ip_ExpCento)

! Deallocate.

call GetMem('VelPart','Free','Real',iVelP,nDim*(nDim+1)/2)
call GetMem('VpoPart','Free','Real',iVpoP,nDim*(nDim+1)/2)

! Howl

return

end subroutine AllenGinsberg
