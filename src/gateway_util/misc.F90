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

subroutine Misc_Seward(iBas,iBas_Aux,iBas_Frag)

use Basis_Info, only: dbsc, iCnttp_Dummy, nCnttp, Shells
use Center_Info, only: dc
use Sizes_of_Seward, only: S
use Gateway_Info, only: RadMax, cdMax, EtMax
use Symmetry_Info, only: nIrrep
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iBas, iBas_Aux, iBas_Frag
#include "Molcas.fh"
integer(kind=iwp) :: iAng, iCnt, iCnttp, iShell, jCnttp, jSh, kCmp, kdc, mc, mdc

!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over distinct shell types

iBas = 0
iBas_Aux = 0
iBas_Frag = 0

iShell = 0
mc = 1
kdc = 0
! Loop over basis sets
iCnttp = 0
do jCnttp=1,nCnttp

  ! Make sure that we process the dummy shell last

  if ((jCnttp == iCnttp_Dummy) .and. (jCnttp /= nCnttp)) then
    iCnttp = iCnttp+2
  else if ((jCnttp == nCnttp) .and. (iCnttp == jCnttp)) then
    iCnttp = iCnttp_Dummy
  else
    iCnttp = iCnttp+1
  end if

  ! Loop over distinct centers

  do iCnt=1,dbsc(iCnttp)%nCntr
    kdc = kdc+1
    mdc = iCnt+dbsc(iCnttp)%mdci
    if (max(mdc,kdc) > MxAtom) then
      call WarningMessage(2,'MxAtom too small:')
      write(u6,*) 'MxAtom=',MxAtom
      write(u6,*) 'Increase mxAtom in Molcas.fh and recompile the code!'
      call Abend()
    end if
    ! Loop over shells associated with this center
    ! Start with s type shells
    jSh = dbsc(iCnttp)%iVal
    do iAng=0,dbsc(iCnttp)%nVal-1
      iShell = iShell+1

      if (Shells(jSh)%nBasis_C > 0) call RdMx(RadMax,Shells(jSh)%Exp,Shells(jSh)%nExp,Shells(jSh)%Cff_c(1,1,1), &
                                              Shells(jSh)%nBasis_C,cdMax,EtMax)
      if (iShell > MxShll) then
        call WarningMessage(2,'iShell > MxShll; Change MxShll in Molcas.fh and recompile the code!')
        call Abend()
      end if
      kCmp = (iAng+1)*(iAng+2)/2
      if (Shells(jSh)%Prjct) kCmp = 2*iAng+1

      if (Shells(jSh)%nBasis /= 0) then
        if (Shells(jSh)%Aux) then
          iBas_Aux = iBas_Aux+Shells(jSh)%nBasis*kCmp*nIrrep/dc(mdc)%nStab
        else if (Shells(jSh)%Frag) then
          iBas_Frag = iBas_Frag+Shells(jSh)%nBasis*kCmp*nIrrep/dc(mdc)%nStab
        else
          iBas = iBas+Shells(jSh)%nBasis*kCmp*nIrrep/dc(mdc)%nStab
        end if
      end if
      jSh = jSh+1
    end do
    mc = mc+nIrrep/dc(mdc)%nStab
  end do

end do
S%nShlls = iShell
if (iBas >= 2*MaxBfn) then
  call WarningMessage(2,'MaxBfn too small')
  write(u6,*) 'Increase 2*MaxBfn to ',iBas
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Misc_Seward
