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

subroutine Drv_AMFI(Label,ip,lOper,nComp,rHrmt,iChO,iAtmNr2,Charge2)

use iSD_data
use Basis_Info
use DKH_Info, only: DKroll
use Symmetry_Info, only: nIrrep

implicit real*8(a-h,o-z)
external Rsv_Tsk
#include "Molcas.fh"
#include "angtp.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "nsd.fh"
#include "setup.fh"
#include "para.fh"
integer, allocatable :: iDel(:)
real*8, allocatable :: SOInt(:)
real*8 Coor(3)
logical EQ, IfTest, Rsv_Tsk
character Label*8
integer ip(nComp), lOper(nComp), iChO(nComp)
integer iAtmNr2(mxdbsc)
real*8 Charge2(mxdbsc)
data IfTest/.false./

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
IfTest = .true.
write(6,*) ' In OneEl: Label',Label
write(6,*) ' In OneEl: nComp'
write(6,'(1X,8I5)') nComp
write(6,*) ' In OneEl: lOper'
write(6,'(1X,8I5)') lOper
write(6,*) ' In OneEl: n2Tri'
do iComp=1,nComp
  ip(iComp) = n2Tri(lOper(iComp))
end do
write(6,'(1X,8I5)') (ip(iComp),iComp=1,nComp)
#endif

Eta_Nuc = Zero

! Allocate memory for symmetry adapted one electron integrals.
! Will just store the unique elements, i.e. low triangular blocks
! and lower triangular elements in the diagonal blocks.

ip(:) = -1
LenTot = 0
do iComp=1,nComp
  ip(iComp) = 1+LenTot
  LenInt = n2Tri(lOper(iComp))
  LenTot = LenTot+LenInt+4
end do
call mma_allocate(SOInt,LenTot,label='SOInt')
SOInt(:) = Zero

! Generate list of shell information

call Nr_Shells(nSkal)

! Check that there are not several instances of the same center.

nCenter = 0
do iSkal=1,nSkal
  nCenter = max(nCenter,iSD(10,iSkal))
  if (iSD(1,iSkal) > Lmax) then
    write(6,*) ' Shells higher than '//Angtp(Lmax)//'-functions not allowed in AMFI.'
    call Quit_OnUserError()
  end if
  if ((iSD(1,iSkal) >= 2) .and. (iand(iSD(9,iSkal),1) /= 1)) then
    write(6,*) ' Only real spherical harmonics allowed'
    write(6,*) ' for AMFI.'
    call Quit_OnUserError()
  end if
end do
Coor(:) = Zero
do iCenter=1,nCenter

  ! Identify which dbsc this center belongs.

  iCnttp = 0
  do iSkal=1,nSkal
    if (iSD(10,iSkal) == iCenter) then
      iCnttp = iSD(13,iSkal)
      iCnt = iSD(14,iSkal)
      Coor(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
    end if
  end do

  ! If iCenter is not a center in the current list of shells that
  ! is to be processed then test the next iCenter.

  if (iCnttp == 0) cycle

  do iSkal=1,nSkal
    if ((iSD(13,iSkal) /= iCnttp) .and. (iSD(10,iSkal) /= iCenter)) then
      jCnttp = iSD(13,iSkal)
      jCnt = iSD(14,iSkal)
      if (EQ(Coor,dbsc(jCnttp)%Coor(1,jCnt))) then
        write(6,*) 'Multiple instances of the same center!'
        write(6,*) 'This is not allowed with AMFI.'
        call Quit_OnUserError()
      end if
    end if
  end do
end do
if ((MolWgh /= 0) .and. (MolWgh /= 2)) then
  write(6,*) ' AMFI integrals not implemented for symmetry adaptation a la MOLECULE'
  call Quit_OnUserError()
end if

! Loop over unique center. Observe that multiple shells of the same
! angular momentum is not allowed.

Lu_AMFI = 21
LUPROP = 22
call molcas_open(Lu_AMFI,'AMFI_INP')
call molcas_binaryopen_vanilla(LUPROP,'AMFI_INT')
nCenter_node = 0

call Init_Tsk(id_Tsk,nCenter)
10 continue
if (.not. Rsv_Tsk(id_Tsk,iCenter)) Go To 11
!do iCenter=1,nCenter
nCenter_node = nCenter_node+1

write(Lu_AMFI,'(A)') ' &AMFI'
if (IfTest) write(6,'(A)') ' &AMFI'

! Find atom type

mdci = 0
do iCnttp=1,nCnttp
  do iCnt=1,dbsc(iCnttp)%nCntr
    mdci = mdci+1
    if (mdci == iCenter) then
      if ((.not. DKroll) .and. (.not. dbsc(iCnttp)%SODK)) write(Lu_AMFI,'(A)') 'Breit-Pauli'
      if (IfTest) then
        if ((.not. DKroll) .and. (.not. dbsc(iCnttp)%SODK)) write(6,'(A)') 'Breit-Pauli'
      end if
      if (iAtmNr2(iCnttp) >= 1) then
        charge_x = dble(iAtmNr2(iCnttp))
      else if ((iAtmNr2(iCnttp) <= 0) .and. (Charge2(iCnttp) == Zero)) then
        charge_x = 0.0d0
      else
        write(6,*) 'Drv_AMFI: Invalid basis!'
        write(6,*) 'iAtmNr=',iAtmNr2(iCnttp)
        write(6,*) 'Charge2=',Charge2(iCnttp)
        call Abend()
      end if
      if (Nuclear_Model == Gaussian_Type) then
        Eta_Nuc = dbsc(iCnttp)%ExpNuc
      end if
      Go To 99
    end if
  end do
end do
99 continue

if (Nuclear_Model == Gaussian_Type) then
  write(Lu_AMFI,'(A)') 'Finite'
  if (IfTest) write(6,'(A)') 'Finite'
  write(Lu_AMFI,*) Eta_Nuc
  if (IfTest) write(6,*) Eta_Nuc
end if

! Generate input for each atom

l_max = -1
do iSkal=1,nSkal
  if (iSD(10,iSkal) == iCenter) l_Max = max(l_Max,iSD(1,iSkal))
end do
if (l_max > LMax) then
  write(6,*) 'AMFI integrals only implemented up to '//Angtp(Lmax)//'-functions.'
  call Quit_OnUserError()
end if

! Check if there are any core orbitals to be deleted

nCore = 0
do l=0,l_max
  do iSkal=1,nSkal
    if ((iSD(10,iSkal) == iCenter) .and. (iSD(1,iSkal) == l)) then

      iCnttp = iSD(13,iSkal)
      if (dbsc(iCnttp)%nSOC /= 0) then
        iShll = dbsc(iCnttp)%iSOC+l
        !jShll = dbsc(iCnttp)%iPrj+l
        !nCore = nCore+Shells(jShll)%nBasis
        nCore = nCore+dbsc(iCnttp)%kDel(l)
      end if
    end if
  end do
end do
if (IfTest) write(6,*) 'nCore: ',nCore

! Set up delete array

if (nCore /= 0) then
  lDel = l_Max+1
  call mma_allocate(iDel,lDel,label='iDel')
  do l=0,l_max
    do iSkal=1,nSkal
      if ((iSD(10,iSkal) == iCenter) .and. (iSD(1,iSkal) == l)) then

        iCnttp = iSD(13,iSkal)
        if (dbsc(iCnttp)%nSOC /= 0) then
          iShll = dbsc(iCnttp)%iSOC+l
          !jShll = dbsc(iCnttp)%iPrj+l
          !iDel(ip_iDel+l) = Shells(jShll)%nBasis
          iDel(1+l) = dbsc(iCnttp)%kDel(l)
        end if
      end if
    end do
  end do
  write(Lu_AMFI,'(A)') 'AIMP'
  write(Lu_AMFI,*) lDel-1,(iDel(i),i=1,lDel)
  if (IfTest) write(6,'(A)') 'AIMP'
  if (IfTest) write(6,*) lDel,(iDel(i),i=1,lDel)
  call mma_deallocate(iDel)
end if

write(Lu_AMFI,'(A)') '     '
write(Lu_AMFI,'(3X,F5.1,I4)') charge_x,l_max
if (IfTest) write(6,*) charge_x,l_max

do l=0,l_max
  do iSkal=1,nSkal
    if ((iSD(10,iSkal) == iCenter) .and. (iSD(1,iSkal) == l)) then

      iCnttp = iSD(13,iSkal)
      if (dbsc(iCnttp)%nSOC == 0) then

        ! Use valence basis

        iShll = dbsc(iCnttp)%iVal+l
        iCase = 2

      else

        ! Use special valence basis in case of a ECP where the
        ! normal valence might not be adequate.

        iShll = dbsc(iCnttp)%iSOC+l
        iCase = 1

      end if
      nBas_x = Shells(iShll)%nBasis
      nExp_x = Shells(iShll)%nExp

      if (IfTest) write(6,*) 'iShll=',iShll
      write(Lu_AMFI,*) nExp_x,nBas_x
      if (IfTest) write(6,*) nExp_x,nBas_x
      write(Lu_AMFI,*) (Shells(iShll)%exp(iExp_x),iExp_x=1,nExp_x)
      if (IfTest) write(6,*) (Shells(iShll)%exp(iExp_x),iExp_x=1,nExp_x)
      do iExp_x=1,nExp_x
        write(Lu_AMFI,*) (Shells(iShll)%Cff_c(iExp_x,iCff_x,iCase),iCff_x=1,nBas_x)
        if (IfTest) write(6,*) (Shells(iShll)%Cff_c(iExp_x,iCff_x,iCase),iCff_x=1,nBas_x)
      end do

    end if
  end do
end do
write(Lu_AMFI,'(A)') 'End of Input'
if (IfTest) write(6,'(A)') 'End of Input'

! Now call AMFI

rewind(Lu_AMFI)
call AMFI(Lu_AMFI,LUPROP,iCenter)

rewind(Lu_AMFI)

Go To 10
11 continue
call Free_Tsk(id_Tsk)
!end do
close(Lu_AMFI)

! Now symmetry adopt.

rewind(LUPROP)
call SymTrafo(LUPROP,ip,lOper,nComp,nBas,nIrrep,Label,MolWgh,SOInt,LenTot)

close(LUPROP)

call mma_deallocate(SOInt)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real(rHrmt)
  call Unused_integer_array(iChO)
end if

end subroutine Drv_AMFI
