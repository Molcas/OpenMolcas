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

subroutine Drv_AMFI(Label,lOper,nComp,iAtmNr2,Charge2)

use AMFI_global, only: Lmax
use iSD_data, only: iSD
use Basis_Info, only: dbsc, Gaussian_Type, MolWgh, nBas, nCnttp, Nuclear_Model, Shells
use DKH_Info, only: DKroll
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
character(len=8), intent(in) :: Label
integer(kind=iwp), intent(in) :: nComp, lOper(nComp), iAtmNr2(mxdbsc)
real(kind=wp), intent(in) :: Charge2(mxdbsc)
#include "angtp.fh"
integer(kind=iwp) :: i, iCase, iCenter, iCff_x, iCnt, iCnttp, id_Tsk, iExp_x, iShll, iSkal, jCnt, jCnttp, l, l_max, lDel, &
                     Lu_AMFI, LUPROP, mdci, nBas_x, nCenter, nCenter_node, nCore, nExp_x, nSkal
real(kind=wp) :: charge_x, Coor(3), Eta_Nuc
logical(kind=iwp) :: EQ
integer(kind=iwp), allocatable :: iDel(:)
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
#define _TEST_ .true.
#else
#define _TEST_ .false.
#endif
logical(kind=iwp), parameter :: IfTest = _TEST_
logical(kind=iwp), external :: Rsv_Tsk

#ifdef _DEBUGPRINT_
write(u6,*) ' In OneEl: Label',Label
write(u6,*) ' In OneEl: nComp'
write(u6,'(1X,8I5)') nComp
write(u6,*) ' In OneEl: lOper'
write(u6,'(1X,8I5)') lOper
#endif

Eta_Nuc = Zero

! Generate list of shell information

call Nr_Shells(nSkal)

! Check that there are not several instances of the same center.

nCenter = 0
do iSkal=1,nSkal
  nCenter = max(nCenter,iSD(10,iSkal))
  if (iSD(1,iSkal) > Lmax) then
    write(u6,*) ' Shells higher than '//Angtp(Lmax)//'-functions not allowed in AMFI.'
    call Quit_OnUserError()
  end if
  if ((iSD(1,iSkal) >= 2) .and. (.not. btest(iSD(9,iSkal),0))) then
    write(u6,*) ' Only real spherical harmonics allowed'
    write(u6,*) ' for AMFI.'
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
        write(u6,*) 'Multiple instances of the same center!'
        write(u6,*) 'This is not allowed with AMFI.'
        call Quit_OnUserError()
      end if
    end if
  end do
end do
if ((MolWgh /= 0) .and. (MolWgh /= 2)) then
  write(u6,*) ' AMFI integrals not implemented for symmetry adaptation a la MOLECULE'
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
do
  if (.not. Rsv_Tsk(id_Tsk,iCenter)) exit
  !do iCenter=1,nCenter
  nCenter_node = nCenter_node+1

  write(Lu_AMFI,'(A)') ' &AMFI'
  if (IfTest) write(u6,'(A)') ' &AMFI'

  ! Find atom type

  mdci = 0
  outer: do iCnttp=1,nCnttp
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdci = mdci+1
      if (mdci == iCenter) then
        if ((.not. DKroll) .and. (.not. dbsc(iCnttp)%SODK)) write(Lu_AMFI,'(A)') 'Breit-Pauli'
        if (IfTest) then
          if ((.not. DKroll) .and. (.not. dbsc(iCnttp)%SODK)) write(u6,'(A)') 'Breit-Pauli'
        end if
        if (iAtmNr2(iCnttp) >= 1) then
          charge_x = real(iAtmNr2(iCnttp),kind=wp)
        else if ((iAtmNr2(iCnttp) <= 0) .and. (Charge2(iCnttp) == Zero)) then
          charge_x = Zero
        else
          write(u6,*) 'Drv_AMFI: Invalid basis!'
          write(u6,*) 'iAtmNr=',iAtmNr2(iCnttp)
          write(u6,*) 'Charge2=',Charge2(iCnttp)
          call Abend()
        end if
        if (Nuclear_Model == Gaussian_Type) Eta_Nuc = dbsc(iCnttp)%ExpNuc
        exit outer
      end if
    end do
  end do outer

  if (Nuclear_Model == Gaussian_Type) then
    write(Lu_AMFI,'(A)') 'Finite'
    if (IfTest) write(u6,'(A)') 'Finite'
    write(Lu_AMFI,*) Eta_Nuc
    if (IfTest) write(u6,*) Eta_Nuc
  end if

  ! Generate input for each atom

  l_max = -1
  do iSkal=1,nSkal
    if (iSD(10,iSkal) == iCenter) l_Max = max(l_Max,iSD(1,iSkal))
  end do
  if (l_max > LMax) then
    write(u6,*) 'AMFI integrals only implemented up to '//Angtp(Lmax)//'-functions.'
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
  if (IfTest) write(u6,*) 'nCore: ',nCore

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
    if (IfTest) write(u6,'(A)') 'AIMP'
    if (IfTest) write(u6,*) lDel,(iDel(i),i=1,lDel)
    call mma_deallocate(iDel)
  end if

  write(Lu_AMFI,'(A)') '     '
  write(Lu_AMFI,'(3X,F5.1,I4)') charge_x,l_max
  if (IfTest) write(u6,*) charge_x,l_max

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

        if (IfTest) write(u6,*) 'iShll=',iShll
        write(Lu_AMFI,*) nExp_x,nBas_x
        if (IfTest) write(u6,*) nExp_x,nBas_x
        write(Lu_AMFI,*) (Shells(iShll)%Exp(iExp_x),iExp_x=1,nExp_x)
        if (IfTest) write(u6,*) (Shells(iShll)%Exp(iExp_x),iExp_x=1,nExp_x)
        do iExp_x=1,nExp_x
          write(Lu_AMFI,*) (Shells(iShll)%Cff_c(iExp_x,iCff_x,iCase),iCff_x=1,nBas_x)
          if (IfTest) write(u6,*) (Shells(iShll)%Cff_c(iExp_x,iCff_x,iCase),iCff_x=1,nBas_x)
        end do

      end if
    end do
  end do
  write(Lu_AMFI,'(A)') 'End of Input'
  if (IfTest) write(u6,'(A)') 'End of Input'

  ! Now call AMFI

  rewind(Lu_AMFI)
  call AMFI(Lu_AMFI,LUPROP,iCenter)

  rewind(Lu_AMFI)

end do
call Free_Tsk(id_Tsk)
!end do
close(Lu_AMFI)

! Now symmetry adapt.

rewind(LUPROP)

call SymTrafo(LUPROP,lOper,nComp,nBas,nIrrep,Label,MolWgh)

close(LUPROP)

return

end subroutine Drv_AMFI
