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

!#define _DEBUGPRINT_
subroutine Def_Shells(iSD,nSD,mSkal)

use Index_Functions, only: nTri_Elem1
use Basis_Info, only: dbsc, iCnttp_dummy, nBas, nCnttp, Shells
use Sizes_of_Seward, only: S
use BasisMode, only: All_Mode, Atomic, Auxiliary_Mode, Basis_Mode, Fragment_Mode, kCnttp, Valence_Mode, With_Auxiliary_Mode, &
                     With_Fragment_Mode
use disp, only: Dirct, IndDsp
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nSD, mSkal
integer(kind=iwp), intent(out) :: iSD(0:nSD,mSkal)
integer(kind=iwp) :: iAng, iAOttp, iCar, iCase, iCmp, iCnt, iCnttp, iComp, iIrrep, iShell, iShell_Set, iShll, iTemp, iTmp, jCnttp, &
                     kSh, mdc, mdci, nBasisi, nDIsp, nExpi, nFunctions, nSkal, ntest
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, j
#endif
logical(kind=iwp), external :: TF

!                                                                      *
!***********************************************************************
!                                                                      *
if ((Basis_Mode /= Valence_Mode) .and. (Basis_Mode /= Auxiliary_Mode) .and. (Basis_Mode /= Fragment_Mode) .and. &
    (Basis_Mode /= With_Auxiliary_Mode) .and. (Basis_Mode /= With_Fragment_Mode) .and. (Basis_Mode /= All_Mode)) then
  call WarningMessage(2,'Def_Shells: Basis_Mode is not defined')
  call Abend()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
iIrrep = 0
nSkal = 0
iAOttp = 0 ! Number of AO functions proceeding a particular shell
S%m2Max = 0

if (.not. Atomic) then
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Molecular setup                                                    *
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  iCnttp = 0
  mdc = 0
  iShell = 0
  do jCnttp=1,nCnttp

    ! Make sure that we process the dummy shell last

    if ((jCnttp == iCnttp_Dummy) .and. (jCnttp /= nCnttp)) then
      iCnttp = iCnttp+2
    else if ((jCnttp == nCnttp) .and. (iCnttp == jCnttp)) then
      iCnttp = iCnttp_Dummy
    else
      iCnttp = iCnttp+1
    end if

    nTest = dbsc(iCnttp)%nVal-1
    mdci = dbsc(iCnttp)%mdci
    do iCnt=1,dbsc(iCnttp)%nCntr
      mdci = mdci+1
      mdc = mdc+1
      iShell_Set = iShell

      do iAng=0,nTest
        iShell = iShell+1
        iShll = dbsc(iCnttp)%iVal+iAng
        nExpi = Shells(iShll)%nExp
        nBasisi = Shells(iShll)%nBasis
        if (Shells(iShll)%Prjct) then
          iCmp = 2*iAng+1
        else
          iCmp = nTri_Elem1(iAng)
        end if
        if (nExpi == 0) cycle
        if (nBasisi == 0) cycle
        if ((Basis_Mode == Valence_Mode) .and. (Shells(iShll)%Aux .or. Shells(iShll)%Frag)) cycle
        if ((Basis_Mode == Auxiliary_Mode) .and. (.not. Shells(iShll)%Aux)) cycle
        if ((Basis_Mode == Fragment_Mode) .and. (.not. Shells(iShll)%Frag)) cycle
        if ((Basis_Mode == With_Auxiliary_Mode) .and. Shells(iShll)%Frag) cycle
        if ((Basis_Mode == With_Fragment_Mode) .and. Shells(iShll)%Aux) cycle

        kSh = dbsc(iCnttp)%iVal+iAng

        nSkal = nSkal+1
        !                                                              *
        !***************************************************************
        !                                                              *
        iSD(0,nSkal) = iShll                  ! Unique shell ind.
        iSD(1,nSkal) = iAng                   ! l value
        iSD(2,nSkal) = iCmp                   ! # of ang. comp.
        iSD(3,nSkal) = nBasisi                ! # of cont. func.
        iSD(4,nSkal) = -1                     ! Not used
        iSD(5,nSkal) = nExpi                  ! # of prim.
        iSD(6,nSkal) = -1                     ! Not used
        iSD(7,nSkal) = iAOttp+(iCnt-1)*dbsc(iCnttp)%lOffAO+Shells(kSh)%kOffAO
        iSD(8,nSkal) = -1
        itemp = 0
        if (Shells(iShll)%Prjct) itemp = itemp+1
        if (Shells(iShll)%Transf) itemp = itemp+2
        iSD(9,nSkal) = itemp                  ! sph., car., cont.
        iSD(10,nSkal) = mdci                  ! Center index
        iSD(11,nSkal) = iShell_Set+iAng+1
        if (dbsc(iCnttp)%pChrg) then
          iSD(12,nSkal) = 1                   ! pseudo charge
        else
          iSD(12,nSkal) = 0                   ! pseudo charge
        end if
        iSD(13,nSkal) = iCnttp
        iSD(14,nSkal) = iCnt

        nDisp = IndDsp(mdci,iIrrep)
        iTmp = 0
        do iCar=0,2
          iComp = 2**iCar
          if (TF(mdci,iIrrep,iComp) .and. (.not. dbsc(iCnttp)%pChrg)) then
            nDisp = nDisp+1
            if (Dirct(nDisp)) then
              iSD(iCar+16,nSkal) = nDisp
              iTmp = ibset(iTmp,iCar)
            else
              iSD(iCar+16,nSkal) = 0
            end if
          else
            iSD(iCar+16,nSkal) = 0
          end if
        end do
        iSD(15,nSkal) = iTmp
        iSD(20,nSkal) = nSkal

        S%m2Max = max(S%m2Max,nExpi**2)
        !                                                              *
        !***************************************************************
        !                                                              *
      end do  ! iAng
    end do    ! iCnt
    iAOttp = iAOttp+dbsc(iCnttp)%lOffAO*dbsc(iCnttp)%nCntr
  end do      ! iCnttp

else
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Atomic setup                                                       *
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *

  iCase = 1
  mdci = 1
  iCnt = 1
  nFunctions = 0

  do
    if (iCase == 1) then

      ! Add the auxiliary basis

      iCnttp = kCnttp
    else

      ! Add the dummy basis

      iCnttp = iCnttp_Dummy
    end if
    nTest = dbsc(iCnttp)%nVal-1

    do iAng=0,nTest
      iShll = dbsc(iCnttp)%iVal+iAng
      nExpi = Shells(iShll)%nExp
      if (nExpi == 0) cycle
      nBasisi = Shells(iShll)%nBasis
      if (nBasisi == 0) cycle
      if (Shells(iShll)%Frag) cycle
      iCmp = nTri_Elem1(iAng)
      if (Shells(iShll)%Prjct) iCmp = 2*iAng+1
      kSh = dbsc(iCnttp)%iVal+iAng

      nSkal = nSkal+1
      !                                                                *
      !*****************************************************************
      !                                                                *
      iSD(0,nSkal) = iShll                  ! Unique shell ind.
      iSD(1,nSkal) = iAng                   ! l value
      iSD(2,nSkal) = iCmp                   ! # of ang. comp.
      iSD(3,nSkal) = nBasisi                ! # of cont. func.
      iSD(4,nSkal) = -1                     ! Not used
      iSD(5,nSkal) = nExpi                  ! # of prim.
      iSD(6,nSkal) = -1                     ! Not used
      iSD(7,nSkal) = iAOttp+Shells(kSh)%kOffAO
      iSD(8,nSkal) = -1                     ! Not used
      itemp = 0
      if (Shells(iShll)%Prjct) itemp = itemp+1
      if (Shells(iShll)%Transf) itemp = itemp+2
      iSD(9,nSkal) = itemp                  ! sph., car., cont.
      iSD(10,nSkal) = mdci                  ! Center index
      iSD(11,nSkal) = iAng+1                ! Not used
      if (dbsc(iCnttp)%pChrg) then
        iSD(12,nSkal) = 1                   ! pseudo charge
      else
        iSD(12,nSkal) = 0                   ! pseudo charge
      end if
      iSD(13,nSkal) = iCnttp
      iSD(14,nSkal) = 1

      iSD(15,nSkal) = 0
      iSD(16,nSkal) = 0
      iSD(17,nSkal) = 0
      iSD(18,nSkal) = 0
      iSD(19,nSkal) = 0
      iSD(20,nSkal) = nSkal

      S%m2Max = max(S%m2Max,nExpi**2)

      if (Shells(iShll)%Prjct) then
        nFunctions = nFunctions+nBasisi*(2*iAng+1)
      else
        nFunctions = nFunctions+nBasisi*nTri_Elem1(iAng)
      end if
      !                                                                *
      !*****************************************************************
      !                                                                *
    end do      !iAng

    iCase = iCase+1
    if ((iCase > 2) .or. (.not. dbsc(iCnttp)%Aux)) exit
  end do

  if (dbsc(iCnttp)%Aux) then
    nBas(0) = 0
  else
    nBas(0) = nFunctions
  end if
# ifdef _DEBUGPRINT_
  write(u6,*) 'in Define_Shells...'
  do i=1,mSkal
    write(u6,*) 'i=',i
    write(u6,'(10I8,/,8I8)') (iSD(j,i),j=0,nSD)
  end do
# endif
end if
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *

end subroutine Def_Shells
