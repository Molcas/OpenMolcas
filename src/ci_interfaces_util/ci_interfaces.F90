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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

module CI_Interfaces

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
private

public :: Mk_H_Psi, Mk_pdms

contains

subroutine Mk_H_Psi(SGS,EXS,CIS,STSYM,nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,ndeta,ndetb,nTU,TU,nTUVX,TUVX)

  use Lucia_Interface, only: Lucia_Util
  use lucia_data, only: Sigma_on_disk
  use citrans, only: citrans_csf2sd, citrans_sd2csf, citrans_sort
  use sguga, only: CIStruct, EXStruct, SGStruct
  use rasscf_global, only: DoFaro
  use faroald, only: gtuvx, htu, my_norb, sigma_update

  type(SGStruct), intent(in) :: SGS
  type(EXStruct), intent(in) :: EXS
  type(CIStruct), intent(in) :: CIS
  integer(kind=iwp), intent(in) :: STSYM, nCSF, ntemp, ndeta, ndetb, nTU, nTUVX
  real(kind=wp), intent(in) :: CI_Vec(nCSF), TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(out) :: Sigma_Vec(nCSF)
  real(kind=wp), intent(inout), target :: ctemp(ntemp), sigtemp(ntemp)
  integer(kind=iwp) :: it, itu, ituvx, iu, iv, ix, ixmax
  real(kind=wp), pointer :: Faroald_PSI(:,:), Faroald_SGM(:,:)

  if (DOFARO) then

    htu(:,:) = Zero
    gtuvx(:,:,:,:) = Zero
    itu = 0
    ituvx = 0
    do it=1,my_norb
      do iu=1,it
        itu = itu+1
        htu(iu,it) = TU(itu)
        htu(it,iu) = TU(itu)
        do iv=1,it
          ixmax = iv
          if (it == iv) ixmax = iu
          do ix=1,ixmax
            ituvx = ituvx+1
            GTUVX(IT,IU,IV,IX) = TUVX(ITUVX)
            GTUVX(IU,IT,IV,IX) = TUVX(ITUVX)
            GTUVX(IT,IU,IX,IV) = TUVX(ITUVX)
            GTUVX(IU,IT,IX,IV) = TUVX(ITUVX)
            GTUVX(IV,IX,IT,IU) = TUVX(ITUVX)
            GTUVX(IX,IV,IT,IU) = TUVX(ITUVX)
            GTUVX(IV,IX,IU,IT) = TUVX(ITUVX)
            GTUVX(IX,IV,IU,IT) = TUVX(ITUVX)
          end do
        end do
      end do
    end do

    Faroald_Psi(1:nDetA,1:nDetB) => ctemp(:)
    Faroald_SGM(1:nDetA,1:nDetB) => sigtemp(:)

    call SG_REORD(SGS,EXS,STSYM,0,CIS%nCSF(STSYM),CI_Vec,ctemp)
    call CITRANS_SORT('C',ctemp,Sigma_Vec)
    Faroald_PSI(:,:) = Zero
    call CITRANS_CSF2SD(Sigma_Vec,Faroald_PSI)
    Faroald_SGM(:,:) = Zero
    call SIGMA_UPDATE(HTU,GTUVX,Faroald_SGM,Faroald_PSI)
    call CITRANS_SD2CSF(Faroald_SGM,Sigma_Vec)
    call CITRANS_SORT('O',Sigma_Vec,ctemp)
    call SG_Reord(SGS,EXS,STSYM,1,CIS%nCSF(STSYM),ctemp,Sigma_Vec)

    nullify(Faroald_Psi)
    nullify(Faroald_SGM)

  else

    ! Convert the CI-vector from CSF to Det. basis
    ! sigtemp is scratch, converted vector is stored in ctemp.
    ! Note that ctemp is of the size of the nDet. basis

    ctemp(1:nCSF) = CI_Vec(1:nCSF)
    sigtemp(:) = Zero
    call csdtvc(ctemp,sigtemp,1,STSym,1)

    ! Calling Lucia to determine the sigma vector
    call Lucia_Util('Sigma', &
                    CI_Vector=ctemp(:), &
                    Sigma_Vector=sigtemp(:), &
                    nTU=size(TU),TU=TU, &
                    nTUVX=size(TUVX),TUVX=TUVX)

    ! Set mark so densi_master knows that the Sigma-vector exists on disk.
    Sigma_on_disk = .true.
    ! Convert the Sigma vector from Det. to CSF basis. Converted vector is
    ! stored in Sigma_vec.
    call CSDTVC(Sigma_Vec,sigtemp,2,stSym,1)

  end if

end subroutine Mk_H_Psi

!#define _SGUGA_VERIFY_
subroutine Mk_pdms(CIVec,nCIVEC,D,SD,P,PA,nD,nP)

  use Lucia_Interface, only: Lucia_Util
  use rasscf_global, only: DoFaro, NAC
  use sguga_states, only: CIS, EXS, SGS
  use general_data, only: STSYM
  use faroald, only: fold_two_pdm, ndeta, ndetb, one_pdm, two_pdm
  use citrans, only: citrans_csf2sd, citrans_sort
# ifdef _SGUGA_VERIFY_
  use gas_data, only: iDoGAS
  use rasscf_global, only: NACPAR, NACPR2
  use general_data, only: NCONF
  use Definitions, only: u6
# endif
  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp), intent(in) :: nCIVEC
  real(kind=wp), intent(inout) :: CIVEC(nCIVEC)
  real(kind=wp), intent(out), optional :: D(:), SD(:), P(:), PA(:)
  integer(kind=iwp), intent(in) :: nD, nP
  real(kind=wp), allocatable :: CIV(:), D_FAROALD(:,:), D_loc(:), Faroald_Psi(:,:), P_Faroald(:,:,:,:), P_loc(:), PA_loc(:), &
                                SD_FAROALD(:,:), SD_loc(:), temp(:)
  integer(kind=iwp), parameter :: iState = 1
# ifdef _SGUGA_VERIFY_
  real(kind=wp) :: Check_D1, Check_P, Check_PA
  real(kind=wp), allocatable :: D_Sguga(:), P_Sguga(:), PA_sguga(:)
# endif

  if (DoFaro) then
    call mma_allocate(D_loc,nD,Label='D_loc')
    call mma_allocate(SD_loc,nD,Label='SD_loc')
    call mma_allocate(P_loc,nP,Label='P_loc')
    call mma_allocate(PA_loc,nP,Label='PA_loc')

    call mma_allocate(CIV,nDetA*nDetB,Label='CIV')
    CIV(:) = Zero
    call mma_allocate(temp,nDetA*nDetB,Label='temp')
    call mma_allocate(Faroald_Psi,nDetA,nDetB,Label='Psi')

    call SG_Reord(SGS(istate),EXS(istate),STSYM,0,CIS(istate)%nCSF(STSYM),CIVEC,CIV)
    Temp(:) = Zero
    call CITRANS_SORT('C',CIV,temp)
    Faroald_Psi(:,:) = Zero
    call CITRANS_CSF2SD(temp,Faroald_PSI)

    call mma_deallocate(CIV)
    call mma_deallocate(temp)

    call mma_allocate(D_Faroald,NAC,NAC)
    call mma_allocate(SD_Faroald,NAC,NAC)
    call One_pdm(Faroald_Psi,D_Faroald,SD_Faroald)
    call Fold2(1,[NAC],D_faroald,D_loc)
    call Fold2(1,[NAC],SD_faroald,SD_loc)

    call mma_allocate(P_Faroald,NAC,NAC,NAC,NAC)
    call two_pdm(Faroald_psi,P_Faroald)
    call Fold_Two_pdm(P_Faroald,P_loc,PA_loc)

    call mma_deallocate(Faroald_Psi)
    call mma_deallocate(P_faroald)
    call mma_deallocate(D_faroald)
    call mma_deallocate(SD_faroald)

    if (present(D)) D(1:nD) = D_loc(1:nD)
    if (present(SD)) SD(1:nD) = SD_loc(1:nD)
    if (present(P)) P(1:nP) = P_loc(1:nP)
    if (present(PA)) PA(1:nP) = PA_loc(1:nP)

    call mma_deallocate(D_loc)
    call mma_deallocate(SD_loc)
    call mma_deallocate(P_loc)
    call mma_deallocate(PA_loc)
  else
    call Lucia_Util('Densi',CI_Vector=CIVEC)
  end if

 ! temporary code to verify the functionality of the SGUGA code and its interface
# ifdef _SGUGA_VERIFY_
  if (.not. iDoGAS) then

    call mma_allocate(CIV,nConf,Label='CIV')
    call SG_Reord(SGS(istate),EXS(istate),STSYM,0,CIS(istate)%nCSF(STSYM),CIVEC,CIV)

    ! Test the one-particle density matrix
    !call TriPrt('D(Lucia)',' ',D,NAC)
    Check_D1 = CheckSum(D,NACPAR)
    !write(u6,*) 'Check_D1=',Check_D1
    call mma_allocate(D_sguga,NAC*(NAC+1)/2)

    call sg_one_pdm(SGS(istate),CIS(istate),EXS(istate),CIV,size(CIV),STSYM,D_sguga,size(D_sguga))
    if (abs(CheckSum(D_sguga,NACPAR)-Check_D1)/size(D_sguga) > 1.0e12_wp) then
      !write(u6,*) 'Check_D1=',Check_D1
      Check_D1 = CheckSum(D_sguga,NACPAR)
      write(u6,*) 'SGUGA error in D1Mat'
      call Abend()
    end if
    call mma_deallocate(D_sguga)

    ! Test the one-particle spin-density matrix
    ! This option is not yet developed for the SGUGA code. To come...

    ! Test the symmetric two-particle density matrix.
    !call TRIPRT('P(Lucia)',' ',P,NACPAR)
    Check_P = CheckSum(P,NACPR2)
    !write(u6,*) 'Check_P=',Check_P
    !call TRIPRT('PA(Lucia)',' ',PA,NACPAR)
    Check_PA = CheckSum(PA,NACPR2)
    !write(u6,*) 'Check_PA=',Check_PA

    !call mma_allocate(P_sguga,NAC**4,Label='P')
    !call sg_two_pdm_full(SGS(istate),CIS(istate),EX(istate)S,CIV,SIZE(CIV),STSYM,P_sguga,NAC)
    !call mma_deallocate(P_sguga)

    call mma_allocate(P_sguga,NACPR2,Label='P')
    call mma_allocate(PA_sguga,NACPR2,Label='PA')

    call sg_two_pdm(SGS(istate),CIS(istate),EXS(istate),CIV,size(CIV),STSYM,P_sguga,PA_sguga,NACPAR*(NACPAR+1)/2)

    !call TRIPRT('P(SGUGA)',' ',P_sguga,NACPAR)
    if (abs(CheckSum(P_sguga,NACPR2)-Check_P)/size(p_sguga) > 1.0e-12_wp) then
      Check_P = CheckSum(P_sguga,NACPR2)
      !write(u6,*) 'Check_P=',Check_P
      write(u6,*) 'SGUGA error in P'
      call Abend()
    end if

    !call TRIPRT('PA(SGUGA)',' ',PA_sguga,NACPAR)
    if (abs(CheckSum(PA_sguga,NACPR2)-Check_PA)/size(p_sguga) > 1.0e-12_wp) then
      Check_PA = CheckSum(PA_sguga,NACPR2)
      !write(u6,*) 'Check_PA=',Check_PA
      write(u6,*) 'SGUGA error in PA'
      call Abend()
    end if

    call mma_deallocate(PA_sguga)
    call mma_deallocate(P_sguga)
    call mma_deallocate(CIV)

  end if
# endif
  ! end temporary code

end subroutine Mk_pdms

function Checksum(A,nA)

  real(kind=wp) :: Checksum
  integer(kind=iwp), intent(in) :: nA
  real(kind=wp), intent(in) :: A(nA)
  integer(kind=iwp) :: i

  Checksum = Zero
  do i=1,nA
    Checksum = Checksum+abs(A(i))/real(i,kind=wp)
  end do

end function Checksum

end module CI_Interfaces
