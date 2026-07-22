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
use definitions, only: wp, iwp

Private

Public :: Mk_H_Psi, Mk_pdms

contains

Subroutine Mk_H_Psi(SGS,EXS,CIS,STSYM,nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,ndeta,ndetb, &
                    nTU,TU,nTUVX,TUVX)
use Lucia_Interface, only: Lucia_Util
use lucia_data, only: Sigma_on_disk
use citrans, only: citrans_csf2sd, citrans_sd2csf, citrans_sort
use sguga, only: SGStruct, EXStruct, CIStruct
use rasscf_global, only: DoFaro
use Constants, only: Zero
use faroald, only: my_norb, sigma_update, htu, gtuvx
Implicit None

type(SGStruct), intent(in) :: SGS
type(EXStruct), intent(in) :: EXS
type(CIStruct), intent(in) :: CIS
integer(kind=iwp), intent(in):: STSYM, nCSF
real(kind=wp), intent(in) :: CI_Vec(nCSF)
real(kind=wp), intent(out) :: Sigma_Vec(nCSF)
integer(kind=iwp), intent(in) :: ntemp,ndeta,ndetb
real(kind=wp), intent(inout), target :: ctemp(ntemp), sigtemp(ntemp)
integer(kind=iwp), intent(in):: nTU, nTUVX
real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)

integer(kind=iwp) :: itu, ituvx, it, iu, iv, ixmax, ix
real(kind=wp), pointer:: Faroald_PSI(:,:), Faroald_SGM(:,:)

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

  Faroald_Psi => Null()
  Faroald_SGM => Null()

else

  ! Convert the CI-vector from CSF to Det. basis
  ! sigtemp is scratch, converted vector is stored in ctemp.
  ! Note that ctemp is of the size of the nDet. basis

  ctemp(1:nCSF) = CI_Vec(1:nCSF)
  sigtemp(:) = Zero
  call csdtvc(ctemp,sigtemp,1,STSym,1)

  ! Calling Lucia to determine the sigma vector
  call Lucia_Util('Sigma',                   &
                  CI_Vector=ctemp(:),        &
                  Sigma_Vector=sigtemp(:),   &
                  nTU=Size(TU),TU=TU,        &
                  nTUVX=Size(TUVX),TUVX=TUVX)

  ! Set mark so densi_master knows that the Sigma-vector exists on disk.
  Sigma_on_disk = .true.
  ! Convert the Sigma vector from Det. to CSF basis. Converted vector is
  ! stored in Sigma_vec.
  call CSDTVC(Sigma_Vec,sigtemp,2,stSym,1)

end if

End Subroutine Mk_H_Psi


!#define _SGUGA_VERIFY_
 Subroutine Mk_pdms(CIVec,nCIVEC,D,SD,P,PA,nD,nP)
 use Lucia_Interface, only: Lucia_Util
 use stdalloc, only: mma_allocate, mma_deallocate
 use rasscf_global, only: DoFaro, NAC
 use sguga_states, only: CIS, SGS, EXS
 use general_data, only: STSYM
 use faroald, only: ndeta, ndetb ,one_pdm, two_pdm, fold_two_pdm
 use citrans, only: citrans_csf2sd, citrans_sort
#ifdef _SGUGA_VERIFY_
 use gas_data, only: iDoGAS
 use rasscf_global, only: NACPAR, NACPR2
 use general_data, only: NCONF
#endif
 use constants, only: Zero

 implicit none
 integer(kind=iwp), intent(in) :: nCIVEC
 real(kind=wp), intent(inout) :: CIVEC(nCIVEC)
 real(kind=wp), intent(out), optional :: D(:), SD(:), P(:), PA(:)
 integer(kind=iwp), intent(in) :: nD, nP

 real(kind=wp), allocatable :: D_loc(:), SD_loc(:), P_loc(:), PA_loc(:)
 real(kind=wp), allocatable :: D_FAROALD(:,:)
 real(kind=wp), allocatable :: SD_FAROALD(:,:)
 real(kind=wp), allocatable :: Faroald_Psi(:,:)
 real(kind=wp), allocatable :: P_Faroald(:,:,:,:)
 real(kind=wp), allocatable :: CIV(:), temp(:)
 integer(kind=iwp), parameter :: iState=1

#ifdef _SGUGA_VERIFY_
real(kind=wp) :: Check_D1, Check_P, Check_PA
real(kind=wp), allocatable :: D_Sguga(:)
real(kind=wp), allocatable :: P_Sguga(:), PA_sguga(:)
#endif

 If (DoFaro) Then
   call mma_allocate(D_loc,nD,Label='D_loc')
   call mma_allocate(SD_loc,nD,Label='SD_loc')
   call mma_allocate(P_loc,nP,Label='P_loc')
   call mma_allocate(PA_loc,nP,Label='PA_loc')

   Call mma_allocate(CIV,nDetA*nDetB,Label='CIV')
   CIV(:)=Zero
   Call mma_allocate(temp,nDetA*nDetB,Label='temp')
   Call mma_allocate(Faroald_Psi,nDetA,nDetB,Label='Psi')

   call SG_Reord(SGS(istate),EXS(istate),STSYM,0,CIS(istate)%nCSF(STSYM),CIVEC,CIV)
   Temp(:)=Zero
   call CITRANS_SORT('C',CIV,temp)
   Faroald_Psi(:,:)=Zero
   call CITRANS_CSF2SD(temp,Faroald_PSI)

   Call mma_deallocate(CIV)
   Call mma_deallocate(temp)

   Call mma_allocate(D_Faroald,NAC,NAC)
   Call mma_allocate(SD_Faroald,NAC,NAC)
   Call One_pdm(Faroald_Psi,D_Faroald,SD_Faroald)
   Call Fold2(1,[NAC],D_faroald,D_loc)
   Call Fold2(1,[NAC],SD_faroald,SD_loc)

   Call mma_allocate(P_Faroald,NAC,NAC,NAC,NAC)
   Call two_pdm(Faroald_psi,P_Faroald)
   Call Fold_Two_pdm(P_Faroald,P_loc,PA_loc)

   Call mma_deallocate(Faroald_Psi)
   Call mma_deallocate(P_faroald)
   Call mma_deallocate(D_faroald)
   Call mma_deallocate(SD_faroald)

   If (Present(D)) D(1:nD)=D_loc(1:nD)
   If (Present(SD)) SD(1:nD)=SD_loc(1:nD)
   If (Present(P)) P(1:nP)=P_loc(1:nP)
   If (Present(PA)) PA(1:nP)=PA_loc(1:nP)

   call mma_deallocate(D_loc)
   call mma_deallocate(SD_loc)
   call mma_deallocate(P_loc)
   call mma_deallocate(PA_loc)
 Else
   call Lucia_Util('Densi',CI_Vector=CIVEC)
 End If

! temporary code to verify the functionality of the SGUGA code and its interface
#ifdef _SGUGA_VERIFY_
        If (.NOT.iDoGAS) Then

          Call mma_allocate(CIV,nConf,Label='CIV')
          call SG_Reord(SGS(istate),EXS(istate),STSYM,0,CIS(istate)%nCSF(STSYM),CIVEC,CIV)

!         Test the one-particle density matrix
!         Call TriPrt('D(Lucia)',' ',D,NAC)
          Check_D1=CheckSum(D,NACPAR)
!         Write (6,*) 'Check_D1=',Check_D1
          Call mma_allocate(D_sguga,NAC*(NAC+1)/2)

          call sg_one_pdm(SGS(istate),CIS(istate),EXS(istate),CIV,SIZE(CIV),STSYM,D_sguga,Size(D_sguga))
          If (ABS(CheckSum(D_sguga,NACPAR)-Check_D1)/SIZE(D_sguga)>1.0e12_wp) Then
!            Write (6,*) 'Check_D1=',Check_D1
             Check_D1=CheckSum(D_sguga,NACPAR)
             Write (6,*) 'SGUGA error in D1Mat'
             Call Abend()
          End If
          Call mma_deallocate(D_sguga)

!         Test the one-particle spin-density matrix
!         This option is not yet developed for the SGUGA code. To come...


!         Test the symmetric two-particle density matrix.
!         call TRIPRT('P(Lucia)',' ',P,NACPAR)
          Check_P=CheckSum(P,NACPR2)
!         Write (6,*) 'Check_P=',Check_P
!         call TRIPRT('PA(Lucia)',' ',PA,NACPAR)
          Check_PA=CheckSum(PA,NACPR2)
!         Write (6,*) 'Check_PA=',Check_PA

!         Call mma_allocate(P_sguga,NAC**4,Label='P')
!         Call sg_two_pdm_full(SGS(istate),CIS(istate),EX(istate)S,CIV,SIZE(CIV),STSYM,P_sguga,NAC)
!         Call mma_deallocate(P_sguga)

          Call mma_allocate(P_sguga,NACPR2,Label='P')
          Call mma_allocate(PA_sguga,NACPR2,Label='PA')

          Call sg_two_pdm(SGS(istate),CIS(istate),EXS(istate),CIV,SIZE(CIV),STSYM,P_sguga,PA_sguga,NACPAR*(NACPAR+1)/2)

!         call TRIPRT('P(SGUGA)',' ',P_sguga,NACPAR)
          If (ABS(CheckSum(P_sguga,NACPR2)-Check_P)/SIZE(p_sguga)>1.0e-12_wp) Then
             Check_P=CheckSum(P_sguga,NACPR2)
!            Write (6,*) 'Check_P=',Check_P
             Write (6,*) 'SGUGA error in P'
             Call Abend()
          End If

!         call TRIPRT('PA(SGUGA)',' ',PA_sguga,NACPAR)
          If (ABS(CheckSum(PA_sguga,NACPR2)-Check_PA)/SIZE(p_sguga)>1.0e-12_wp) Then
             Check_PA=CheckSum(PA_sguga,NACPR2)
!            Write (6,*) 'Check_PA=',Check_PA
             Write (6,*) 'SGUGA error in PA'
             Call Abend()
          End If

          Call mma_deallocate(PA_sguga)
          Call mma_deallocate(P_sguga)
          Call mma_deallocate(CIV)

        END IF
#endif
! end temporary code

 End Subroutine Mk_pdms

 Function Checksum(A,nA)
 real(kind=wp) :: Checksum
 integer(kind=iwp), intent(in):: nA
 real(kind=wp), intent(in):: A(nA)
 integer(kind=iwp) :: i
 Checksum=0.0_wp
 Do i = 1, nA
    Checksum = Checksum + Abs(A(i))/real(i,kind=wp)
 End Do
 End Function Checksum

End module CI_Interfaces
