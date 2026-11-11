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

#define _SANITY_CHECK_
#define _SGUGA_VERIFY_
module CI_Interfaces
use sguga, only: SG_Free, SG_ReOrd
use Lucia_Interface, only: Lucia_Util
use faroald, only: my_norb, sigma_update, htu, gtuvx, ndeta, ndetb ,transition_one_pdm, one_pdm, two_pdm, fold_two_pdm, &
                   traci_faroald
use citrans, only: citrans_csf2sd, citrans_sd2csf, citrans_sort
use rasscf_global, only: DoFaro
use general_data, only: ISPIN, STSYM, nDet, nConf, nAsh, nActEl, NAC=>NLEV
use spinfo, only: NDTASM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use definitions, only: wp, iwp, u6
#ifdef _SGUGA_VERIFY_
use sguga, only: CIS, SGS, EXS
use general_data, only: iDoGAS, nRsPrt
#endif

#include "intent.fh"

Private

Public :: Mk_Tra_Psi, Mk_H_Psi, Mk_H_Psi_Det, Mk_PDMS, Mk_T1DM, Mk_CI_Diag, CI_Timer, CI_Initialize, CI_Close

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

Subroutine Mk_H_Psi(iState,STSYM,nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,ndeta,ndetb,nTU,TU,nTUVX,TUVX)
use lucia_data, only: Sigma_on_disk, IREFSM
use cands, only: ICSM, ISSM
Implicit None

integer(kind=iwp), intent(in):: iState, STSYM, nCSF
real(kind=wp), intent(in) :: CI_Vec(nCSF)
real(kind=wp), intent(out) :: Sigma_Vec(nCSF)
integer(kind=iwp), intent(in) :: ntemp,ndeta,ndetb
real(kind=wp), intent(inout), target :: ctemp(ntemp), sigtemp(ntemp)
integer(kind=iwp), intent(in):: nTU, nTUVX
real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)

integer(kind=iwp) :: itu, ituvx, it, iu, iv, ixmax, ix
real(kind=wp), pointer:: Faroald_PSI(:,:), Faroald_SGM(:,:)
integer(kind=iwp) :: IREFSM_SAVE, nDet

#ifdef _SGUGA_VERIFY_
real(kind=wp) :: Check_Href, Check_H
real(kind=wp), allocatable ::  SG_PSI(:), SG_SGM(:)
#endif

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

  call SG_REORD(iState,STSYM,0,nCSF,CI_Vec,ctemp)
  call CITRANS_SORT('C',ctemp,Sigma_Vec)
  Faroald_PSI(:,:) = Zero
  call CITRANS_CSF2SD(Sigma_Vec,Faroald_PSI)
  Faroald_SGM(:,:) = Zero
  call SIGMA_UPDATE(HTU,GTUVX,Faroald_SGM,Faroald_PSI)
  call CITRANS_SD2CSF(Faroald_SGM,Sigma_Vec)
  call CITRANS_SORT('O',Sigma_Vec,ctemp)
  call SG_Reord(iState,STSYM,1,nCSF,ctemp,Sigma_Vec)

  Faroald_Psi => Null()
  Faroald_SGM => Null()

else

  IREFSM_SAVE=IREFSM
  IREFSM=STSYM
  ICSM=STSYM
  ISSM=STSYM

  nDET = NDTASM(IREFSM)

  ! Convert the CI-vector from CSF to Det. basis
  ! sigtemp is scratch, converted vector is stored in ctemp.
  ! Note that ctemp is of the size of the nDet. basis

  ctemp(1:nCSF) = CI_Vec(1:nCSF)
  sigtemp(:) = Zero
  call csdtvc(ctemp,sigtemp,1,STSym,0)
  ctemp(1:nDet)=sigtemp(1:nDet)

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
  call CSDTVC(Sigma_Vec,sigtemp,2,STSym,0)
  sigtemp(1:nCSF)=Sigma_Vec(1:nCSF)

  IREFSM=IREFSM_SAVE
  ICSM=IREFSM
  ISSM=IREFSM

end if

#ifdef _SGUGA_VERIFY_
If (.NOT.iDoGAS .and. nRsPrt==1) Then
   Call mma_allocate(SG_PSI,nCSF,Label='SG_PSI')
   call SG_Reord(iState,STSYM,0,nCSF,CI_VEC,SG_PSI)
   Check_Href=CheckSum(Sigma_Vec,nCSF)
   Call mma_allocate(SG_SGM,nCSF,Label='SG_SGM')
   Call sg_h_psi(SGS(iState),CIS(iState),EXS(iState),SG_PSI,nCSF,STSYM,SG_SGM,TUVX,nTUVX,TU,nTU)
   call SG_Reord(iState,STSYM,1,nCSF,SG_SGM,SG_PSI)
   Check_H   =CheckSum(SG_PSI,nCSF)
   If (Abs(Check_HRef-Check_H)/nCSF>1.0E-12_wp) Then
      Write (u6,*) 'SGUGA error in H|Psi>'
      Write (u6,*) 'Check_HRef=',Check_HRef
      Write (u6,*) 'Check_H   =',Check_H
      Call RecPrt('Sigma_Vec',' ',Sigma_Vec,1,nCSF)
      Call RecPrt('SG_SGM',' ',SG_PSI,1,nCSF)
      Call Abend()
   End If
   Sigma_Vec(1:nCSF) = SG_PSI(1:nCSF)
   Call mma_deallocate(SG_SGM)
   Call mma_deallocate(SG_PSI)
End if
#endif

End Subroutine Mk_H_Psi

Subroutine Mk_H_Psi_Det(istate,STSYM,nDet,CI_Vec,Sigma_Vec,nTU,TU,nTUVX,TUVX)
use lucia_data, only: Sigma_on_disk, IREFSM
use sd_csf_codec_module, only: sd_csf_codec, sd_to_csf, csf_to_sd
use sguga, only: EXS
use cands, only: ICSM, ISSM
Implicit None

integer(kind=iwp), intent(in):: istate, STSYM, nDet
real(kind=wp), intent(inout),target :: CI_Vec(nDet)
real(kind=wp), intent(inout),target :: Sigma_Vec(nDet)
integer(kind=iwp), intent(in):: nTU, nTUVX
real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)

real(kind=wp), allocatable, target:: temp1(:), temp2(:)  ! Local scratch
real(kind=wp), pointer:: Faroald_PSI(:,:), Faroald_SGM(:,:)

integer(kind=iwp) :: itu, ituvx, it, iu, iv, ixmax, ix, nCSF
integer(kind=iwp) :: iDet
integer(kind=iwp) :: IREFSM_SAVE

#ifdef _SGUGA_VERIFY_
!integer(kind=iwp) :: iDet
#ifdef _Compute_Full_CI_Matrix
real(kind=wp), allocatable :: Full(:,:)
#endif
#endif

!Note that we Lucia code doesn't explicily spin adapt the SD. Rather it just process the input coefficients.
! If we project the starting coefficients to the spin-adapted iCSF basis and then expand it back to the SD basis
! we have effectivly generated a SD coefficient vector that is implicitly spin-adapted. Since the CASVB code generates
! starting CI vectors that are not spin-adapted there will be a difference in the convergence patterns.

#ifdef _SGUGA_VERIFY_
#ifdef _Compute_Full_CI_Matrix
! Reference computation of the full CI matrix.
! in CSF basis
nCSF=CIS(istate)%nCSF(STSYM)
Call mma_allocate(Full,nCSF,nCSF,Label='Full')
Call sg_ci_full(SGS(istate),CIS(istate),EXS(istate),nCSF,STSym, &
                TUVX,nTUVX, &
                TU,nTU, &
                Full)
Call mma_deallocate(Full)

Call mma_allocate(Full,nDet,nDet,Label='Full')
Call ci_full_det(nDet,STSym, &
                TUVX,nTUVX, &
                TU,nTU, &
                Full)
Call mma_deallocate(Full)
! in SD basis
#endif
#endif

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

  ! We need to convert this to the correct form! The input is in LUCIA format.

  nCSF=CIS(istate)%nCSF(STSym) ! ToDo: This depends on the SGUGA part of the code.
  if (nCSF>nDet) Then
     Write(u6,*) 'Mk_H_Psi_Det: nCSF>nDet'
     Call Abend()
  End If

  Call mma_allocate(temp1,nDet,Label='temp1')
  Call mma_allocate(temp2,nDet,Label='temp2')
  Faroald_Psi(1:nDetA,1:nDetB) => temp2(:)
  Faroald_SGM(1:nDetA,1:nDetB) => temp1(:)

  temp2(:)=CI_Vec(:)
  ! Convert from LUCIA SD basis to Lucia CSF basis
  temp1(:)=Zero
  Call SD_CSF_Codec(temp2,temp1,SD_to_CSF,STSym)
  ! Now use the FAROALD machinary to go to FAROALD SD basis
  call SG_REORD(iState,STSYM,0,nCSF,temp1,temp2)
  call CITRANS_SORT('C',temp2,temp1)
  Faroald_PSI(:,:) = Zero
  call CITRANS_CSF2SD(temp1,Faroald_PSI)


  call SIGMA_UPDATE(HTU,GTUVX,Faroald_SGM,Faroald_PSI)

  call CITRANS_SD2CSF(Faroald_SGM,Sigma_Vec)
  call CITRANS_SORT('O',Sigma_Vec,temp2)
  call SG_Reord(iState,STSYM,1,nCSF,temp2,temp1)

  ! Convert from LUCIA CSF basis to LUCIA SD basis
  Call SD_CSF_Codec(Sigma_Vec,temp1,CSF_to_SD,STSym)

  Faroald_Psi => Null()
  Faroald_SGM => Null()
  Call mma_deallocate(temp2)
  Call mma_deallocate(temp1)

else

  IREFSM_SAVE=IREFSM
  IREFSM=STSYM
  ICSM=STSYM
  ISSM=STSYM

  ! Convert the CI-vector from CSF to Det. basis
  ! sigtemp is scratch, converted vector is stored in ctemp.
  ! Note that ctemp is of the size of the nDet. basis

  ! Calling Lucia to determine the sigma vector
  call Lucia_Util('Sigma',                   &
                  CI_Vector=CI_Vec(:),        &
                  Sigma_Vector=Sigma_Vec(:),   &
                  nTU=Size(TU),TU=TU,        &
                  nTUVX=Size(TUVX),TUVX=TUVX)

  ! Set mark so densi_master knows that the Sigma-vector exists on disk.
  Sigma_on_disk = .true.
  ! Convert the Sigma vector from Det. to CSF basis. Converted vector is
  ! stored in Sigma_vec.

  IREFSM=IREFSM_SAVE
  ICSM=IREFSM
  ISSM=IREFSM
endif

#ifdef _SGUGA_VERIFY_
  If (SGS(iState)%nRsPrt /= 1) Then
     Write (6,*) ' This option is only available for CASSCF'
     Call abend()
  End If
  If (SGS(iState)%nSym /= 1) Return
  nCSF=CIS(istate)%nCSF(STSym)
  Call mma_allocate(temp1,nDet,Label='temp1')
  Call mma_allocate(temp2,nDet,Label='temp2')

  temp2(:)=CI_Vec(:)
  ! Convert from LUCIA SD basis to Lucia CSF basis
  temp1(:)=Zero
  Call SD_CSF_Codec(temp2,temp1,SD_to_CSF,STSYM)
  ! Convert from Lucia CSF basis to SGUGA CSF basis
  !     iMode=0 --> from SGA to split graph GUGA order                   *
  !     iMode=1 --> from split graph GUGA to SGA order                   *
  call SG_Reord(iState,STSYM,0,nCSF,temp1,temp2)

  Call sg_h_psi(SGS(iState),CIS(iState),EXS(iState),temp2,nCSF,STSYM,temp1,TUVX,nTUVX,TU,nTU)

  ! Convert from SGUGA CSF basis to Lucia CSF basis
  call SG_Reord(iState,STSYM,1,nCSF,temp1,temp2)
  ! Convert from Lucia CSF basis to Lucia SD Basis
   call SD_CSF_Codec(temp1,temp2,CSF_to_SD,STSym)

  Do iDet = 1, nDet
     If (Abs(temp1(iDet)-Sigma_vec(iDet))>1.0e-10_wp) Then
        Write (u6,*) 'Mk_H_Psi_Det: Abs(temp1(iDet)-Sigma_vec(iDet))>1.0e-10_wp)'
        Write (u6,*) ' iDet=',iDet
        Write (u6,*) 'temp1(iDet)=    ',temp1(iDet)
        Write (u6,*) 'Sigma_vec(iDet)=', Sigma_vec(iDet)
        Call Abend()
     End If
  End Do

  Sigma_Vec(1:nDet)=temp1(1:nDet)
  Call mma_deallocate(temp1)
  Call mma_deallocate(temp2)
#endif

#ifdef _SGUGA_VERIFY_

#ifdef _Compute_Full_CI_Matrix
Contains
subroutine ci_full_det(nDet,STSym, &
                       TUVX,nTUVX, &
                       TU,nTU, &
                       Full)
use constants, only: Zero, One
use lucia_data, only: IREFSM
use cands, only: ICSM, ISSM

  integer(kind=iwp), intent(in):: nDet,STSym,nTUVX,nTU
  real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)
  real(kind=wp), intent(inout):: Full(nDet,nDet)

  integer(kind=iwp) IREFSM_SAVE, iDet, jDet

  real(kind=wp), allocatable :: CI_VEC(:), Sigma_VEC(:)

  Call mma_allocate(CI_Vec,nDet,Label='CI_Vec')
  Call mma_allocate(Sigma_Vec,nDet,Label='Sigma_Vec')

  IREFSM_SAVE=IREFSM
  IREFSM=STSYM
  ICSM=STSYM
  ISSM=STSYM

  Do iDet = 1, nDet
     CI_Vec(:) = Zero
     CI_Vec(iDet) = One

  ! Calling Lucia to determine the sigma vector
     Sigma_Vec(:) = Zero
     call Lucia_Util('Sigma',                   &
                     CI_Vector=CI_Vec(:),        &
                     Sigma_Vector=Sigma_Vec(:),   &
                     nTU=nTU,TU=TU,        &
                     nTUVX=nTUVX,TUVX=TUVX)

     Do jDet = 1, nDet
        Full(iDet,jDet) = Sigma_Vec(jDet)
     End Do

  End Do
  IREFSM=IREFSM_SAVE
  ICSM=IREFSM
  ISSM=IREFSM

  Call mma_deallocate(CI_Vec)
  Call mma_deallocate(Sigma_Vec)

  Call RecPrt('CI Matrix in SD basis',' ',Full,ndet,ndet)

end subroutine ci_full_det
#endif
#endif

End Subroutine Mk_H_Psi_Det

!***********************************************************************************************************************************
!***********************************************************************************************************************************

Subroutine Mk_T1DM(Bra_Vec, Ket_Vec, nVec, T1DM, nT1DM, T1SDM)
use Lucia_Data, only: DTmp, DSTmp
#ifdef _SGUGA_VERIFY_
use general_data, only: iDoGAS
#endif
implicit none
integer(kind=iwp), intent(in) :: nVec, nT1DM
real(kind=wp), intent(_IN_):: Bra_Vec(nVec), Ket_Vec(nVec)
real(kind=wp), intent(out):: T1DM(nT1DM)
real(kind=wp), optional, intent(out):: T1SDM(nT1DM)

real(kind=wp), allocatable:: T1SDM_loc(:), CIV(:), Bra_SD(:,:), Ket_SD(:,:), temp(:)
integer(kind=iwp), parameter :: iState=1
logical(kind=iwp) :: free_D, free_DS

If (DoFaro) Then
   Call mma_allocate(T1SDM_loc,nT1DM,Label='T1SDM')

   Call mma_allocate(CIV,nDetA*nDetB,Label='CIV')
   Call mma_allocate(temp,nDetA*nDetB,Label='temp')
   Call mma_allocate(Bra_SD,nDetA,nDetB,Label='Psi')
   Call mma_allocate(Ket_SD,nDetA,nDetB,Label='Psi')

   CIV(:)=Zero
   call SG_Reord(iState,STSYM,0,nVec,Bra_Vec,CIV)
   Temp(:)=Zero
   call CITRANS_SORT('C',CIV,temp)
   Bra_SD(:,:)=Zero
   call CITRANS_CSF2SD(temp,Bra_SD)

   CIV(:)=Zero
   call SG_Reord(iState,STSYM,0,nVec,Ket_Vec,CIV)
   Temp(:)=Zero
   call CITRANS_SORT('C',CIV,temp)
   Ket_SD(:,:)=Zero
   call CITRANS_CSF2SD(temp,Ket_SD)

!  Don't ask me why we need to do it the other way?
!  Call Transition_One_PDM(Bra_SD,Ket_SD,T1DM,T1SDM_loc)
   Call Transition_One_PDM(Ket_SD,Bra_SD,T1DM,T1SDM_loc)

   If (present(T1SDM)) T1SDM(1:nT1DM) = T1SDM_loc(1:nT1DM)

   Call mma_deallocate(Ket_SD)
   Call mma_deallocate(Bra_SD)
   Call mma_deallocate(temp)
   Call mma_deallocate(CIV)
   Call mma_deallocate(T1SDM_loc)
Else
   Free_D =.False.
   Free_DS=.False.
   If (.Not.Allocated(DTMP)) Then
      Free_D =.True.
      Call mma_allocate(DTMP,nT1DM,Label='DTMP')
   End If
   If (.Not.Allocated(DSTMP)) Then
      Free_DS=.True.
      Call mma_allocate(DSTMP,nT1DM,Label='DSTMP')
   End If

   call Lucia_Util('Densi',CI_Vector=Ket_Vec(:),RVec=Bra_Vec(:))
   T1DM(1:nT1DM)=DTMP(1:nT1DM)
   If (present(T1SDM)) T1SDM(1:nT1DM) = DSTMP(1:nT1DM)

   If (Free_D ) Call mma_deallocate(DTMP )
   If (Free_DS) Call mma_deallocate(DSTMP)
End If
#ifdef _SGUGA_VERIFY_
If (.NOT.iDoGAS) Then
End If
#endif

End Subroutine Mk_T1DM

!***********************************************************************************************************************************
!***********************************************************************************************************************************

Subroutine Mk_pdms(CIVec,nCIVEC,D,SD,P,PA,nD,nP)
use lucia_data, only: DStmp, Dtmp, PAtmp, PTmp
#ifdef _SGUGA_VERIFY_
use rasscf_global, only: NACPAR, NACPR2
#endif

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
   Call mma_allocate(Temp,nDetA*nDetB,Label='temp')
   Call mma_allocate(Faroald_Psi,nDetA,nDetB,Label='Psi')

   call SG_Reord(istate,STSYM,0,nCIVEC,CIVEC,CIV)
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
   Call mma_allocate(DTMP,nD,Label='DTMP')
   Call mma_allocate(DSTMP,nD,Label='DSTMP')
   Call mma_allocate(PTMP,nP,Label='PTMP')
   Call mma_allocate(PATMP,nP,Label='PATMP')
   call Lucia_Util('Densi',CI_Vector=CIVEC)
   If (Present(D)) D(1:nD)=DTMP(1:nD)
   If (Present(SD)) SD(1:nD)=DSTMP(1:nD)
   If (Present(P)) P(1:nP)=PTMP(1:nP)
   If (Present(PA)) PA(1:nP)=PATMP(1:nP)
   Call mma_deallocate(DTMP)
   Call mma_deallocate(DSTMP)
   Call mma_deallocate(PTMP)
   Call mma_deallocate(PATMP)
 End If

! temporary code to verify the functionality of the SGUGA code and its interface
#ifdef _SGUGA_VERIFY_
        If (.NOT.iDoGAS) Then

          Call mma_allocate(CIV,nCIVEC,Label='CIV')
          call SG_Reord(iState,STSYM,0,nCIVEC,CIVEC,CIV)

!         Test the one-particle density matrix
          Check_D1=CheckSum(D,NACPAR)
          Call mma_allocate(D_sguga,NAC*(NAC+1)/2)

          call sg_one_pdm(SGS(istate),CIS(istate),EXS(istate),CIV,SIZE(CIV),STSYM,D_sguga,Size(D_sguga))
          If (ABS(CheckSum(D_sguga,NACPAR)-Check_D1)/SIZE(D_sguga)>1.0e-12_wp) Then
             Check_D1=CheckSum(D_sguga,NACPAR)
             Write (u6,*) 'SGUGA error in D1Mat'
             Call Abend()
          End If
          Call mma_deallocate(D_sguga)

!         Test the one-particle spin-density matrix
!         This option is not yet developed for the SGUGA code. To come...

!         Test the symmetric two-particle density matrix.
          Check_P=CheckSum(P,NACPR2)
          Check_PA=CheckSum(PA,NACPR2)

!         Call mma_allocate(P_sguga,NAC**4,Label='P')
!         Call sg_two_pdm_full(SGS(istate),CIS(istate),EX(istate)S,CIV,SIZE(CIV),STSYM,P_sguga,NAC)
!         Call mma_deallocate(P_sguga)

          Call mma_allocate(P_sguga,NACPR2,Label='P')
          Call mma_allocate(PA_sguga,NACPR2,Label='PA')

          Call sg_two_pdm(SGS(istate),CIS(istate),EXS(istate),CIV,SIZE(CIV),STSYM,P_sguga,PA_sguga,NACPAR*(NACPAR+1)/2)

          If (ABS(CheckSum(P_sguga,NACPR2)-Check_P)/SIZE(p_sguga)>1.0e-12_wp) Then
             Check_P=CheckSum(P_sguga,NACPR2)
             Write (u6,*) 'SGUGA error in P'
             Call Abend()
          End If

          If (ABS(CheckSum(PA_sguga,NACPR2)-Check_PA)/SIZE(p_sguga)>1.0e-12_wp) Then
             Check_PA=CheckSum(PA_sguga,NACPR2)
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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine CI_Timer(Opt)
use lucia_data, only: TDENSI, TSIGMA

character(Len=*) Opt

Select case (Opt(1:3))
Case ('Ini')
   !SVC: lucia timers
   tsigma(:) = Zero
   tdensi(:) = Zero
Case ('Pri')
  If (DoFaro) Then
  Else
    call PrtTim()
  End If
Case Default
  Write (u6,*) 'CI_Timer called with illegal option:',Opt
  Call Abend()
End Select

end subroutine CI_Timer

!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine CI_Initialize(DO_SGUGA,DO_FAROALD)
use sguga, only: SG_Test_EpqErs_Ref
use general_data, only: STSym
use Molcas, only: MxSym
use spinfo, only: ncsasm, ndtasm

implicit none
logical(kind=iwp), intent(in) :: DO_SGUGA,DO_FAROALD
integer(kind=iwp), parameter :: iState=1
character(len=8) :: Env

! Construct the determinant tables
! Initialize LUCIA and determinant control
call Lucia_Util('Ini')

! to get number of CSFs for GAS
! and number of determinants to store
nconf = MaxVal(ncsasm(1:mxsym))
nDet = MaxVal(ndtasm(1:mxsym))

! Initiate the SGUGA environment conditional to all flags
If (Do_SGUGA) Then
#ifdef _SGUGA_VERIFY_
    write(u6,'(1X,A)') '**EXPERIMENTAL**'
    write(u6,'(1X,A)') 'CI backend is SGUGA instead of LUCIA.'
    write(u6,'(1X,A)') '**EXPERIMENTAL**'
#endif
   call SG_Setup_RASSCF()
   ! Exhaustive verification of the SGUGA coupling-coefficient codes,
   ! only on request (used by the test/sguga group):
   !   MOLCAS_SGUGA_VERIFY=YES  <I|Epq|J> and <I|EpqErs|J>
   !   MOLCAS_SGUGA_VERIFY=EPQ  <I|Epq|J> only
   Env=' '
   call get_environment_variable('MOLCAS_SGUGA_VERIFY',Env)
   call UpCase(Env)
   select case (Env)
     case ('YES')
       call SG_Test_EpqErs_Ref(iState,STSym)
     case ('EPQ')
       call SG_Test_EpqErs_Ref(iState,STSym,Epq_Only=.true.)
   end select
End If

! faroald initializations
if (DO_FAROALD) then
  write(u6,'(1X,A)') '**EXPERIMENTAL**'
  write(u6,'(1X,A)') 'CI backend is FAROALD instead of LUCIA.'
  write(u6,'(1X,A)') '**EXPERIMENTAL**'
  call FAROALD_INIT(NACTEL,NASH(1),ISPIN)
  call CITRANS_INIT(NACTEL,NASH(1),ISPIN)
end if

End subroutine CI_Initialize

!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine CI_Close(iState)
integer(kind=iwp), Intent(in):: iState

! All of these are safe to call even of their corresponding environment have not been activated.

Call Lucia_Util('Close')

call faroald_free()
call citrans_free()

Call SG_free(iState)

end subroutine CI_Close

!***********************************************************************************************************************************
!***********************************************************************************************************************************

! Activate verification of the SGUGA to computed CI diagonal in the CSF basis.
! The analytic result is computed by sg_hdiag_fast (using sg_EpqEqp_diag).
! sg_hdiag_ref computes the diagonal brute force by setting the CI element iCSF to 1 and compute the matrix element
! <Psi|H|Psi> using the routine sg_h_psi.
!#define _SGUGA_VERIFY_EXTRA_

subroutine Mk_CI_Diag(CSFDIA,nConf,TU,nTU,TUVX,nTUVX)
use Lucia_Interface, only: Lucia_Util
use Lucia_data, only: SDREO
use spinfo, only: NCNFTP, NCSFTP, NDTFTP, NTYP
use faroald, only: hdiag, ncomb_tot
!use faroald, only: nPat, hdiag, ncomb_tot
#ifdef _SGUGA_VERIFY_EXTRA_
use sguga, only: SGS, EXS, CIS
#endif

integer(kind=iwp), intent(in) :: nConf, nTU, nTUVX
real(kind=wp), intent(out):: CSFDIA(nConf)
real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)

real(kind=wp), allocatable :: DDIA(:), DIAG(:), CSFDIA_Faroald(:)
integer(kind=iwp):: IPRINT=0

integer(kind=iwp) :: i
#ifdef _SGUGA_VERIFY_EXTRA_
integer(kind=iwp), parameter :: iState=1
real(kind=wp), allocatable :: TMP(:), CSFDIA_SGUGA(:), CSFDIA_FAST(:)
real(kind=wp) :: Err
#endif

! COMPUTE CI DIAGONAL IN DETERMINANT BASIS

call Lucia_Util('Diag',nTU=nTU,TU=TU,nTUVX=nTUVX,TUVX=TUVX)

call mma_allocate(DDIA,NDET,label='DETDIA')

call get_diag(DDIA,ndet)

! Todo: ndet comes from general_data!! Check this!
! TRANSFORM CI DIAGONAL FROM DET TO CSF BASIS
call CSDIAG(NCONF,ndet,CSFDIA,DDIA,NCNFTP(1,STSYM),NTYP,SDREO,NDTFTP,NCSFTP,IPRINT)

If (DoFaro) Then
! If (nDet/=nComb_tot) Then
!    Write (u6,*) 'Warning: nDet/=nComb'
!    Write (u6,*) 'iSpin=',iSpin
!    Write (u6,*) 'nConf=',nConf
!    Write (u6,*) 'nDet=',nDet
!    Write (u6,*) 'nPat=',nPat
!    Write (u6,*) 'nComb_tot',nComb_tot
!    Call Abend()
! End If
!Call RecPrt('DDIA(LUCIA)',' ',DDIA,1,nDet)
!Call RecPrt('CSFDIA(Lucia)',' ',CSFDIA,1,nConf)

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

  Call mma_allocate(Diag,nComb_tot,Label='Diag')
  Call mma_allocate(CSFDIA_Faroald,nConf,Label='CSFDIA_F')
  Diag(:)=Zero
  Call hdiag(htu,gtuvx,diag,CSFDIA_Faroald,nConf)
! Call RecPrt('Diag(FAROALD)',' ',Diag,1,nComb_tot)
  Do i=1,nComb_tot
     If (abs(ddia(i)-diag(i))/real(nComb_tot,wp)>1.0e-10_wp) Then
        Write (u6,*) 'i, abs(ddia(i)-diag(i))'
        Write (u6,*) i, abs(ddia(i)-diag(i))
        Write (u6,*) ddia(i),diag(i)
        Call RecPrt('DDIA(LUCIA)',' ',DDIA,1,nDet)
        Call RecPrt('Diag(FAROALD)',' ',Diag,1,nComb_tot)
        Call abend()
     End If
  End Do
  Do i=1,nConf
     If (abs(CSFDIA(i)-CSFDIA_Faroald(i))/real(nConf,wp)>1.0e-10_wp) Then
        Write (u6,*) 'i, abs(CSFDIA(i)-CSFDIA_Faroald(i))', i, abs(CSFDIA(i)-CSFDIA_Faroald(i))
        Write (u6,*) CSFDIA(i),CSFDIA_Faroald(i)
        Call RecPrt('CSFDIA(LUCIA)',' ',CSFDIA,1,nConf)
        Call RecPrt('CSFDIA(Faroald)',' ',CSFDIA_Faroald,1,nConf)
        Call abend()
     End If
  End Do
! Call RecPrt('CSFDIA(Faroald)',' ',CSFDIA_Faroald,1,nConf)
  CSFDIA(:)= CSFDIA_Faroald(:)

  Call mma_deallocate(CSFDIA_Faroald)
  Call mma_deallocate(Diag)

End If

! DEALLOCATE LOCAL MEMORY
call mma_deallocate(DDIA)

#ifdef _SGUGA_VERIFY_EXTRA_
If (.NOT.iDoGAS .and. nRsPrt==1) Then
Call mma_allocate(TMP,nConf,Label='CSFDIA_S')
Call sg_hdiag_ref(SGS(iState),CIS(iState),EXS(iState), &
                  nConf,STSYM, &
                  TUVX,nTUVX, &
                  TU,nTU, &
                  TMP)
!Call RecPrt('TMP',' ',TMP,1,nConf)
Call mma_allocate(CSFDIA_SGUGA,nConf,Label='CSFDIA_S')
Call SG_ReOrd(iState,STSYM,1,nConf,TMP,CSFDIA_SGUGA)
Call mma_deallocate(TMP)
! Note the SG_ReOrd introduce a sign change of some elements which we have to remove.
CSFDIA_SGUGA(:)=-ABS(CSFDIA_SGUGA(:))
Call RecPrt('CSFDIA(SGUGA_ref)',' ',CSFDIA_SGUGA,1,nConf)

Call mma_allocate(TMP,nConf,Label='CSFDIA_S')
Call sg_hdiag_fast(SGS(iState),CIS(iState),EXS(iState), &
                   nConf,STSYM, &
                   TUVX,nTUVX, &
                   TU,nTU, &
                   TMP)
!Call RecPrt('TMP',' ',TMP,1,nConf)
Call mma_allocate(CSFDIA_FAST,nConf,Label='CSFDIA_S')
Call SG_ReOrd(iState,STSYM,1,nConf,TMP,CSFDIA_FAST)
Call mma_deallocate(TMP)
! Note the SG_ReOrd introduce a sign change of some elements which we have to remove.
CSFDIA_FAST(:)=-ABS(CSFDIA_FAST(:))
Call RecPrt('CSFDIA(SGUGA_fast)',' ',CSFDIA_FAST,1,nConf)

Do i=1,nConf
   Err = Abs(CSFDIA_SGUGA(i)-CSFDIA_fast(i))
   If (Err > 1.0D-12) Then
      Write(6,*) 'CSF=',i
      Write(6,*) 'REF =',CSFDIA_SGUGA(i)
      Write(6,*) 'FAST=',CSFDIA_fast(i)
      Write(6,*) 'ERR =',Err
   End If
End Do

Call mma_deallocate(CSFDIA_SGUGA)
Call mma_deallocate(CSFDIA_FAST)
endif

#endif

End subroutine Mk_CI_Diag

subroutine Mk_Tra_Psi(iState,TU,nTU,CI_Vector,nCSF)

use sguga, only: sg_traci

integer(kind=iwp), intent(in) ::  iState,nTU, nCSF
real(kind=wp), intent(in) :: TU(nTU)
real(kind=wp), intent(inout) :: CI_Vector(nCSF)

real(kind=wp), allocatable :: Faroald_Psi(:,:), ctemp(:)
real(kind=wp), allocatable :: dtemp(:)

integer(kind=iwp) :: iCSF
real(kind=wp), allocatable :: CI_Faroald(:)
#ifdef _SGUGA_VERIFY_
real(kind=wp), allocatable :: CI_SGUGA(:)
#endif

#ifdef _SGUGA_VERIFY_
  Call mma_allocate(CI_SGUGA,nCSF,Label='CI_SGUGA')
  CI_SGUGA(:) = CI_Vector(:)
#endif

If (DoFaro) Then
  Call mma_allocate(CI_Faroald,nCSF,Label='CI_Faroald')
  CI_Faroald(:) = CI_Vector(:)

  Call mma_allocate(Faroald_Psi,nDetA,nDetB,Label='Faroald_Psi')
  Call mma_allocate(ctemp,Max(nCSF,nDetA*nDetB),Label='ctemp')
  Call mma_allocate(dtemp,nDetA*nDetB,Label='dtemp')

  call SG_REORD(iState,STSYM,0,nCSF,CI_Vector,ctemp)
  call CITRANS_SORT('C',ctemp,dTemp)
  Faroald_PSI(:,:) = Zero
  call CITRANS_CSF2SD(dTemp,Faroald_PSI)

  Call traci_faroald(TU,nTU,Faroald_Psi)

  call CITRANS_SD2CSF(Faroald_Psi,dTemp)
  call CITRANS_SORT('O',dTemp,ctemp)
! call SG_Reord(iState,STSYM,1,nCSF,cTemp,CI_Vector)
  call SG_Reord(iState,STSYM,1,nCSF,cTemp,CI_Faroald)

  Call mma_deallocate(dTemp)
  Call mma_deallocate(ctemp)
  Call mma_deallocate(Faroald_Psi)
End If
!Else

!  Note that the CI vector is transformed to Lucia format inside the routine! This
!  in difference to the other calls to Lucia_util.
   Call Lucia_util('TRACI', &
                   CI_Vector=CI_Vector(:), &
                   TU=TU(:), nTU=nTU)
!End If

if (DoFaro) Then
Do iCSF = 1, nCSF
   If (Abs(CI_Vector(iCSF)-CI_Faroald(iCSF))/nCSF>1.0E-12_wp) Then
      Write (6,*) 'Mk_Tra_Psi error, iCSF=',iCSF
      Write (6,*) 'CI_Vector(iCSF) =',CI_Vector(iCSF)
      Write (6,*) 'CI_Faroald(iCSF)=',CI_Faroald(iCSF)
      Call Abend()
   End If
End Do
Call mma_deallocate(CI_Faroald)
Endif

#ifdef _SGUGA_VERIFY_
If (.NOT.iDoGAS .and. nRsPrt==1) Then

!  Call SG_TRACI(iState,STSYM,TU,nTU,CI_SGUGA,nCSF)
   Call mma_allocate(ctemp,nCSF,Label='ctemp')
   call SG_REORD(iState,STSYM,0,nCSF,CI_SGUGA,ctemp)
   Call SG_TRACI(iState,STSYM,TU,nTU,ctemp,nCSF)
   call SG_REORD(iState,STSYM,1,nCSF,ctemp,CI_SGUGA)
   Call mma_deallocate(ctemp)

   Do iCSF = 1, nCSF
      If (Abs(CI_Vector(iCSF)-CI_SGUGA(iCSF))/nCSF>1.0E-12_wp) Then
         Write (6,*) 'Mk_Tra_Psi(SGUGA) error, iCSF=',iCSF
         Write (6,*) 'CI_Vector(iCSF) =',CI_Vector(iCSF)
         Write (6,*) 'CI_SGUGA(iCSF)=',CI_SGUGA(iCSF)
         Call Abend()
      End If
   End Do
End If
Call mma_deallocate(CI_SGUGA)
#endif

End subroutine Mk_Tra_Psi

End module CI_Interfaces
