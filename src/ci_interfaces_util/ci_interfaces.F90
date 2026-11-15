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
use sguga, only: SG_Free, SG_ReOrd, SG_Lucia_Walk, SG_Verify_Lucia_Walk
use Lucia_Interface, only: Lucia_Util
use faroald, only: my_norb, sigma_update, htu, gtuvx, ndeta, ndetb ,transition_one_pdm, one_pdm, two_pdm, fold_two_pdm, &
                   traci_faroald
use citrans, only: citrans_conf_dets, citrans_csf2sd, citrans_csf_rank, citrans_sd2csf, citrans_sort
use rasscf_global, only: DoFaro, DoSGEx
use general_data, only: ISPIN, STSYM, nDet, nConf, nAsh, nActEl, NAC=>NLEV
use spinfo, only: NDTASM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use definitions, only: wp, iwp, u6
use sguga, only: CIS, SGS, EXS
use general_data, only: iDoGAS, nRsPrt
use CI_Verify_Options, only: CI_Verify

#include "intent.fh"

Private

Public :: Mk_Tra_Psi, Mk_H_Psi, Mk_H_Psi_Det, Mk_PDMS, Mk_T1DM, Mk_IHI, Mk_IHJ, Mk_Cnf_Layout
Public :: CI_Timer, CI_Initialize, CI_Close, CI_Verify

! CI_Verify (module CI_Verify_Options) is re-exported, see there for the
! run-time verification options (MOLCAS_CI_VERIFY).

! Mk_IHJ is verified against Mk_H_Psi only once per run.
logical(kind=iwp) :: IHJ_Verified = .false.
! sg_h_psi_ras is compared with sg_h_psi (CAS) only once per run; for RAS
! the difference to the CI back end is printed once (checked every time).
logical(kind=iwp) :: HPsi_RAS_Verified = .false.
! The SGUGA CI transformation (Mk_Tra_Psi, SIGMA) is checked on every call,
! its maximum error printed once.
logical(kind=iwp) :: TraPsi_Verified = .false.
! SGEX given but SGUGA not available: warned once.
logical(kind=iwp) :: SGEx_Warned = .false.

contains

!> SGUGA covers the CI space of RASSCF: CAS or RAS, but not GAS (no graph)
!> and not HEXS/DEXS (I_ELIMINATE_GAS), which remove configurations that the
!> graph keeps.
function SG_Avail() result(Avail)
use spinfo, only: I_ELIMINATE_GAS
logical(kind=iwp) :: Avail
Avail = (.not. iDoGAS) .and. (I_ELIMINATE_GAS == 0)
end function SG_Avail

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

integer(kind=iwp) :: IREFSM_SAVE, nDet

real(kind=wp) :: Check_Href, Check_H, DMax
real(kind=wp), allocatable ::  SG_PSI(:), SG_SGM(:), SG_RAS(:)

if (DOFARO) then

  call Faroald_Ints(TU,TUVX)
  call Faroald_Sigma(iState,STSYM,nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,ndeta,ndetb)

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

! CI development (MOLCAS_CI_VERIFY=SIGMA): compare with the SGUGA sigma
! vector, which is then used. With RAS partitions it comes from
! sg_h_psi_ras (no intermediate states outside the RAS space); for CAS from
! sg_h_psi, checked once per run against sg_h_psi_ras. Not with HEXS/DEXS
! (I_ELIMINATE_GAS): the configurations they remove are not part of the
! SGUGA graph.
If (CI_Verify('SIGMA') .and. SG_Avail()) Then
   Call mma_allocate(SG_PSI,nCSF,Label='SG_PSI')
   call SG_Reord(iState,STSYM,0,nCSF,CI_VEC,SG_PSI)
   Check_Href=CheckSum(Sigma_Vec,nCSF)
   Call mma_allocate(SG_SGM,nCSF,Label='SG_SGM')
   if (nRsPrt == 1) then
     Call sg_h_psi(SGS(iState),CIS(iState),EXS(iState),SG_PSI,nCSF,STSYM,SG_SGM,TUVX,nTUVX,TU,nTU)
     if (.not. HPsi_RAS_Verified) then
       call mma_allocate(SG_RAS,nCSF,Label='SG_RAS')
       call sg_h_psi_ras(SGS(iState),CIS(iState),EXS(iState),SG_PSI,nCSF,STSYM,SG_RAS,TUVX,nTUVX,TU,nTU)
       DMax = maxval(abs(SG_RAS(:)-SG_SGM(:)))
       write(u6,'(1X,A,ES10.3)') 'Mk_H_Psi: sg_h_psi_ras vs sg_h_psi, max. abs. difference',DMax
       if (DMax > 1.0e-12_wp*max(One,maxval(abs(SG_SGM(:))))) then
         write(u6,*) 'SGUGA error in sg_h_psi_ras'
         call Abend()
       end if
       call mma_deallocate(SG_RAS)
       HPsi_RAS_Verified = .true.
     end if
   else
     Call sg_h_psi_ras(SGS(iState),CIS(iState),EXS(iState),SG_PSI,nCSF,STSYM,SG_SGM,TUVX,nTUVX,TU,nTU)
   end if
   call SG_Reord(iState,STSYM,1,nCSF,SG_SGM,SG_PSI)
   if (nRsPrt /= 1) then
     DMax = maxval(abs(SG_PSI(:)-Sigma_Vec(:)))
     if (.not. HPsi_RAS_Verified) &
       write(u6,'(1X,A,ES10.3)') 'Mk_H_Psi: sg_h_psi_ras vs CI back end (RAS), max. abs. difference',DMax
     HPsi_RAS_Verified = .true.
     if (DMax > 1.0e-10_wp*max(One,maxval(abs(Sigma_Vec(:))))) then
       write(u6,'(1X,A,ES10.3)') 'SGUGA error in sg_h_psi_ras (RAS), max. abs. difference',DMax
       call Abend()
     end if
   end if
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

End Subroutine Mk_H_Psi

!> Unpack the active one- and two-electron integrals (triangular packing)
!> into the square arrays htu and gtuvx of Faroald.
subroutine Faroald_Ints(TU,TUVX)
real(kind=wp), intent(in) :: TU(:), TUVX(:)
integer(kind=iwp) :: itu, ituvx, it, iu, iv, ixmax, ix

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

end subroutine Faroald_Ints

!> Faroald sigma vector Sigma_Vec = H CI_Vec in Lucia's CSF basis, with the
!> integrals set by Faroald_Ints. ctemp and sigtemp are scratch of at least
!> nDetA*nDetB elements (Faroald determinant basis).
subroutine Faroald_Sigma(iState,STSYM,nCSF,CI_Vec,Sigma_Vec,ctemp,sigtemp,ntemp,ndeta,ndetb)
integer(kind=iwp), intent(in) :: iState, STSYM, nCSF, ntemp, ndeta, ndetb
real(kind=wp), intent(in) :: CI_Vec(nCSF)
real(kind=wp), intent(out) :: Sigma_Vec(nCSF)
real(kind=wp), intent(inout), target :: ctemp(ntemp), sigtemp(ntemp)
real(kind=wp), pointer :: Faroald_PSI(:,:), Faroald_SGM(:,:)

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

end subroutine Faroald_Sigma

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

integer(kind=iwp) :: nCSF
integer(kind=iwp) :: iDet
integer(kind=iwp) :: IREFSM_SAVE

#ifdef _Compute_Full_CI_Matrix
real(kind=wp), allocatable :: Full(:,:)
#endif

!Note that we Lucia code doesn't explicily spin adapt the SD. Rather it just process the input coefficients.
! If we project the starting coefficients to the spin-adapted iCSF basis and then expand it back to the SD basis
! we have effectivly generated a SD coefficient vector that is implicitly spin-adapted. Since the CASVB code generates
! starting CI vectors that are not spin-adapted there will be a difference in the convergence patterns.

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

if (DOFARO) then

  call Faroald_Ints(TU,TUVX)

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

! CI development (MOLCAS_CI_VERIFY=SIGMA): compare with the SGUGA sigma
! vector (CASSCF without symmetry only), which is then used. The
! determinant (CASVB) route needs the SD/CSF codec, which exists for CAS
! only; with RAS partitions the check is skipped.
If (CI_Verify('SIGMA')) Then
  If (SGS(iState)%nRsPrt /= 1) Then
     Return
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
End If

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

End Subroutine Mk_H_Psi_Det

!***********************************************************************************************************************************
!***********************************************************************************************************************************

Subroutine Mk_T1DM(Bra_Vec, Ket_Vec, nVec, T1DM, nT1DM, T1SDM)
use Lucia_Data, only: DTmp, DSTmp
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

End Subroutine Mk_T1DM

!***********************************************************************************************************************************
!***********************************************************************************************************************************

Subroutine Mk_pdms(CIVec,nCIVEC,D,SD,P,PA,nD,nP)
use lucia_data, only: DStmp, Dtmp, PAtmp, PTmp
use rasscf_global, only: NACPAR, NACPR2

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

real(kind=wp) :: Check_D1, Check_P, Check_PA
real(kind=wp), allocatable :: D_Sguga(:)
real(kind=wp), allocatable :: P_Sguga(:), PA_sguga(:)

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

! CI development (MOLCAS_CI_VERIFY=SIGMA): verify the SGUGA density matrices
! against the CI back end. D, P and PA are optional, so only those present
! are compared.
        If (CI_Verify('SIGMA') .and. .NOT.iDoGAS) Then

          Call mma_allocate(CIV,nCIVEC,Label='CIV')
          call SG_Reord(iState,STSYM,0,nCIVEC,CIVEC,CIV)

!         Test the one-particle density matrix
          If (Present(D)) Then
          Check_D1=CheckSum(D,NACPAR)
          Call mma_allocate(D_sguga,NAC*(NAC+1)/2)

          call sg_one_pdm(SGS(istate),CIS(istate),EXS(istate),CIV,SIZE(CIV),STSYM,D_sguga,Size(D_sguga))
          If (ABS(CheckSum(D_sguga,NACPAR)-Check_D1)/SIZE(D_sguga)>1.0e-12_wp) Then
             Check_D1=CheckSum(D_sguga,NACPAR)
             Write (u6,*) 'SGUGA error in D1Mat'
             Call Abend()
          End If
          Call mma_deallocate(D_sguga)
          End If

!         Test the one-particle spin-density matrix
!         This option is not yet developed for the SGUGA code. To come...

!         Test the symmetric two-particle density matrix.
          If (Present(P) .and. Present(PA)) Then
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
          End If
          Call mma_deallocate(CIV)

        END IF

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

! Construct the determinant tables
! Initialize LUCIA and determinant control
call Lucia_Util('Ini')

! to get number of CSFs for GAS
! and number of determinants to store
nconf = MaxVal(ncsasm(1:mxsym))
nDet = MaxVal(ndtasm(1:mxsym))

! Initiate the SGUGA environment conditional to all flags
If (Do_SGUGA) Then
   if (CI_Verify('SIGMA')) then
      write(u6,'(1X,A)') '**EXPERIMENTAL**'
      write(u6,'(1X,A)') 'CI back end cross-checked against SGUGA (MOLCAS_CI_VERIFY=SIGMA).'
      write(u6,'(1X,A)') '**EXPERIMENTAL**'
   end if
   call SG_Setup_RASSCF()
   ! Exhaustive verification of the SGUGA coupling-coefficient codes,
   ! only on request (used by the test/csf_sd groups), see CI_Verify.
   if (CI_Verify('YES')) then
      call SG_Test_EpqErs_Ref(iState,STSym)
   else if (CI_Verify('EPQ')) then
      call SG_Test_EpqErs_Ref(iState,STSym,Epq_Only=.true.)
   end if
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

!> Diagonal of the CI Hamiltonian in the CSF basis as used by the CI solver
!> (preconditioner, selection of the explicit subspace in PHPCSF): for each
!> configuration the average of its determinant diagonals, assigned to all
!> of its CSFs (Lucia, CSDIAG). With DoFaro the same average from Faroald
!> (hdiag), compared with Lucia and then used.
!>
!> With MOLCAS_CI_VERIFY=SIGMA (single RAS partition, no GAS) the exact
!> diagonal <I|H|I> is computed by SGUGA (sg_hdiag_fast, for small spaces
!> checked against sg_hdiag_ref) and compared (Verify_IHI): exactly for the
!> single-determinant configurations, where the average is exact, and
!> otherwise by the configurations that the selection of PHPCSF would pick
!> with either diagonal. The latter is reported, not enforced, since the
!> average is an approximation.
subroutine Mk_IHI(CSFDIA,nConf,TU,nTU,TUVX,nTUVX)
use Lucia_Interface, only: Lucia_Util
use Lucia_data, only: SDREO
use spinfo, only: MINOP, NCNFTP, NCSFTP, NDTFTP, NTYP
use faroald, only: hdiag, ncomb_tot
use Constants, only: One

integer(kind=iwp), intent(in) :: nConf, nTU, nTUVX
real(kind=wp), intent(out):: CSFDIA(nConf)
real(kind=wp), intent(in):: TU(nTU), TUVX(nTUVX)

real(kind=wp), allocatable :: DDIA(:), DIAG(:), CSFDIA_Faroald(:)
integer(kind=iwp):: IPRINT=0

integer(kind=iwp) :: i
real(kind=wp), parameter :: ThrDia = 1.0e-10_wp

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

  call Faroald_Ints(TU,TUVX)

  Call mma_allocate(Diag,nComb_tot,Label='Diag')
  Call mma_allocate(CSFDIA_Faroald,nConf,Label='CSFDIA_F')
  Diag(:)=Zero
  Call hdiag(htu,gtuvx,diag,CSFDIA_Faroald,nConf)
! Call RecPrt('Diag(FAROALD)',' ',Diag,1,nComb_tot)
  Do i=1,nComb_tot
     If (abs(ddia(i)-diag(i)) > ThrDia*max(One,abs(ddia(i)))) Then
        Write (u6,*) 'i, abs(ddia(i)-diag(i))'
        Write (u6,*) i, abs(ddia(i)-diag(i))
        Write (u6,*) ddia(i),diag(i)
        Call RecPrt('DDIA(LUCIA)',' ',DDIA,1,nDet)
        Call RecPrt('Diag(FAROALD)',' ',Diag,1,nComb_tot)
        Call abend()
     End If
  End Do
  Do i=1,nConf
     If (abs(CSFDIA(i)-CSFDIA_Faroald(i)) > ThrDia*max(One,abs(CSFDIA(i)))) Then
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

! CI development (MOLCAS_CI_VERIFY=SIGMA): compare with the exact diagonal.
if (CI_Verify('SIGMA') .and. SG_Avail()) call Verify_IHI()

! CI development (MOLCAS_CI_VERIFY=HDIAG): use the exact diagonal instead.
if (CI_Verify('HDIAG') .and. SG_Avail()) then
  call Exact_Diag(CSFDIA)
  write(u6,'(1X,A)') 'Mk_IHI: exact SGUGA diagonal used (MOLCAS_CI_VERIFY=HDIAG)'
end if

contains

  subroutine Verify_IHI()
    use general_data, only: NSEL
    use sguga, only: SGStruct, CIStruct, EXStruct
    integer(kind=iwp), parameter :: iState = 1, MxRef = 500, MxPrt = 10
    integer(kind=iwp) :: i, iTyp, iCnf, k, nCnf, nC, iC, nSD, nSDErr, nSelL, nSelMin, nSelMean, nBoth, nPrt, iRk(3), iRkMin, &
                         iRkMax
    integer(kind=iwp), allocatable :: nCSFCnf(:), nOpCnf(:)
    real(kind=wp) :: Dev, MaxDev
    real(kind=wp), allocatable :: DSG(:), DRef(:), DEx(:), KeyL(:), KeyMin(:), KeyMean(:)
    logical(kind=iwp), allocatable :: SelL(:), SelMin(:), SelMean(:)
    interface
      subroutine sg_hdiag_fast(SGS,CIS,EXS,nCSFs,CSFSym,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri,Diag)
        import :: SGStruct, CIStruct, EXStruct, iwp, wp
        type(SGStruct), intent(in) :: SGS
        type(CIStruct), intent(in) :: CIS
        type(EXStruct), intent(inout) :: EXS
        integer(kind=iwp), intent(in) :: nCSFs, CSFSym, nTUVX_Tri, nTU_Tri
        real(kind=wp), intent(in) :: TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
        real(kind=wp), intent(out) :: Diag(nCSFs)
      end subroutine sg_hdiag_fast
      subroutine sg_hdiag_ref(SGS,CIS,EXS,nCSFs,CSFSym,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri,Diag)
        import :: SGStruct, CIStruct, EXStruct, iwp, wp
        type(SGStruct), intent(in) :: SGS
        type(CIStruct), intent(in) :: CIS
        type(EXStruct), intent(inout) :: EXS
        integer(kind=iwp), intent(in) :: nCSFs, CSFSym, nTUVX_Tri, nTU_Tri
        real(kind=wp), intent(in) :: TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
        real(kind=wp), intent(out) :: Diag(nCSFs)
      end subroutine sg_hdiag_ref
    end interface

    ! exact diagonal in the SGUGA order
    call mma_allocate(DSG,nConf,Label='DSG')
    call sg_hdiag_fast(SGS(iState),CIS(iState),EXS(iState),nConf,STSYM,TUVX,nTUVX,TU,nTU,DSG)
    if (nConf <= MxRef) then
      call mma_allocate(DRef,nConf,Label='DRef')
      call sg_hdiag_ref(SGS(iState),CIS(iState),EXS(iState),nConf,STSYM,TUVX,nTUVX,TU,nTU,DRef)
      MaxDev = Zero
      do i=1,nConf
        Dev = abs(DSG(i)-DRef(i))
        MaxDev = max(MaxDev,Dev)
        if (Dev > ThrDia*max(One,abs(DRef(i)))) then
          write(u6,*) 'Mk_IHI: sg_hdiag_fast differs from sg_hdiag_ref for CSF',i,DSG(i),DRef(i)
          call Abend()
        end if
      end do
      write(u6,'(1X,A,I0,A,ES10.3)') 'Mk_IHI: sg_hdiag_fast verified against sg_hdiag_ref, ',nConf,' CSFs, max. error ',MaxDev
      call mma_deallocate(DRef)
    end if

    ! to the Lucia CSF order
    call mma_allocate(DEx,nConf,Label='DEx')
    call Diag_to_Lucia(DSG,DEx)
    call mma_deallocate(DSG)

    ! configuration layout (as Mk_Cnf_Layout) with the number of open shells
    ! (NDTFTP counts spin combinations, not determinants, for MS = 0)
    nCnf = sum(NCNFTP(1:NTYP,STSYM))
    call mma_allocate(nCSFCnf,nCnf,Label='nCSFCnf')
    call mma_allocate(nOpCnf,nCnf,Label='nOpCnf')
    k = 0
    do iTyp=1,NTYP
      do iCnf=1,NCNFTP(iTyp,STSYM)
        k = k+1
        nCSFCnf(k) = NCSFTP(iTyp)
        nOpCnf(k) = iTyp-1+MINOP
      end do
    end do

    ! keys per configuration: Lucia (the average, as PHPCSF uses it) and
    ! the exact lowest and mean CSF energies; exact check for single
    ! determinants: closed shells, or all open shells parallel (2S open
    ! shells, the high-spin determinant is the only one with MS = S)
    call mma_allocate(KeyL,nCnf,Label='KeyL')
    call mma_allocate(KeyMin,nCnf,Label='KeyMin')
    call mma_allocate(KeyMean,nCnf,Label='KeyMean')
    nSD = 0
    nSDErr = 0
    MaxDev = Zero
    iC = 0
    do k=1,nCnf
      nC = nCSFCnf(k)
      KeyL(k) = CSFDIA(iC+1)
      KeyMin(k) = minval(DEx(iC+1:iC+nC))
      KeyMean(k) = sum(DEx(iC+1:iC+nC))/real(nC,kind=wp)
      if (nOpCnf(k) == ISPIN-1) then
        nSD = nSD+1
        Dev = abs(CSFDIA(iC+1)-DEx(iC+1))
        MaxDev = max(MaxDev,Dev)
        if (Dev > ThrDia*max(One,abs(DEx(iC+1)))) then
          nSDErr = nSDErr+1
          if (nSDErr <= MxPrt) write(u6,*) 'Mk_IHI: single-determinant configuration',k,' CI solver:',CSFDIA(iC+1), &
                                           ' exact:',DEx(iC+1)
        end if
      end if
      iC = iC+nC
    end do
    if (nSDErr > 0) then
      write(u6,*) 'Mk_IHI: diagonal of single-determinant configurations differs from SGUGA'
      call Abend()
    end if
    write(u6,'(1X,A,I0,A,ES10.3)') 'Mk_IHI: diagonal verified against SGUGA for ',nSD,' single-determinant configurations, &
                                   &max. error ',MaxDev

    ! selection of PHPCSF with each key
    call mma_allocate(SelL,nCnf,Label='SelL')
    call mma_allocate(SelMin,nCnf,Label='SelMin')
    call mma_allocate(SelMean,nCnf,Label='SelMean')
    call Sel_PHP(KeyL,SelL)
    call Sel_PHP(KeyMin,SelMin)
    call Sel_PHP(KeyMean,SelMean)
    nSelL = count(SelL)
    nSelMin = count(SelMin)
    nSelMean = count(SelMean)
    write(u6,'(1X,A,I0,A,I0,A,I0,A,I0,A)') 'Mk_IHI: explicit subspace (NSEL = ',NSEL,' of ',nConf,' CSFs, ',nCnf, &
                                           ' configurations) from the CI-solver diagonal: ',nSelL,' configurations'
    if (NSEL >= nConf) write(u6,'(1X,A)') '        (the whole space, the comparison is trivial)'
    nBoth = count(SelL .and. SelMin)
    write(u6,'(1X,A,I0,A,I0,A,I0,A)') '        exact lowest CSF energy per configuration: ',nSelMin,' configurations, ', &
                                      nBoth,' in common, ',count(SelL .neqv. SelMin),' differ'
    nBoth = count(SelL .and. SelMean)
    write(u6,'(1X,A,I0,A,I0,A,I0,A)') '        exact mean CSF energy per configuration:   ',nSelMean,' configurations, ', &
                                      nBoth,' in common, ',count(SelL .neqv. SelMean),' differ'
    ! The selection stops at the first configuration that does not fit into
    ! NSEL CSFs, so near the cutoff the sizes of the configurations decide;
    ! what matters is whether the differing configurations are ranked near
    ! the cutoff (rank = position in the order of increasing energy).
    if (any(SelL .neqv. SelMin) .or. any(SelL .neqv. SelMean)) then
      write(u6,'(1X,A)') '        differing configurations (selected by CI solver/exact lowest/exact mean, CSFs, '// &
                         'open shells, rank by CI-solver key, by exact lowest, by exact mean):'
      nPrt = 0
      iRkMin = huge(iRkMin)
      iRkMax = 0
      do k=1,nCnf
        if ((SelL(k) .eqv. SelMin(k)) .and. (SelL(k) .eqv. SelMean(k))) cycle
        iRk(1) = Cnf_Rank(KeyL,k)
        iRk(2) = Cnf_Rank(KeyMin,k)
        iRk(3) = Cnf_Rank(KeyMean,k)
        iRkMin = min(iRkMin,minval(iRk))
        iRkMax = max(iRkMax,maxval(iRk))
        nPrt = nPrt+1
        if (nPrt <= MxPrt) write(u6,'(9X,I8,3L2,2I6,3I8)') k,SelL(k),SelMin(k),SelMean(k),nCSFCnf(k),nOpCnf(k),iRk(:)
      end do
      write(u6,'(1X,A,I0,A,I0,A,I0,A)') '        differing configurations are ranked ',iRkMin,' to ',iRkMax, &
                                        ' (',nSelL,' selected by the CI solver)'
    end if

    call mma_deallocate(SelMean)
    call mma_deallocate(SelMin)
    call mma_deallocate(SelL)
    call mma_deallocate(KeyMean)
    call mma_deallocate(KeyMin)
    call mma_deallocate(KeyL)
    call mma_deallocate(nOpCnf)
    call mma_deallocate(nCSFCnf)
    call mma_deallocate(DEx)

  end subroutine Verify_IHI

  ! A diagonal in the SGUGA CSF order to the Lucia order: a permutation,
  ! without the CSF phases (SG_ReOrd of the index vector: Idx(ICSF) = +-ISG).
  subroutine Diag_to_Lucia(DSG,DL)
    real(kind=wp), intent(in) :: DSG(nConf)
    real(kind=wp), intent(out) :: DL(nConf)
    integer(kind=iwp) :: j
    real(kind=wp), allocatable :: Idx(:), Map(:)

    call mma_allocate(Map,nConf,Label='Map')
    call mma_allocate(Idx,nConf,Label='Idx')
    do j=1,nConf
      Map(j) = real(j,kind=wp)
    end do
    call SG_ReOrd(1,STSYM,1,nConf,Map,Idx)
    do j=1,nConf
      DL(j) = DSG(nint(abs(Idx(j))))
    end do
    call mma_deallocate(Idx)
    call mma_deallocate(Map)

  end subroutine Diag_to_Lucia

  ! Exact diagonal <I|H|I> in the Lucia CSF order (sg_hdiag_fast).
  subroutine Exact_Diag(DL)
    use sguga, only: SGStruct, CIStruct, EXStruct
    real(kind=wp), intent(out) :: DL(nConf)
    real(kind=wp), allocatable :: DSG(:)
    interface
      subroutine sg_hdiag_fast(SGS,CIS,EXS,nCSFs,CSFSym,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri,Diag)
        import :: SGStruct, CIStruct, EXStruct, iwp, wp
        type(SGStruct), intent(in) :: SGS
        type(CIStruct), intent(in) :: CIS
        type(EXStruct), intent(inout) :: EXS
        integer(kind=iwp), intent(in) :: nCSFs, CSFSym, nTUVX_Tri, nTU_Tri
        real(kind=wp), intent(in) :: TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
        real(kind=wp), intent(out) :: Diag(nCSFs)
      end subroutine sg_hdiag_fast
    end interface

    call mma_allocate(DSG,nConf,Label='DSG')
    call sg_hdiag_fast(SGS(1),CIS(1),EXS(1),nConf,STSYM,TUVX,nTUVX,TU,nTU,DSG)
    call Diag_to_Lucia(DSG,DL)
    call mma_deallocate(DSG)

  end subroutine Exact_Diag

  ! Position of Key(k) in the order of increasing Key (ties by index).
  function Cnf_Rank(Key,k) result(iRank)
    real(kind=wp), intent(in) :: Key(:)
    integer(kind=iwp), intent(in) :: k
    integer(kind=iwp) :: iRank, l

    iRank = 1
    do l=1,size(Key)
      if ((Key(l) < Key(k)) .or. ((Key(l) == Key(k)) .and. (l < k))) iRank = iRank+1
    end do

  end function Cnf_Rank

  ! Configurations selected by PHPCSF for the explicit subspace with the
  ! configuration energies Key: the lowest first (ties within Acc in the
  ! order of the configurations), until the next one does not fit into NSEL
  ! CSFs.
  subroutine Sel_PHP(Key,Sel)
    use general_data, only: NSEL
    real(kind=wp), intent(in) :: Key(:)
    logical(kind=iwp), intent(out) :: Sel(size(Key))
    integer(kind=iwp) :: k, kMin, nP, nPCnf, nCMin
    real(kind=wp) :: XMax, XMin
    real(kind=wp), allocatable :: Wrk(:)
    real(kind=wp), parameter :: Acc = 1.0e-13_wp
    integer(kind=iwp), allocatable :: nCSFCnf(:)

    call Mk_Cnf_Layout(STSYM,nCSFCnf)
    call mma_allocate(Wrk,size(Key),Label='Wrk')
    Wrk(:) = Key(:)
    XMax = maxval(abs(Key))
    Sel(:) = .false.
    nP = 0
    nPCnf = 0
    do
      XMin = XMax+One
      kMin = 0
      nCMin = 0
      do k=1,size(Key)
        if (Wrk(k)+Acc < XMin) then
          XMin = Wrk(k)
          kMin = k
          nCMin = nCSFCnf(k)
        end if
      end do
      if ((kMin == 0) .or. (nP+nCMin > NSEL)) exit
      Sel(kMin) = .true.
      nP = nP+nCMin
      nPCnf = nPCnf+1
      Wrk(kMin) = XMax+One
      if (nPCnf >= size(Key)) exit
    end do
    call mma_deallocate(Wrk)
    call mma_deallocate(nCSFCnf)

  end subroutine Sel_PHP

End subroutine Mk_IHI

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!> Layout of the CSF basis of symmetry IREFSM by configuration:
!> nCSFCnf(k) is the number of CSFs of the k-th configuration, the
!> configurations being in the order of the CSF basis (the CSFs of one
!> configuration are consecutive). Allocated here, to be freed by the caller
!> with mma_deallocate. (Lucia back end: configuration types from spinfo.)
subroutine Mk_Cnf_Layout(IREFSM,nCSFCnf)
use spinfo, only: NCNFTP, NCSFTP, NTYP

integer(kind=iwp), intent(in) :: IREFSM
integer(kind=iwp), allocatable, intent(out) :: nCSFCnf(:)
integer(kind=iwp) :: iTyp, iCnf, k

call mma_allocate(nCSFCnf,sum(NCNFTP(1:NTYP,IREFSM)),Label='nCSFCnf')
k = 0
do iTyp=1,NTYP
  do iCnf=1,NCNFTP(iTyp,IREFSM)
    k = k+1
    nCSFCnf(k) = NCSFTP(iTyp)
  end do
end do

end subroutine Mk_Cnf_Layout

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!> Explicit CI Hamiltonian <I|H|J> over the CSFs of the configurations
!> IPCNF(1:NPCNF) of symmetry IREFSM (configuration indices as in
!> Mk_Cnf_Layout), packed as a lower triangle: the CSFs are ordered
!> configuration by configuration in the order given, the CSFs of each
!> configuration consecutively. TU and TUVX are the one- and two-electron
!> integrals in triangular packing (as for Mk_IHI); ECORE is added on the
!> diagonal. The Hamiltonian is defined by TU and TUVX alone, as for the
!> sigma vector (Mk_H_Psi): no scaling of the exchange (ExFac), which would
!> only be consistent if built into the integrals by the caller.
!>
!> Back ends:
!>  - Lucia: CNHCN for each pair of configurations (determinant block,
!>    transformed with DTOC).
!>  - Faroald (DoFaro): for each pair of configurations the determinant
!>    block (Slater-Condon rules, Faroald_HDet) transformed with the spin
!>    tables of citrans, C_L^T H C_R.
!>  - SGUGA (RASSCF keyword SGEX; no GAS, single RAS partition): <I|H|J>
!>    for each pair of CSFs from the occupation difference (SG_MatEl_H),
!>    pairs of configurations more than a double excitation apart skipped.
!> With MOLCAS_CI_VERIFY=SIGMA, once per run, the block in use is compared
!> with the other back ends (the Faroald block also with columns H e_J of
!> the Faroald sigma vector) and checked against Mk_H_Psi; the CPU time of
!> each back end is reported.
subroutine Mk_IHJ(PHP,nPHP,IPCNF,NPCNF,IREFSM,ECORE,TU,nTU,TUVX,nTUVX,NTEST)
use Index_Functions, only: nTri_Elem
use Lucia_data, only: Conf_Occ, DFTP, DTOC, IREOTS, Sigma_on_disk
use csfbas, only: NAEL, NBEL
use spinfo, only: NCSFTP, NDTFTP, NTYP
use Constants, only: One

integer(kind=iwp), intent(in) :: nPHP, NPCNF, IPCNF(NPCNF), IREFSM, nTU, nTUVX
real(kind=wp), intent(out) :: PHP(nPHP)
real(kind=wp), intent(in) :: ECORE, TU(nTU), TUVX(nTUVX)
integer(kind=iwp), intent(inout) :: NTEST
integer(kind=iwp) :: NPCSF
integer(kind=iwp), allocatable :: IPCSF(:), nSelCnf(:)
real(kind=wp) :: T0, T1, W
real(kind=wp), allocatable :: PHP_Ref(:)
real(kind=wp), parameter :: Thr = 1.0e-10_wp
logical(kind=iwp) :: SG_OK, UseSG, Verify
character(len=8) :: Lab

! CSF indices (in the CSF basis of symmetry IREFSM) of the selection.
call Sel_CSFs()
NPCSF = size(IPCSF)
if (nTri_Elem(NPCSF) > nPHP) then
  write(u6,*) 'Mk_IHJ: PHP too small for the selected CSFs:',nTri_Elem(NPCSF),nPHP
  call Abend()
end if

! SGUGA is available without GAS and for a single RAS partition (SGS is not
! set up with GAS, so test iDoGAS first).
SG_OK = SG_Avail()
if (DoSGEx .and. (.not. SG_OK) .and. (.not. SGEx_Warned)) then
  write(u6,'(1X,A)') 'Mk_IHJ: SGEX ignored, SGUGA is not available with GAS, HEXS or DEXS; using the CI back end'
  SGEx_Warned = .true.
end if
UseSG = DoSGEx .and. SG_OK
Verify = CI_Verify('SIGMA') .and. (.not. IHJ_Verified)
! the walks of the Lucia CSFs, once per run
if (Verify .and. SG_OK) call SG_Verify_Lucia_Walk(1,IREFSM,nConf)

call CWTime(T0,W)
if (UseSG) then
  Lab = 'SGUGA'
  call IHJ_SGUGA(PHP)
else if (DoFaro) then
  Lab = 'Faroald'
  call IHJ_Faroald(PHP)
else
  Lab = 'Lucia'
  call IHJ_Lucia(PHP)
end if
call CWTime(T1,W)

! CI development (MOLCAS_CI_VERIFY=SIGMA): compare the block, once per run,
! with the other back ends and with exact columns H e_J from Mk_H_Psi.
if (Verify) then
  write(u6,'(1X,3A,I0,A,F10.3,A)') 'Mk_IHJ: ',trim(Lab),' block, ',NPCSF,' CSFs, CPU ',T1-T0,' s'
  call mma_allocate(PHP_Ref,nTri_Elem(NPCSF),Label='PHP_Ref')
  if (UseSG .or. DoFaro) call Ref_Block('Lucia')
  if (DoFaro) then
    if (UseSG) call Ref_Block('Faroald')
    call Ref_Block('Faroald H e_J')
  end if
  if (SG_OK) then
    if (.not. UseSG) call Ref_Block('SGUGA')
    call Verify_SG_Explicit()
  else
    write(u6,'(1X,A)') 'Mk_IHJ: SGUGA block skipped (not available with GAS, HEXS or DEXS)'
  end if
  call mma_deallocate(PHP_Ref)
  IHJ_Verified = .true.
  call Verify_IHJ()
end if

call mma_deallocate(IPCSF)
call mma_deallocate(nSelCnf)

contains

  ! Reference block of back end LabRef in PHP_Ref, compared with PHP.
  subroutine Ref_Block(LabRef)
    character(len=*), intent(in) :: LabRef
    real(kind=wp) :: T2, T3, W2

    call CWTime(T2,W2)
    select case (LabRef)
      case ('Lucia')
        call IHJ_Lucia(PHP_Ref)
      case ('Faroald')
        call IHJ_Faroald(PHP_Ref)
      case ('Faroald H e_J')
        call IHJ_Faroald_Col(PHP_Ref)
      case ('SGUGA')
        call IHJ_SGUGA(PHP_Ref)
    end select
    call CWTime(T3,W2)
    write(u6,'(1X,3A,F10.3,A)') 'Mk_IHJ: ',LabRef,' block, CPU ',T3-T2,' s'
    call Compare_IHJ(trim(Lab),LabRef,PHP_Ref)

  end subroutine Ref_Block

  ! IPCSF: the CSFs of the selected configurations, in the order of PHP;
  ! nSelCnf: the number of CSFs of each selected configuration.
  subroutine Sel_CSFs()
    integer(kind=iwp) :: iCnf, iOff, kCnf, I, ICSF
    integer(kind=iwp), allocatable :: nCSFCnf(:), iCSF1(:)

    call Mk_Cnf_Layout(IREFSM,nCSFCnf)
    call mma_allocate(iCSF1,size(nCSFCnf),Label='iCSF1')
    iOff = 1
    do iCnf=1,size(nCSFCnf)
      iCSF1(iCnf) = iOff
      iOff = iOff+nCSFCnf(iCnf)
    end do
    call mma_allocate(IPCSF,sum(nCSFCnf(IPCNF(1:NPCNF))),Label='IPCSF')
    call mma_allocate(nSelCnf,NPCNF,Label='nSelCnf')
    ICSF = 0
    do kCnf=1,NPCNF
      nSelCnf(kCnf) = nCSFCnf(IPCNF(kCnf))
      do I=1,nCSFCnf(IPCNF(kCnf))
        ICSF = ICSF+1
        IPCSF(ICSF) = iCSF1(IPCNF(kCnf))+I-1
      end do
    end do
    call mma_deallocate(iCSF1)
    call mma_deallocate(nCSFCnf)

  end subroutine Sel_CSFs

  ! Lucia back end: CNHCN for each pair of selected configurations.
  subroutine IHJ_Lucia(H)
    real(kind=wp), intent(out) :: H(nTri_Elem(NPCSF))
    integer(kind=iwp) :: NEL, MXDT, MXCS, nScr, ICNL, ICNR, ILTYP, IRTYP, NCSFL, NCSFR, IILB, IIRB, IIL, IIR, IIRMAX, DL, SL, DR, &
                         SR
    integer(kind=iwp), allocatable :: ICONFL(:), ICONFR(:)
    real(kind=wp), allocatable :: HONE(:,:), CNHCNM(:), Scr(:)

    NEL = NAEL+NBEL
    MXDT = maxval(NDTFTP(1:NTYP))
    MXCS = maxval(NCSFTP(1:NTYP))

    ! One-body matrix in square storage, as CNHCN wants it.
    call mma_allocate(HONE,NAC,NAC,Label='HONE')
    call SQUARE(TU,HONE,NAC,1,NAC)

    ! Scratch for CNHCN: determinant strings and signs of both configurations,
    ! CNFSTR work space, the determinant block, DIHDJ_MOLCAS work space and the
    ! half-transformed block.
    nScr = 2*MXDT*(NEL+1)+NEL*(MXDT+1)+MXDT*MXDT+4*NAC+MXDT+MXCS*MXDT+64
    call mma_allocate(Scr,nScr,Label='IHJ_Scr')
    call mma_allocate(CNHCNM,MXCS*MXCS,Label='CNHCNM')
    call mma_allocate(ICONFL,NEL,Label='ICONFL')
    call mma_allocate(ICONFR,NEL,Label='ICONFR')

    IILB = 1
    do ICNL=1,NPCNF
      call GETCNF(ICONFL,ILTYP,IPCNF(ICNL),Conf_Occ(IREFSM)%A,IREFSM,NEL)
      NCSFL = NCSFTP(ILTYP)
      call Lucia_Occ(ICONFL,ILTYP,DL,SL)
      IIRB = 1
      do ICNR=1,ICNL
        call GETCNF(ICONFR,IRTYP,IPCNF(ICNR),Conf_Occ(IREFSM)%A,IREFSM,NEL)
        NCSFR = NCSFTP(IRTYP)
        call Lucia_Occ(ICONFR,IRTYP,DR,SR)
        if (Cnf_Exc(DL,SL,DR,SR) > 2) then
          CNHCNM(1:NCSFL*NCSFR) = Zero
        else
          call CNHCN(ICONFL,ILTYP,ICONFR,IRTYP,CNHCNM,Scr,NAEL,NBEL,ECORE,HONE,DFTP,DTOC,NAC,TUVX,NTEST,One,IREOTS)
        end if
        do IIL=1,NCSFL
          IIRMAX = merge(IIL,NCSFR,IILB == IIRB)
          do IIR=1,IIRMAX
            H(nTri_Elem(IILB-1+IIL-1)+IIRB-1+IIR) = CNHCNM((IIR-1)*NCSFL+IIL)
          end do
        end do
        IIRB = IIRB+NCSFR
      end do
      IILB = IILB+NCSFL
    end do

    call mma_deallocate(ICONFR)
    call mma_deallocate(ICONFL)
    call mma_deallocate(CNHCNM)
    call mma_deallocate(Scr)
    call mma_deallocate(HONE)

  end subroutine IHJ_Lucia

  ! Doubly (D) and singly (S) occupied orbitals, as bit strings (orbital t in
  ! bit t-1), of the Lucia configuration ICONF of type ITYP (GETCNF form:
  ! the closed orbitals first, then the open ones).
  subroutine Lucia_Occ(ICONF,ITYP,D,S)
    use spinfo, only: MINOP
    integer(kind=iwp), intent(in) :: ICONF(*), ITYP
    integer(kind=iwp), intent(out) :: D, S
    integer(kind=iwp) :: IOP, ICL, I

    IOP = ITYP-1+MINOP
    ICL = (NAEL+NBEL-IOP)/2
    D = 0
    S = 0
    do I=1,ICL
      D = ibset(D,ICONF(I)-1)
    end do
    do I=ICL+1,ICL+IOP
      S = ibset(S,ICONF(I)-1)
    end do

  end subroutine Lucia_Occ

  ! Faroald back end: for each pair of selected configurations L, R the
  ! determinant block <D_L|H|D_R> (Faroald_HDet) is transformed to the CSFs
  ! with the spin tables of citrans, C_L^T H C_R; pairs of configurations
  ! more than a double excitation apart are skipped (Cnf_Exc). The selected
  ! CSFs of the CI vector (Lucia's CSF basis) are mapped to the
  ! configuration order of citrans, with the phases of SG_ReOrd and
  ! citrans_sort, as in Faroald_Sigma (Map_Sel).
  subroutine IHJ_Faroald(H)
    real(kind=wp), intent(out) :: H(nTri_Elem(NPCSF))
    integer(kind=iwp) :: I, J, kCnf, ICNL, ICNR, IILB, IIRB, IIL, IIR, IIRMAX, k1, nc, nd, nDTot, nCTot, nL, nR, ndL, ndR
    integer(kind=iwp), allocatable :: kSel(:), iSel(:), iDOff(:), iCOff(:), nDCnf(:), DA(:), DB(:), da1(:), db1(:), &
                                      DOcc(:), SOcc(:)
    real(kind=wp), allocatable :: C(:), cf(:,:), HDet(:,:), T(:,:), Blk(:,:)

    call Faroald_Ints(TU,TUVX)

    ! selected CSF J <-> citrans CSF kSel(J), with phase iSel(J)
    call mma_allocate(kSel,NPCSF,Label='kSel')
    call mma_allocate(iSel,NPCSF,Label='iSel')
    call Map_Sel(kSel,iSel)

    ! Determinants and CSF coefficients of the selected configurations.
    call mma_allocate(nDCnf,NPCNF,Label='nDCnf')
    call mma_allocate(iDOff,NPCNF,Label='iDOff')
    call mma_allocate(iCOff,NPCNF,Label='iCOff')
    call mma_allocate(DOcc,NPCNF,Label='DOcc')
    call mma_allocate(SOcc,NPCNF,Label='SOcc')
    nDTot = 0
    nCTot = 0
    IILB = 0
    do kCnf=1,NPCNF
      call citrans_conf_dets(kSel(IILB+1),k1,nc,nd,da1,db1,cf)
      if (nc /= nSelCnf(kCnf)) then
        write(u6,*) 'Mk_IHJ: different number of CSFs in Lucia and citrans:',nSelCnf(kCnf),nc
        call Abend()
      end if
      nDCnf(kCnf) = nd
      ! doubly and singly occupied orbitals of the configuration
      DOcc(kCnf) = iand(da1(1),db1(1))
      SOcc(kCnf) = ieor(da1(1),db1(1))
      iDOff(kCnf) = nDTot
      iCOff(kCnf) = nCTot
      nDTot = nDTot+nd
      nCTot = nCTot+nd*nc
      IILB = IILB+nc
      call mma_deallocate(da1)
      call mma_deallocate(db1)
      call mma_deallocate(cf)
    end do
    call mma_allocate(DA,nDTot,Label='DA')
    call mma_allocate(DB,nDTot,Label='DB')
    call mma_allocate(C,nCTot,Label='C')
    IILB = 0
    do kCnf=1,NPCNF
      call citrans_conf_dets(kSel(IILB+1),k1,nc,nd,da1,db1,cf)
      DA(iDOff(kCnf)+1:iDOff(kCnf)+nd) = da1(:)
      DB(iDOff(kCnf)+1:iDOff(kCnf)+nd) = db1(:)
      do I=1,nc
        J = kSel(IILB+I)-k1+1
        if ((J < 1) .or. (J > nc)) then
          write(u6,*) 'Mk_IHJ: CSF of the selection not in its configuration:',IPCSF(IILB+I)
          call Abend()
        end if
        C(iCOff(kCnf)+(I-1)*nd+1:iCOff(kCnf)+I*nd) = iSel(IILB+I)*cf(:,J)
      end do
      IILB = IILB+nc
      call mma_deallocate(da1)
      call mma_deallocate(db1)
      call mma_deallocate(cf)
    end do
    call mma_deallocate(iSel)
    call mma_deallocate(kSel)

    ! Blocks C_L^T <D_L|H|D_R> C_R.
    nd = maxval(nDCnf)
    nc = maxval(nSelCnf)
    call mma_allocate(HDet,nd,nd,Label='HDet')
    call mma_allocate(T,nd,nc,Label='T')
    call mma_allocate(Blk,nc,nc,Label='Blk')
    IILB = 1
    do ICNL=1,NPCNF
      nL = nSelCnf(ICNL)
      ndL = nDCnf(ICNL)
      IIRB = 1
      do ICNR=1,ICNL
        nR = nSelCnf(ICNR)
        ndR = nDCnf(ICNR)
        if (Cnf_Exc(DOcc(ICNL),SOcc(ICNL),DOcc(ICNR),SOcc(ICNR)) > 2) then
          Blk(1:nL,1:nR) = Zero
        else
          do J=1,ndR
            do I=1,ndL
              HDet(I,J) = Faroald_HDet(DA(iDOff(ICNL)+I),DB(iDOff(ICNL)+I),DA(iDOff(ICNR)+J),DB(iDOff(ICNR)+J))
            end do
          end do
          call DGEMM_('N','N',ndL,nR,ndR,One,HDet,nd,C(iCOff(ICNR)+1:iCOff(ICNR)+ndR*nR),ndR,Zero,T,nd)
          call DGEMM_('T','N',nL,nR,ndL,One,C(iCOff(ICNL)+1:iCOff(ICNL)+ndL*nL),ndL,T,nd,Zero,Blk,nc)
        end if
        do IIL=1,nL
          IIRMAX = merge(IIL,nR,IILB == IIRB)
          do IIR=1,IIRMAX
            H(nTri_Elem(IILB-1+IIL-1)+IIRB-1+IIR) = Blk(IIL,IIR)
          end do
        end do
        IIRB = IIRB+nR
      end do
      IILB = IILB+nL
    end do
    do I=1,NPCSF
      H(nTri_Elem(I)) = H(nTri_Elem(I))+ECORE
    end do

    call mma_deallocate(Blk)
    call mma_deallocate(T)
    call mma_deallocate(HDet)
    call mma_deallocate(C)
    call mma_deallocate(DB)
    call mma_deallocate(DA)
    call mma_deallocate(SOcc)
    call mma_deallocate(DOcc)
    call mma_deallocate(iCOff)
    call mma_deallocate(iDOff)
    call mma_deallocate(nDCnf)

  end subroutine IHJ_Faroald

  ! Selected CSF J (Lucia CSF IPCSF(J)) <-> citrans CSF kSel(J), with phase
  ! iSel(J): cinew(kSel(J)) = iSel(J)*ciold(IPCSF(J)) for the reordering
  ! SG_ReOrd (mode 0) followed by citrans_sort (mode 'C'). For a single RAS
  ! partition from the step vectors of the selected CSFs alone, otherwise
  ! by reordering the index vector 1..nConf (Map_All). With
  ! MOLCAS_CI_VERIFY=SIGMA the two are compared (once per run).
  subroutine Map_Sel(kSel,iSel)
    use Molcas, only: MxAct
    integer(kind=iwp), intent(out) :: kSel(NPCSF), iSel(NPCSF)
    integer(kind=iwp) :: J, IP, iPh
    integer(kind=iwp) :: IWALK(MxAct)
    integer(kind=iwp), allocatable :: kRef(:), iRef(:)

    if (SGS(1)%nRsPrt /= 1) then
      call Map_All(kSel,iSel)
      return
    end if
    do J=1,NPCSF
      call SG_Lucia_Walk(1,IREFSM,IPCSF(J),IWALK,IP)
      call citrans_csf_rank(IWALK(1:NAC),kSel(J),iPh)
      iSel(J) = IP*iPh
    end do

    if (CI_Verify('SIGMA') .and. (.not. IHJ_Verified)) then
      call mma_allocate(kRef,NPCSF,Label='kRef')
      call mma_allocate(iRef,NPCSF,Label='iRef')
      call Map_All(kRef,iRef)
      if (any(kSel(:) /= kRef(:)) .or. any(iSel(:) /= iRef(:))) then
        write(u6,*) 'Mk_IHJ: mapping of the selected CSFs to citrans differs from SG_ReOrd + citrans_sort'
        call Abend()
      end if
      write(u6,'(1X,A,I0,A)') 'Mk_IHJ: mapping of ',NPCSF,' selected CSFs to citrans verified'
      call mma_deallocate(iRef)
      call mma_deallocate(kRef)
    end if

  end subroutine Map_Sel

  ! Map_Sel by reordering the full index vector 1..nConf.
  subroutine Map_All(kSel,iSel)
    integer(kind=iwp), intent(out) :: kSel(NPCSF), iSel(NPCSF)
    integer(kind=iwp) :: I, K
    integer(kind=iwp), allocatable :: kAll(:), iAll(:)
    real(kind=wp), allocatable :: Idx(:), Map(:)

    call mma_allocate(Idx,nConf,Label='Idx')
    call mma_allocate(Map,nConf,Label='Map')
    do I=1,nConf
      Map(I) = real(I,kind=wp)
    end do
    call SG_REORD(1,IREFSM,0,nConf,Map,Idx)
    call CITRANS_SORT('C',Idx,Map)
    call mma_allocate(kAll,nConf,Label='kAll')
    call mma_allocate(iAll,nConf,Label='iAll')
    do K=1,nConf
      I = nint(abs(Map(K)))
      kAll(I) = K
      iAll(I) = nint(sign(One,Map(K)))
    end do
    kSel(:) = kAll(IPCSF(:))
    iSel(:) = iAll(IPCSF(:))
    call mma_deallocate(iAll)
    call mma_deallocate(kAll)
    call mma_deallocate(Map)
    call mma_deallocate(Idx)

  end subroutine Map_All

  ! SGUGA back end: <I|H|J> = P_I P_J <I_G|H|J_G> + delta_IJ ECORE, with
  ! I_G the GUGA CSF of the Lucia CSF I and P_I the phase of SG_ReOrd
  ! (mode 0: c_G(I_G) = P_I c(I)); the block <I_G|H|J_G> from SG_Block_H
  ! (pairs of configurations more than a double excitation apart skipped).
  subroutine IHJ_SGUGA(H)
    use Molcas, only: MxAct
    use sguga, only: SG_Block_H, SG_MatEl_Init, SG_NUM
    real(kind=wp), intent(out) :: H(nTri_Elem(NPCSF))
    integer(kind=iwp) :: I, J
    integer(kind=iwp) :: IWALK(MxAct)
    integer(kind=iwp), allocatable :: iG(:), iPG(:)

    call SG_MatEl_Init(1,IREFSM)
    call mma_allocate(iG,NPCSF,Label='iG')
    call mma_allocate(iPG,NPCSF,Label='iPG')
    do J=1,NPCSF
      call SG_Lucia_Walk(1,IREFSM,IPCSF(J),IWALK,iPG(J))
      iG(J) = SG_NUM(SGS(1),EXS(1),IWALK)
    end do
    call SG_Block_H(IREFSM,NPCNF,nSelCnf,NPCSF,iG,TU,nTU,TUVX,nTUVX,H)
    do I=1,NPCSF
      do J=1,I
        H(nTri_Elem(I-1)+J) = iPG(I)*iPG(J)*H(nTri_Elem(I-1)+J)
      end do
      H(nTri_Elem(I)) = H(nTri_Elem(I))+ECORE
    end do
    call mma_deallocate(iPG)
    call mma_deallocate(iG)

  end subroutine IHJ_SGUGA

  ! CI development (SIGMA): SGUGA's own configuration list and explicit
  ! subspace (SG_Mk_Cnf_List, SG_Select_Explicit, SG_Explicit_Block) against
  ! the CI back end: (1) each configuration of the back end is exactly one
  ! SGUGA configuration (same CSFs), (2) the two selections are compared
  ! (the keys differ: the configuration average of the back end, the mean
  ! exact CSF diagonal of SGUGA), (3) the SGUGA block over the configurations
  ! selected by the back end has the eigenvalues of the block in use.
  subroutine Verify_SG_Explicit()
    use Constants, only: Two
    use sguga, only: SGCnfList, SG_Mk_Cnf_List, SG_Free_Cnf_List, SG_Select_Explicit, SG_Explicit_Block, SGStruct, CIStruct, &
                     EXStruct
    type(SGCnfList) :: Cnf
    integer(kind=iwp) :: I, k, kG, iOffL, MxDim, nGSel, nGCSF, nBoth, nLSel
    integer(kind=iwp), allocatable :: nCSFL(:), LtoG(:), CnfOfG(:), LtoGCnf(:), GSel(:), GCSF(:), GtoL(:)
    logical(kind=iwp), allocatable :: SelG(:), SelL(:)
    real(kind=wp) :: Dev
    real(kind=wp), allocatable :: Idx(:), Map(:), DSG(:), HG(:), EG(:), EL(:)
    interface
      subroutine sg_hdiag_fast(SGS,CIS,EXS,nCSFs,CSFSym,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri,Diag)
        import :: SGStruct, CIStruct, EXStruct, iwp, wp
        type(SGStruct), intent(in) :: SGS
        type(CIStruct), intent(in) :: CIS
        type(EXStruct), intent(inout) :: EXS
        integer(kind=iwp), intent(in) :: nCSFs, CSFSym, nTUVX_Tri, nTU_Tri
        real(kind=wp), intent(in) :: TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
        real(kind=wp), intent(out) :: Diag(nCSFs)
      end subroutine sg_hdiag_fast
    end interface

    call SG_Mk_Cnf_List(1,IREFSM,Cnf)
    ! back-end CSF -> split-graph CSF, and the SGUGA configuration of each
    call mma_allocate(Idx,nConf,Label='Idx')
    call mma_allocate(Map,nConf,Label='Map')
    do I=1,nConf
      Idx(I) = real(I,kind=wp)
    end do
    call SG_ReOrd(1,IREFSM,0,nConf,Idx,Map)
    call mma_allocate(LtoG,nConf,Label='LtoG')
    do I=1,nConf
      LtoG(nint(abs(Map(I)))) = I
    end do
    call mma_deallocate(Map)
    call mma_deallocate(Idx)
    call mma_allocate(CnfOfG,nConf,Label='CnfOfG')
    do kG=1,Cnf%nCnf
      CnfOfG(Cnf%iCSF(Cnf%iOff(kG)+1:Cnf%iOff(kG+1))) = kG
    end do
    ! (1) configurations of the back end -> SGUGA configurations, a bijection
    call Mk_Cnf_Layout(IREFSM,nCSFL)
    if (size(nCSFL) /= Cnf%nCnf) then
      write(u6,*) 'Mk_IHJ: number of configurations, CI back end and SGUGA:',size(nCSFL),Cnf%nCnf
      call Abend()
    end if
    call mma_allocate(LtoGCnf,size(nCSFL),Label='LtoGCnf')
    call mma_allocate(SelG,Cnf%nCnf,Label='SelG')
    SelG(:) = .false.
    iOffL = 0
    do k=1,size(nCSFL)
      kG = CnfOfG(LtoG(iOffL+1))
      if (any(CnfOfG(LtoG(iOffL+1:iOffL+nCSFL(k))) /= kG) .or. (nCSFL(k) /= Cnf%iOff(kG+1)-Cnf%iOff(kG)) .or. SelG(kG)) then
        write(u6,*) 'Mk_IHJ: configuration',k,' of the CI back end is not one SGUGA configuration'
        call Abend()
      end if
      SelG(kG) = .true.
      LtoGCnf(k) = kG
      iOffL = iOffL+nCSFL(k)
    end do
    write(u6,'(1X,A,I0,A)') 'Mk_IHJ: SGUGA configuration list verified, ',Cnf%nCnf,' configurations as in the CI back end'
    ! (2) the SGUGA selection with the exact diagonal
    MxDim = nint((sqrt(real(8*nPHP+1,kind=wp))-One)/Two)
    call mma_allocate(DSG,nConf,Label='DSG')
    call sg_hdiag_fast(SGS(1),CIS(1),EXS(1),nConf,IREFSM,TUVX,nTUVX,TU,nTU,DSG)
    call mma_allocate(GSel,Cnf%nCnf,Label='GSel')
    call mma_allocate(GCSF,max(1,MxDim),Label='GCSF')
    call SG_Select_Explicit(Cnf,DSG,MxDim,nGSel,GSel,nGCSF,GCSF)
    call mma_deallocate(DSG)
    call mma_allocate(SelL,Cnf%nCnf,Label='SelL')
    SelG(:) = .false.
    SelL(:) = .false.
    SelG(GSel(1:nGSel)) = .true.
    SelL(LtoGCnf(IPCNF(1:NPCNF))) = .true.
    nBoth = count(SelG .and. SelL)
    nLSel = count(SelL)
    write(u6,'(1X,A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'Mk_IHJ: SGUGA selection (mean exact CSF diagonal): ',nGSel, &
          ' configurations, ',nGCSF,' CSFs; CI back end: ',nLSel,' configurations; ',nBoth,' in common, ', &
          count(SelG .neqv. SelL),' differ'
    ! (3) SGUGA block over the configurations selected by the back end
    call mma_allocate(GtoL,max(1,NPCNF),Label='GtoL')
    GtoL(1:NPCNF) = LtoGCnf(IPCNF(1:NPCNF))
    call mma_allocate(HG,nTri_Elem(NPCSF),Label='HG')
    call SG_Explicit_Block(1,Cnf,NPCNF,GtoL,ECORE,TU,nTU,TUVX,nTUVX,HG)
    call mma_allocate(EG,NPCSF,Label='EG')
    call mma_allocate(EL,NPCSF,Label='EL')
    call Block_Eigenvalues(HG,NPCSF,EG)
    HG(:) = PHP(1:nTri_Elem(NPCSF))
    call Block_Eigenvalues(HG,NPCSF,EL)
    Dev = maxval(abs(EG(:)-EL(:)))
    write(u6,'(1X,A,I0,A,ES10.3)') 'Mk_IHJ: SG_Explicit_Block, eigenvalues verified against the block in use, ',NPCSF, &
                                   ' CSFs, max. error ',Dev
    if (Dev > Thr*max(One,maxval(abs(EL(:))))) then
      write(u6,*) 'Mk_IHJ: eigenvalues of the SGUGA explicit block differ'
      call Abend()
    end if
    call mma_deallocate(EL)
    call mma_deallocate(EG)
    call mma_deallocate(HG)
    call mma_deallocate(GtoL)
    call mma_deallocate(SelL)
    call mma_deallocate(GCSF)
    call mma_deallocate(GSel)
    call mma_deallocate(SelG)
    call mma_deallocate(LtoGCnf)
    call mma_deallocate(nCSFL)
    call mma_deallocate(CnfOfG)
    call mma_deallocate(LtoG)
    call SG_Free_Cnf_List(Cnf)

  end subroutine Verify_SG_Explicit

  ! Eigenvalues (increasing) of a lower-triangle packed symmetric matrix;
  ! the matrix is overwritten.
  subroutine Block_Eigenvalues(A,n,E)
    integer(kind=iwp), intent(in) :: n
    real(kind=wp), intent(inout) :: A(nTri_Elem(n))
    real(kind=wp), intent(out) :: E(n)
    integer(kind=iwp) :: I
    real(kind=wp), allocatable :: V(:,:)
    call mma_allocate(V,n,n,Label='V')
    V(:,:) = Zero
    do I=1,n
      V(I,I) = One
    end do
    call NIdiag_new(A,V,n,n)
    call JACORD(A,V,n,n)
    do I=1,n
      E(I) = A(nTri_Elem(I))
    end do
    call mma_deallocate(V)
  end subroutine Block_Eigenvalues

  ! Faroald reference (CI development): H(I,J) = (H e_J)(I) from the Faroald
  ! sigma vector of the unit vector e_J, J over the selected CSFs.
  subroutine IHJ_Faroald_Col(H)
    real(kind=wp), intent(out) :: H(nTri_Elem(NPCSF))
    integer(kind=iwp) :: I, J, nTmp
    real(kind=wp), allocatable :: EJ(:), HEJ(:), cTmp(:), sTmp(:)

    nTmp = nDetA*nDetB
    call mma_allocate(cTmp,nTmp,Label='cTmp')
    call mma_allocate(sTmp,nTmp,Label='sTmp')
    call mma_allocate(EJ,nConf,Label='EJ')
    call mma_allocate(HEJ,nConf,Label='HEJ')

    call Faroald_Ints(TU,TUVX)
    EJ(:) = Zero
    do J=1,NPCSF
      EJ(IPCSF(J)) = One
      call Faroald_Sigma(1,IREFSM,nConf,EJ,HEJ,cTmp,sTmp,nTmp,nDetA,nDetB)
      EJ(IPCSF(J)) = Zero
      do I=J,NPCSF
        H(nTri_Elem(I-1)+J) = HEJ(IPCSF(I))
      end do
      H(nTri_Elem(J)) = H(nTri_Elem(J))+ECORE
    end do

    call mma_deallocate(HEJ)
    call mma_deallocate(EJ)
    call mma_deallocate(sTmp)
    call mma_deallocate(cTmp)

  end subroutine IHJ_Faroald_Col

  ! Compare PHP (back end Lab) with the reference block HRef (back end LabRef).
  subroutine Compare_IHJ(Lab,LabRef,HRef)
    character(len=*), intent(in) :: Lab, LabRef
    real(kind=wp), intent(in) :: HRef(nTri_Elem(NPCSF))
    integer(kind=iwp) :: I, J, IJ
    real(kind=wp) :: ErrMax

    ErrMax = Zero
    do I=1,NPCSF
      do J=1,I
        IJ = nTri_Elem(I-1)+J
        ErrMax = max(ErrMax,abs(PHP(IJ)-HRef(IJ)))
        if (abs(PHP(IJ)-HRef(IJ)) > Thr) then
          write(u6,*) 'Mk_IHJ: <I|H|J> of ',Lab,' differs from ',LabRef
          write(u6,*) 'I, J, CSF(I), CSF(J) =',I,J,IPCSF(I),IPCSF(J)
          write(u6,*) Lab,' =',PHP(IJ)
          write(u6,*) LabRef,' =',HRef(IJ)
          call Abend()
        end if
      end do
    end do
    write(u6,'(1X,5A,I0,A,ES10.3)') 'Mk_IHJ: ',Lab,' verified against ',LabRef,', ',NPCSF,' CSFs, max. error ',ErrMax

  end subroutine Compare_IHJ

  ! Compare PHP with columns H e_J computed by Mk_H_Psi (the sigma vector of
  ! the CI back end, with SIGMA the one of SGUGA) for at most mxCol CSFs J
  ! of the selection, evenly spaced.
  subroutine Verify_IHJ()
    integer(kind=iwp), parameter :: mxCol = 32
    integer(kind=iwp) :: I, J, JJ, nCol, nTmp
    real(kind=wp), allocatable :: EJ(:), HEJ(:), cTmp(:), sTmp(:)
    real(kind=wp) :: Ref, Err, ErrMax
    logical(kind=iwp) :: Sigma_Save

    if (DoFaro) then
      nTmp = nDetA*nDetB
    else
      nTmp = nDet
    end if
    call mma_allocate(cTmp,nTmp,Label='cTmp')
    call mma_allocate(sTmp,nTmp,Label='sTmp')
    call mma_allocate(EJ,nConf,Label='EJ')
    ! CSDTVC may write beyond nConf (cf. David5)
    call mma_allocate(HEJ,max(nConf,nDet),Label='HEJ')

    Sigma_Save = Sigma_on_disk
    nCol = min(NPCSF,mxCol)
    ErrMax = Zero
    do JJ=1,nCol
      J = 1+((JJ-1)*(NPCSF-1))/max(nCol-1,1)
      EJ(:) = Zero
      EJ(IPCSF(J)) = One
      call Mk_H_Psi(1,IREFSM,nConf,EJ,HEJ,cTmp,sTmp,nTmp,nDetA,nDetB,nTU,TU,nTUVX,TUVX)
      do I=1,NPCSF
        Ref = HEJ(IPCSF(I))
        if (I == J) Ref = Ref+ECORE
        Err = abs(PHP(nTri_Elem(max(I,J)-1)+min(I,J))-Ref)
        ErrMax = max(ErrMax,Err)
        if (Err > Thr) then
          write(u6,*) 'Mk_IHJ: <I|H|J> differs from (H e_J)(I)'
          write(u6,*) 'I, J, CSF(I), CSF(J) =',I,J,IPCSF(I),IPCSF(J)
          write(u6,*) 'Mk_IHJ   =',PHP(nTri_Elem(max(I,J)-1)+min(I,J))
          write(u6,*) 'Mk_H_Psi =',Ref
          call Abend()
        end if
      end do
    end do
    Sigma_on_disk = Sigma_Save
    write(u6,'(1X,A,I0,A,I0,A,ES10.3)') 'Mk_IHJ verified against Mk_H_Psi: ',nCol,' of ',NPCSF,' columns, max. error ',ErrMax

    call mma_deallocate(HEJ)
    call mma_deallocate(EJ)
    call mma_deallocate(sTmp)
    call mma_deallocate(cTmp)

  end subroutine Verify_IHJ

end subroutine Mk_IHJ

!> Excitation level between two configurations given by their doubly (D)
!> and singly (S) occupied orbitals as bit strings: the number of electrons
!> to move, Sum_t max(0,n_t(L)-n_t(R)). H couples only configurations with
!> an excitation level of at most 2.
function Cnf_Exc(DL,SL,DR,SR) result(nExc)
integer(kind=iwp) :: nExc
integer(kind=iwp), intent(in) :: DL, SL, DR, SR

! n_t(L)-n_t(R) = 2 (doubly -> empty), 1 (doubly -> singly, singly -> empty)
nExc = 2*popcnt(iand(DL,not(ior(DR,SR))))+popcnt(iand(DL,SR))+popcnt(iand(SL,not(ior(DR,SR))))

end function Cnf_Exc

!> <aL bL|H|aR bR> between Faroald determinants (alpha and beta strings as
!> bit strings, orbital t in bit t-1, alpha string before the beta string)
!> by the Slater-Condon rules, with the integrals htu and gtuvx (chemists'
!> notation) of Faroald, H = sum_tu h_tu E_tu + 1/2 sum_tuvx g_tuvx
!> (E_tu E_vx - delta_uv E_tx), as in sigma_update.
function Faroald_HDet(aL,bL,aR,bR) result(E)
use Constants, only: Half

real(kind=wp) :: E
integer(kind=iwp), intent(in) :: aL, bL, aR, bR
integer(kind=iwp) :: dA, dB, nA, nB, i, j, p, q, r, s, aI

dA = ieor(aL,aR)
dB = ieor(bL,bR)
nA = popcnt(dA)/2
nB = popcnt(dB)/2
E = Zero
if (nA+nB > 2) return

if (nA+nB == 0) then
  do i=1,my_norb
    if (btest(aL,i-1)) then
      E = E+htu(i,i)
      do j=1,my_norb
        if (btest(aL,j-1)) E = E+Half*(gtuvx(i,i,j,j)-gtuvx(i,j,j,i))
        if (btest(bL,j-1)) E = E+gtuvx(i,i,j,j)
      end do
    end if
    if (btest(bL,i-1)) then
      E = E+htu(i,i)
      do j=1,my_norb
        if (btest(bL,j-1)) E = E+Half*(gtuvx(i,i,j,j)-gtuvx(i,j,j,i))
      end do
    end if
  end do
else if (nA == 1 .and. nB == 0) then
  call Single(aL,aR,bR)
else if (nA == 0 .and. nB == 1) then
  call Single(bL,bR,aR)
else if (nA == 2) then
  call Double(aL,aR)
else if (nB == 2) then
  call Double(bL,bR)
else
  ! one alpha (p <- q) and one beta (r <- s) excitation
  p = Orb(iand(dA,aL),1)
  q = Orb(iand(dA,aR),1)
  r = Orb(iand(dB,bL),1)
  s = Orb(iand(dB,bR),1)
  E = Par(aR,p,q)*Par(bR,r,s)*gtuvx(p,q,r,s)
end if

contains

  ! single excitation p <- q in the string xL/xR of one spin, yR the string
  ! of the other spin
  subroutine Single(xL,xR,yR)
    integer(kind=iwp), intent(in) :: xL, xR, yR
    integer(kind=iwp) :: dX, iC, j

    dX = ieor(xL,xR)
    p = Orb(iand(dX,xL),1)
    q = Orb(iand(dX,xR),1)
    iC = iand(xL,xR)
    E = htu(p,q)
    do j=1,my_norb
      if (btest(iC,j-1)) E = E+gtuvx(p,q,j,j)-gtuvx(p,j,j,q)
      if (btest(yR,j-1)) E = E+gtuvx(p,q,j,j)
    end do
    E = Par(xR,p,q)*E

  end subroutine Single

  ! double excitation p,r <- q,s in the string xL/xR of one spin
  subroutine Double(xL,xR)
    integer(kind=iwp), intent(in) :: xL, xR
    integer(kind=iwp) :: dX

    dX = ieor(xL,xR)
    p = Orb(iand(dX,xL),1)
    r = Orb(iand(dX,xL),2)
    q = Orb(iand(dX,xR),1)
    s = Orb(iand(dX,xR),2)
    aI = ibset(ibclr(xR,q-1),p-1)
    E = Par(xR,p,q)*Par(aI,r,s)*(gtuvx(p,q,r,s)-gtuvx(p,s,r,q))

  end subroutine Double

  ! the n-th orbital (bit) set in x
  function Orb(x,n) result(k)
    integer(kind=iwp) :: k
    integer(kind=iwp), intent(in) :: x, n
    integer(kind=iwp) :: m

    m = 0
    do k=1,my_norb
      if (btest(x,k-1)) then
        m = m+1
        if (m == n) return
      end if
    end do
    k = 0

  end function Orb

  ! sign of a+_t a_u acting on the string x: the parity of the number of
  ! orbitals occupied in x strictly between t and u
  function Par(x,t,u) result(iSgn)
    integer(kind=iwp) :: iSgn
    integer(kind=iwp), intent(in) :: x, t, u
    integer(kind=iwp) :: lo, hi

    lo = min(t,u)
    hi = max(t,u)
    iSgn = 1-2*poppar(ibits(x,lo,hi-lo-1))

  end function Par

end function Faroald_HDet

subroutine Mk_Tra_Psi(iState,TU,nTU,CI_Vector,nCSF)

use sguga, only: sg_traci

integer(kind=iwp), intent(in) ::  iState,nTU, nCSF
real(kind=wp), intent(in) :: TU(nTU)
real(kind=wp), intent(inout) :: CI_Vector(nCSF)

real(kind=wp), allocatable :: Faroald_Psi(:,:), ctemp(:)
real(kind=wp), allocatable :: dtemp(:)

integer(kind=iwp) :: iCSF
real(kind=wp), allocatable :: CI_Faroald(:)
real(kind=wp), allocatable :: CI_SGUGA(:)
logical(kind=iwp) :: Verify_SG

! CI development (MOLCAS_CI_VERIFY=SIGMA): also transform with SGUGA and compare.
Verify_SG = CI_Verify('SIGMA') .and. SG_Avail()
If (Verify_SG) Then
  Call mma_allocate(CI_SGUGA,nCSF,Label='CI_SGUGA')
  CI_SGUGA(:) = CI_Vector(:)
End If

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

If (Verify_SG) Then

!  Call SG_TRACI(iState,STSYM,TU,nTU,CI_SGUGA,nCSF)
   Call mma_allocate(ctemp,nCSF,Label='ctemp')
   call SG_REORD(iState,STSYM,0,nCSF,CI_SGUGA,ctemp)
   Call SG_TRACI(iState,STSYM,TU,nTU,ctemp,nCSF)
   call SG_REORD(iState,STSYM,1,nCSF,ctemp,CI_SGUGA)
   Call mma_deallocate(ctemp)
   If (.not. TraPsi_Verified) Then
      Write (u6,'(1X,A,I0,A,ES10.3)') 'Mk_Tra_Psi: SGUGA CI transformation verified against the CI back end, ',nCSF, &
                                     ' CSFs, max. error ',maxval(abs(CI_Vector(:)-CI_SGUGA(:)))
      TraPsi_Verified = .true.
   End If

   Do iCSF = 1, nCSF
      If (Abs(CI_Vector(iCSF)-CI_SGUGA(iCSF))/nCSF>1.0E-12_wp) Then
         Write (6,*) 'Mk_Tra_Psi(SGUGA) error, iCSF=',iCSF
         Write (6,*) 'CI_Vector(iCSF) =',CI_Vector(iCSF)
         Write (6,*) 'CI_SGUGA(iCSF)=',CI_SGUGA(iCSF)
         Call Abend()
      End If
   End Do
   Call mma_deallocate(CI_SGUGA)
End If

End subroutine Mk_Tra_Psi

End module CI_Interfaces
