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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine GrdCls(IRETURN,nState,UEFF,U0,H0)

use PrintLevel, only: VERBOSE
use caspt2_global, only: CLag, CLagFull, do_csf, do_nac, DPT2_AO_tot, DPT2_tot, DPT2C_AO_tot, DPT2C_tot, DPT2Canti_tot, FIFA_all, &
                         FIFASA_all, FIMO_all, iPrGlb, iRoot1, iRoot2, iTasks_grad, LuAPT2, LUGRAD, LuPT2, LUSTD, nOLag, nWLag, &
                         OLag, OLagFull, OMGDER, SLag, TraFro, WLag
use caspt2_module, only: DWTYPE, IfChol, IFDW, IFMSCOUP, IFRMS, IFXMS, iRlxRoot, NBAST, NBSQT, NBTRI, NCONF, nFroT, NROOTS, &
                         REFENE, RFPERT, Zeta
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, King
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IRETURN, nState
real(kind=wp), intent(inout) :: UEFF(nState,nState), U0(nState,nState)
real(kind=wp), intent(in) :: H0(nState,nState)
integer(kind=iwp) :: i, iGo, ilStat, j, jlStat
real(kind=wp) :: CPE, CPTF0, CPTF10, CPUT, Scal, TIOE, TIOTF0, TIOTF10, WALLT
logical(kind=iwp) :: DEB, Found
character(len=16) :: mstate1
real(kind=wp), allocatable :: CI1(:,:), HEFF1(:,:), WRK1(:,:), WRK2(:,:)
integer(kind=iwp), external :: isFreeUnit
logical(kind=iwp), external :: RF_On

!! In case convergence of CASPT2 equation failed
!! Call this subroutine just deallocate memory
if (IRETURN == 0) then

  call mma_allocate(HEFF1,nState,nState,Label='HEFF1')
  call mma_allocate(WRK1,nState,nState,Label='WRK1')
  call mma_allocate(WRK2,nState,nState,Label='WRK2')

  !! Add XMS specific terms
  !! Note that CLagFull is in natural CSF basis,
  !! so everything in this subroutine has to be done in natural
  SLag(:,:) = Zero
  WRK2(:,:) = Zero
  if (IFDW .and. (zeta >= Zero)) then
    !! Construct Heff[1] in XMS basis
    HEFF1(:,:) = Zero
    do ilStat=1,nState
      HEFF1(ilStat,ilStat) = REFENE(ilStat)
    end do
    call DGEMM_('T','N',nState,nState,nState,One,U0,nState,HEFF1,nState,Zero,WRK1,nState)
    call DGEMM_('N','N',nState,nState,nState,One,WRK1,nState,U0,nState,Zero,HEFF1,nState)

    !! Derivative of Heff[1] in XMS basis
    !! It is transformed with U0, so the contribution has to be
    !! considered when we construct the auxiliary density in the
    !! XMS-specific term
    call DWDER(OMGDER,HEFF1,SLag)
    call DGEMM_('N','N',nState,nState,nState,One,U0,nState,SLag,nState,Zero,WRK2,nState)
    call DGEMM_('N','T',nState,nState,nState,One,WRK2,nState,U0,nState,Zero,WRK1,nState)

    WRK2(:,:) = Zero
    do ilStat=1,nState
      if (DWTYPE == 1) then
        WRK2(ilStat,ilStat) = SLag(ilStat,ilStat)
      else if ((DWTYPE == 2) .or. (DWTYPE == 3)) then
        WRK2(ilStat,:) = SLag(ilStat,:)
      end if
      if (.not. do_nac) then
        WRK1(ilStat,1:ilStat-1) = WRK1(ilStat,1:ilStat-1)+WRK1(1:ilStat-1,ilStat)
        WRK1(1:ilStat-1,ilStat) = Zero
      end if
    end do
    SLag(:,:) = WRK1(:,:)
  end if

  if (IFXMS .or. IFRMS .or. (IFMSCOUP .and. do_nac .and. do_csf)) then
    if (.not. IFXMS .and. (.not. IFRMS)) then
      !! For MS-CASPT2, only the second term in eq.(68)
      call unitmat(U0,nState)
    end if

    call TIMING(CPTF0,CPE,TIOTF0,TIOE)
    call XMS_Grad(H0,U0,UEFF,WRK2)
    call TIMING(CPTF10,CPE,TIOTF10,TIOE)
    CPUT = CPTF10-CPTF0
    WALLT = TIOTF10-TIOTF0
    if (IPRGLB >= VERBOSE) write(u6,'(a,2f10.2)') ' XMS_Grad: CPU/WALL TIME=',cput,wallt
  end if

  !! Now, compute the state Lagrangian and do some projections
  !! If PCM, we have to obtain internal state rotation parameters
  !! self-consistently, so we should not call this subroutine
  if (.not. RF_On()) call CLagFinal(nConf,nRoots,nState,CLagFull,SLag)

  !! Add MS-CASPT2 contributions
  if (IFMSCOUP) then
    do ilStat=1,nState
      do jlStat=1,ilStat
        if (do_nac) then
          if ((.not. IFXMS) .and. (.not. IFRMS) .and. (ilstat /= jlstat)) cycle

          Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)+UEFF(jlStat,iRoot1)*UEFF(ilStat,iRoot2)
          Scal = Scal*Half
          SLag(ilStat,jlStat) = SLag(ilStat,jlStat)+Scal

          if (ilStat /= jlStat) SLag(jlStat,ilStat) = SLag(jlStat,ilStat)+Scal
        else
          if (IFXMS .or. IFRMS) then
            Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
            if (ilStat /= jlStat) Scal = Scal*Two
            SLag(ilStat,jlStat) = SLag(ilStat,jlStat)+Scal
          else if (ilStat == jlStat) then
            Scal = UEFF(ilStat,iRoot1)*UEFF(jlStat,iRoot2)
            SLag(ilStat,jlStat) = SLag(ilStat,jlStat)+Scal
          end if
        end if
      end do
    end do
  end if

  !! Subtract the original rhs_sa or rhs_nac contribution
  !! For MS-type CASPT2, CASSCF part has to be determined by UEFF
  if (IFMSCOUP .and. (iRoot1 == iRoot2)) then
    ilStat = max(iRoot1,iRoot2)
    jlStat = min(iRoot1,iRoot2)
    SLag(ilStat,jlStat) = SLag(ilStat,jlStat)-One
  end if

  !! Finalize the first-order transition(-like) density matrix
  !! for the CSF derivative term
  if (do_nac) then
    if (do_csf) then
      call CnstAntiC(DPT2Canti_tot,UEFF,U0)
    else
      !! Clear just in case
      DPT2Canti_tot(:) = Zero
    end if
  end if

  !! Back-transform the CI Lagrangian
  !! It is in the XMS basis, so it has to be transformed to
  !! CASSCF basis to be used in Z-vector
  !! No need to do this for SLag.
  if (IFXMS .or. IFRMS) then
    call mma_allocate(CI1,nConf,nState,Label='CI1')
    call DGEMM_('N','T',nConf,nState,nState,One,CLagFull,nConf,U0,nState,Zero,CI1,nConf)
    CLagFull(:,:) = CI1(:,:)
    call mma_deallocate(CI1)
  end if

  !! Compute true unrelaxed properties for MS-CASPT2
  if ((.not. do_nac) .and. ifmscoup) call PRPCTL(1,UEFF,U0,nState)

  LuPT2 = isFreeUnit(LuPT2)
  call Molcas_Open(LuPT2,'PT2_Lag')

  DEB = .false.
  !! configuration Lagrangian (read in RHS_PT2)
  if (DEB) call RecPrt('CLagFull','',CLagFull,nConf,nState)
  do j=1,nState
    do i=1,nConf
      write(LuPT2,*) CLagFull(i,j)
    end do
  end do

  !! orbital Lagrangian (read in RHS_PT2)
# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    if (.not. King()) OLagFull(:) = Zero
    call GADGOP(OLagFull,nOLag,'+')
  end if
# endif
  if (DEB) call RecPrt('OLagFull','',OLagFull,nBasT,nBasT)
  do i=1,nOLag
    write(LuPT2,*) OLagFull(i)
  end do

        !! state Lagrangian (read in RHS_PT2)
  if (DEB) call RecPrt('SLag','',SLag,nState,nState)
  do j=1,nState
    do i=1,nState
      write(LuPT2,*) SLag(i,j)
    end do
  end do

  !! renormalization contributions (read in OUT_PT2)
# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    if (.not. King()) WLag(:) = Zero
    call GADGOP(WLag,nWLag,'+')
  end if
# endif
  if (DEB) call TriPrt('WLag','',WLag,nBast)
  do i=1,nWLag ! = NBTRI
    write(LuPT2,*) WLag(i)
  end do

  !! D^PT2 in MO (read in OUT_PT2)
# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    if (.not. King()) DPT2_tot(:) = Zero
    call GADGOP(DPT2_tot,NBSQT,'+')
  end if
# endif
  if (DEB) call RecPrt('DPT2','',DPT2_tot,nBast,nBast)
  do i=1,NBSQT
    write(LuPT2,*) DPT2_tot(i)
  end do

  !! D^PT2(C) in MO (read in OUT_PT2)
# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    if (.not. King()) DPT2C_tot(:) = Zero
    call GADGOP(DPT2C_tot,NBSQT,'+')
  end if
# endif
  if (DEB) call RecPrt('DPT2C','',DPT2C_tot,nBast,nBast)
  do i=1,NBSQT
    write(LuPT2,*) DPT2C_tot(i)
  end do

  !! NAC
  if (do_nac) then
    do i=1,NBSQT
      write(LuPT2,*) DPT2Canti_tot(i)
    end do
  end if

  !! D^PT2 in AO (not used?)
  !! Recompute DPT2AO and DPT2AO for PCM
  if (RFpert .and. IFMSCOUP) call Recompute_DPT2AO(NBSQT,DPT2_tot,DPT2C_tot,DPT2_AO_tot,DPT2C_AO_tot)
  if (DEB) call TriPrt('DPT2_AO_tot','',DPT2_AO_tot,nBast)
  do i=1,NBSQT
    write(LuPT2,*) DPT2_AO_tot(i)
  end do

  !! D^PT2(C) in AO (not used?)
  if (DEB) call TriPrt('DPT2C_AO','',DPT2C_AO_tot,nBast)
  do i=1,NBSQT
    write(LuPT2,*) DPT2C_AO_tot(i)
  end do

  if (RFpert) then
    !! For CASPT2/PCM gradient
    !! The implicit derivative contributions have not been
    !! considered in the CASPT2 module
    if (ifmscoup) then
      ! I do not remember why this should be halved
      DPT2_AO_tot(1:NBTRI) = DPT2_AO_tot(1:NBTRI)+Half*DPT2C_AO_tot(1:NBTRI)
    else
      DPT2_AO_tot(1:NBTRI) = DPT2_AO_tot(1:NBTRI)+DPT2C_AO_tot(1:NBTRI)
    end if
    call Put_dArray('D1aoVar',DPT2_AO_tot,NBTRI)
  else
    !! not sure this is OK
    call Put_dArray('D1aoVar',DPT2_AO_tot,0)
  end if

  ! close gradient files
  close(LuPT2)

  call mma_deallocate(HEFF1)
  call mma_deallocate(WRK1)
  call mma_deallocate(WRK2)
end if

call mma_deallocate(DPT2_tot)
call mma_deallocate(DPT2C_tot)
call mma_deallocate(DPT2_AO_tot)
call mma_deallocate(DPT2C_AO_tot)

call mma_deallocate(CLag)
call mma_deallocate(CLagFull)
call mma_deallocate(OLag)
call mma_deallocate(OLagFull)
call mma_deallocate(SLag)
call mma_deallocate(WLag)

call mma_deallocate(FIMO_all)
call mma_deallocate(FIFA_all)

if (IFXMS .or. IFRMS) call mma_deallocate(FIFASA_all)
if (IFDW .and. (zeta >= Zero)) call mma_deallocate(OMGDER)
if (do_nac) call mma_deallocate(DPT2Canti_tot)

 ! Prepare for MCLR
iGo = 3
call Put_iScalar('SA ready',iGo)
! overwrites whatever was set in CASSCF with the relax
! root that was chosen in CASPT2
if (do_nac) then
  !write(u6,*) 'NAC'
  !write(u6,*) 'CASSCF/Original = ',iRoot1,iRoot2
  call Put_iScalar('Relax CASSCF root',iRoot1)
  call Put_iScalar('Relax Original root',iRoot2)
  call Qpg_cArray('MCLR Root',Found,I)
  if (Found) then
    call Get_cArray('MCLR Root',mstate1,16)
    if ((mstate1(8:8) == '0') .and. (mstate1(16:16) == '0')) then
      !! NAC states have not been specified
      write(mstate1,'(1X,I7,1X,I7)') iRoot1,iRoot2
      call Put_cArray('MCLR Root',mstate1,16)
    end if
  else
    mstate1 = '****************'
    call Put_cArray('MCLR Root',mstate1,16)
  end if
else
  !write(u6,*) 'GRD'
  !write(u6,*) 'CASSCF/Original = ',irlxroot,irlxroot
  call Put_iScalar('Relax CASSCF root',irlxroot)
  call Put_iScalar('Relax Original root',irlxroot)
  mstate1 = '****************'
  call Put_cArray('MCLR Root',mstate1,16)
end if

call ModDip()

!! Close files
call DaClos(LUSTD)
if (IfChol) call DaClos(LUAPT2)
call DaClos(LUGRAD)

if (nFroT /= 0) call mma_deallocate(TraFro)
call mma_deallocate(iTasks_grad)

return

end subroutine GrdCls
