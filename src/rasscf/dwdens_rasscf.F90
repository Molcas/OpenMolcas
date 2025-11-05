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
! Copyright (C) 2025, Yoshio Nishimoto                                 *
!***********************************************************************

! Construct a dynamically-weighted density which is used for reaction field

Subroutine DWDens_RASSCF(CMO,D1A,RCT_FS,IFINAL)

  use Constants, only: Zero, One
  use rasscf_global, only: DoDMRG, ITER, NAC, NACPAR, NACPR2, nRoots, IADR15, Ener
  use definitions, only: iwp,wp
  use DWSol, only: DWSol_wgt, W_SOLV
  use gas_data, only: iDoGAS
  use general_data, only: JOBIPH, NACTEL, NCONF
  use gugx, only: SGS
  use lucia_data, only: PAtmp, Pscr, PTmp, DStmp, Dtmp
  use Lucia_Interface, only: Lucia_Util
  use stdalloc, only: mma_allocate, mma_deallocate
      use sxci, only: IDXSX

#ifdef _DMRG_
  use lucia_data, only: RF1, RF2
  use rasscf_global, only: TwoRDM_qcm
#endif

  implicit none

  real(kind=wp), intent(in) :: CMO(*)
  real(kind=wp), intent(inout) :: D1A(*), RCT_FS(*)
  integer(kind=iwp), intent(in) :: IFINAL

  integer(kind=iwp) :: i, iDisk, iOpt, ITERcurr, jDisk
  real(kind=wp) :: rdum(1), wgt

  real(kind=wp), allocatable :: CIVEC(:), DA_ave(:), DS_ave(:), DX(:)

  Call mma_allocate(DA_ave,NAC**2,Label='DA_ave')
  Call mma_allocate(DS_ave,NAC**2,Label='DS_ave')
  Call mma_allocate(DX,NACPAR,Label='DX')
  DA_ave(1:NACPAR) = Zero
  DS_ave(1:NACPAR) = Zero

  ITERcurr = 1
  if (ITER.ne.1) ITERcurr = ITER-1
  call DWSol_wgt(2,ENER(:,ITERcurr))

  if (iFinal.eq.0 .or. iFinal.eq.1) then
    jDisk = IADR15(3)
    Do i=1,nRoots
      wgt = W_SOLV(i)
      Call DDaFile(JOBIPH,2,DX    ,NACPAR,jDisk)
      Call DDaFile(JOBIPH,2,RCT_FS,NACPAR,jDisk)
      Call DDaFile(JOBIPH,0,rdum,NACPR2,jDisk)
      Call DDaFile(JOBIPH,0,rdum,NACPR2,jDisk)
      if (wgt < 1.0e-10_wp) cycle
      DA_ave(1:NACPAR) = DA_ave(1:NACPAR) + wgt*DX(1:NACPAR)
      DS_ave(1:NACPAR) = DS_ave(1:NACPAR) + wgt*RCT_FS(1:NACPAR)
    End Do
  else if (iFinal.eq.2) then
    Call mma_allocate(CIVEC,NCONF,Label='CIVEC')
    Call mma_allocate(Dtmp,NAC**2,Label='Dtmp')
    Call mma_allocate(DStmp,NAC**2,Label='DStmp')
    Call mma_allocate(Ptmp,NACPR2,Label='Ptmp')

    iDisk = IADR15(4)
    do i = 1, nRoots
      wgt = W_SOLV(i)
      If (NACTEL.EQ.0) THEN
        CIVEC(1)=One
      Else
        if(.not.(doDMRG))then
          iOpt=2
          ! load back one CI vector at the time
          Call DDafile(JOBIPH,iOpt,CIVEC,nConf,iDisk)
        End If
      end if

! compute density matrices

      If ( NAC.ge.1 ) Then
        If (NACTEL.eq.0) THEN
           Dtmp(:)=Zero
           DStmp(:)=Zero
           Ptmp(:)=Zero
        Else
          if(doDMRG)then
#ifdef _DMRG_
            ! copy the DMs from d1rf/d2rf for ipcmroot
            Dtmp(1:NACPAR) = rf1(1:NACPAR)
            if (twordm_qcm) then
              Ptmp(1:NACPR2) = rf2(1:NACPR2)
            end if
            DStmp(:)=Zero
#endif
          else
            Call mma_allocate(PAtmp,NACPR2,Label='PAtmp')
            Call mma_allocate(Pscr,NACPR2,Label='Pscr')
            CALL Lucia_Util('Densi',CI_Vector=CIVEC(:))
            If (SGS%IFRAS.GT.2 .OR. iDoGAS) Then
              Call CISX(IDXSX,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
            End If
            Call mma_deallocate(Pscr)
            Call mma_deallocate(PAtmp)
          end if ! doDMRG/doBLOK or CI
        End If
      Else
        Dtmp(:)=Zero
        DStmp(:)=Zero
        Ptmp(:)=Zero
      End If
      DA_ave(1:NACPAR) = DA_ave(1:NACPAR) + wgt*Dtmp(1:NACPAR)
      DS_ave(1:NACPAR) = DS_ave(1:NACPAR) + wgt*DStmp(1:NACPAR)
    End Do
    Call mma_deallocate(DStmp)
    Call mma_deallocate(Dtmp)
    Call mma_deallocate(Ptmp)
    Call mma_deallocate(CIVEC)
  End If
!
! Construc D-ACTIVE AND D-INACTIVE IN AO BASIS
!
  DX(1:NACPAR) = DS_ave(1:NACPAR)
  Call DBLOCK(DX)
  CALL Get_D1A_RASSCF(CMO,DX,RCT_FS)

  DX(1:NACPAR) = DA_ave(1:NACPAR)
  Call DBLOCK(DX)
  CALL Get_D1A_RASSCF(CMO,DX,D1A)

  Call mma_deallocate(DA_ave)
  Call mma_deallocate(DS_ave)
  Call mma_deallocate(DX)

  return

End Subroutine DWDens_RASSCF
