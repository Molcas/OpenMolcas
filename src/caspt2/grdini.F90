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

subroutine GrdIni()

use Index_Functions, only: nTri_Elem
use general_data, only: NASH
use caspt2_global, only: CLag, CLagFull, do_lindep, do_nac, DPT2_AO_tot, DPT2_tot, DPT2C_AO_tot, DPT2C_tot, DPT2Canti_tot, &
                         FIFA_all, FIFASA_all, FIMO_all, idBoriMat, idSDMat, iStpGrd, LuAPT2, LuCMOPT2, LuGAMMA, &
                         LUGRAD, LuPT2, LUSTD, nCLag, nOLag, nWLag, OLag, OLagFull, OMGDER, SLag, TraFro, WLag
use caspt2_module, only: IfChol, IFDW, IFRMS, IFXMS, MAXIT, NASUP, NBAS, NBSQT, NBTRI, NCONF, NFROT, NISH, NSTATE, &
                         NSYM, ZETA
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iCase, idSD, idSD_, idSDer, iost, iSym, LENGTH, lRealName, MaxLen, nAS, NS
logical(kind=iwp) :: Exists, is_error
character(len=4096) :: RealName
character(len=128) :: FileName
real(kind=wp), allocatable :: WRK(:)
integer(kind=iwp), external :: isFreeUnit

iStpGrd = 1

! Define (initial) logical unit numbers for gradients files
LUPT2 = 17 ! MCLR
LUGAMMA = 65 ! ERI derivatives ALASKA or 3-center RI/CD
LUCMOPT2 = 66 ! Back-transform ALASKA
LUSTD = 67 ! S and T derivatives in CASPT2
LUAPT2 = 68 ! A_PT2, 2-center derivatives for RI/CD
LUGRAD = 69 ! CASPT2 gradient and property

call PrgmTranslate('GAMMA',RealName,lRealName)
if (IfChol) then
  LENGTH = nBas(1)
else
  LENGTH = nIsh(1)+nAsh(1)
end if
LuGAMMA = isFreeUnit(LuGAMMA)
call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,LENGTH**2*8,'REPLACE',is_error)
close(LuGAMMA)

if (.not. IfChol) then
  call PrgmTranslate('CMOPT2',RealName,lRealName)
  LuCMOPT2 = isFreeUnit(LuCMOPT2)
  call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.false.,1,'REPLACE',is_error)
  close(LuCMOPT2)
end if

call DANAME_MF_wa(LUSTD,'LUSTD')
if (IfChol) call DANAME_MF_wa(LUAPT2,'LUAPT2')

!! Check if this is not the first (MS-)CASPT2 call
!! If this is the first call, compute CASPT2 energies;
!! otherwise read many things from the PT2GRD file.
!! PT2GRD file is always deleted when a RASSCF is finished
!! This is written in Driver/rasscf.prgm.src, but I'm not sure
!! this is the right way to manipulate files...
FileName = 'PT2GRD'
call f_inquire(FileName,Exists)
if (Exists) iStpGrd = 0
! Not sure none/mf/mf_wa
call DANAME_MF_wa(LUGRAD,'LUPT2GRD')

!! Allocate lagrangian terms
! CLag and SLag should allocate for nRoots and not nState,
! but for the time being we only support the case nState=nRoots
nCLag = nconf*nState
nOLag = NBSQT
nWLag = NBTRI

call mma_allocate(DPT2_tot,NBSQT,Label='DPT2_tot')
call mma_allocate(DPT2C_tot,NBSQT,Label='DPT2C_tot')
call mma_allocate(DPT2_AO_tot,NBSQT,Label='DPT2_AO_tot')
call mma_allocate(DPT2C_AO_tot,NBSQT,Label='DPT2C_AO_tot')
DPT2_tot(:) = Zero
DPT2C_tot(:) = Zero
DPT2_AO_tot(:) = Zero
DPT2C_AO_tot(:) = Zero

!! Some Lagrangians for each state are constructed in CLag or
!! OLag. The full (sum over all states, in particular for
!! MS-CASPT2) configuration and orbital Lagrangians are then
!! constructed in OLagFull. For SS-CASPT2, CLag and
!! CLagFull, for instance, will be identical.
call mma_allocate(CLag,nconf,nState,Label='CLAG')
call mma_allocate(CLagFull,nconf,nState,Label='CLAGFULL')
call mma_allocate(OLag,nOLag,Label='OLAG')
call mma_allocate(OLagFull,nOLag,Label='OLAGFULL')
call mma_allocate(SLag,nState,nState,Label='SLAG')
call mma_allocate(WLag,nWLag,Label='WLAG')
CLag(:,:) = Zero
CLagFull(:,:) = Zero
OLag(:) = Zero
OLagFull(:) = Zero
SLag(:,:) = Zero
WLag(:) = Zero
!write(u6,*) 'nclag,nolag'
!write(u6,*) nclag,nolag

call mma_allocate(FIMO_all,NBSQT,Label='FIMO_all')
call mma_allocate(FIFA_all,NBSQT,Label='FIFA_all')
FIMO_all(:) = Zero
FIFA_all(:) = Zero

!! FIFASA is constructed with state-averaged density always
!! FIFA   can be state-specific or dynamically weighted
!! FIMO   is uniquely determined, but the basis can be
!!        either natural or quasi-canonical
if (IFXMS .or. IFRMS) then
  call mma_allocate(FIFASA_all,NBSQT,Label='FIFASA_all')
  FIFASA_all(:) = Zero
  !norbi = norb(1)
end if

if (IFDW .and. (zeta >= Zero)) then
  call mma_allocate(OMGDER,nState,nState,Label='OMGDER')
  OMGDER(:,:) = Zero
end if

if (do_nac) then
  call mma_allocate(DPT2Canti_tot,NBSQT,Label='DPT2Canti_tot')
  DPT2Canti_tot(:) = Zero
end if

MaxLen = maxval(nASUP(1:nSym,1:11)**2)

call mma_allocate(WRK,MaxLen,Label='WRK')
WRK(:) = Zero

idSD = 1
!write (u6,*) 'iflindep = ',iflindeplag
if (do_lindep) then
  do iCase=1,11
    do iSym=1,nSym
      idBoriMat(iSym,iCase) = idSD
      NAS = NASUP(ISYM,ICASE)
      NS = nTri_Elem(NAS)
      call DDAFILE(LuSTD,0,WRK,NS,idSD)
      idSD_ = idBoriMat(iSym,iCase)
      call DDAFILE(LuSTD,1,WRK,NS,idSD_)
    end do
  end do
end if

if (MAXIT /= 0) then
  do iCase=1,11
    do iSym=1,nSym
      idSDMat(iSym,iCase) = idSD
      nAS = nASUP(iSym,iCase)
      call DDAFILE(LuSTD,0,WRK,nAS*nAS,idSD)
      idSDer = idSDMat(iSym,iCase)
      ! idSDMat(iSym,iCase))
      call DDAFILE(LuSTD,1,WRK,nAS*nAS,idSDer)
    end do
  end do
end if
call mma_deallocate(WRK)

if (nFroT /= 0) call mma_allocate(TraFro,nFroT**2,Label='TraFro')
!! iTasks_grad is allocated in MKFG3 where its exact size is determined
!! and in SavGradParams (mode=2) when it is restored from disk.
!call mma_allocate(iTasks_grad,nTri_Elem(ntri2)+ntri1**2,Label='Tasks_grad')
!iTasks_grad(:) = 0

end subroutine GrdIni
