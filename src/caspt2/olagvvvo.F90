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

subroutine OLagVVVO(iSym,NBSQT,lT2AO,MaxVec_PT2,DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,T2AO,DIA,DI,FIFA,FIMO,A_PT2)

use iSD_data, only: iSD
use caspt2_global, only: LuAPT2, LuCMOPT2, LuGAMMA, OLag
use caspt2_global, only: CMOPT2
use caspt2_module, only: IFDW, IFMSCOUP, IFRMS, IFSADREF, IFXMS, iRlxRoot, JSTATE, NASH, NBAS, NBAST, NBMX, NFRO, NISH, NSSH, &
                         NSTATE, NSYM
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
use caspt2_global, only: nOLag
use caspt2_module, only: NFROT
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iSym, NBSQT, lT2AO, MaxVec_PT2
real(kind=wp), intent(in) :: DPT2AO(NBSQT), DPT2CAO(NBSQT), T2AO(lT2AO), DIA(NBSQT), DI(NBSQT)
real(kind=wp), intent(inout) :: FPT2AO(NBSQT), FPT2CAO(NBSQT), FIFA(NBSQT), FIMO(NBSQT), A_PT2(MaxVec_PT2**2)
integer(kind=iwp) :: i, iBas, iBas0, id, iOcc, iost, IRC, iRec, iSh, iSymA, iSymB, iSymI, iSymJ, j, jBas, jBas0, jOcc, jSh, &
                     KEEP(8), loc1, loc2, lRealName, MaxShlAO, nBasI, nBasJ, nBasX(8), nDiff, nocc, nOrbA, nSkal, nSSDM, nSymX
logical(kind=iwp) :: DoCholesky, DoRys, is_error, Square
character(len=4096) :: RealName
integer(kind=iwp), allocatable :: iOffAO(:)
real(kind=wp), allocatable :: T_hbf(:,:,:,:), vLag(:), WRK1(:), WRK2(:)
integer(kind=iwp), external :: isFreeUnit

! ----- (VV|VO)

! Compute L_{pq} = (pa|jb) * T_{qj}^{ab}, in particular for
! (p,q) \in  (virt, inact+act). This operation involves
! (VV|VO) integrals which are stored in neither disk nor memory.
! The back-transformation is applied to occupied orbital derivatives
! of two-electron integrals that have two virtual indices (F, G, H).
! In principle, the algorithm is to avoid (VV|VO) integrals,
! i.e. U_{pq} for (p,q) = (vir, inact+act), but can also be applied
! to U_{pq} for (p,q) = (all, inact+act).

call mma_allocate(vLag,nBasT*nBasT,Label='vLag')
call mma_allocate(WRK1,nBasT*nBasT,Label='WRK1')
call GetOrd(IRC,Square,nSymX,nBasX,KEEP)

vLag(:) = Zero

call DecideOncholesky(DoCholesky)
if (DoCholesky) then
  !! No need to save CMOPT2. Just save A_PT2 and B_PT2.
  !! First, save A_PT2 in LuCMOPT2
  if (IFMSCOUP .and. (jState /= 1)) then
    call mma_allocate(WRK2,MaxVec_PT2**2,Label='WRK2')

    ! read A_PT2 from LUAPT2
    id = 0
    call ddafile(LUAPT2,2,WRK2,MaxVec_PT2**2,id)

    A_PT2(1:MaxVec_PT2**2) = A_PT2(1:MaxVec_PT2**2)+WRK2(1:MaxVec_PT2**2)
    call mma_deallocate(WRK2)
  end if

  ! For SS-CASPT2 I should write A_PT2 on disk only
  ! for the correct iRlxRoot
  if ((jState == iRlxRoot) .or. IFMSCOUP) then
    ! write A_PT2 in LUAPT2
    id = 0
    call ddafile(LUAPT2,1,A_PT2,MaxVec_PT2**2,id)
  end if

  !rewind(LuGamma)
  call PrgmTranslate('GAMMA',RealName,lRealName)
  LuGAMMA = isFreeUnit(LuGAMMA)
  call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,nBas(iSym)**2*8,'OLD',is_error)
end if
!! 2) Compute ERI (mu rho | nu sigma)
!! 3) Quarter-transformation of ERI
!!    (mu rho | j sigma) = sum_{j} C_{nu j} (mu rho | nu sigma)
!! 4) Contract with AO amplitude
!!    L_{mu i} = sum_{j,rho sigma} T_{ij}^{mu nu}*(mu rho|j sigma)
!! the third argument is used wrongly, but should be no problem

!! D_{mu nu} := DPT2AO
!! FPT2AO = ((mu nu|rho sigma)-(mu rho|nu sigma)/4
!!          -(mu sigma|nu rho)/4)*D_{rho sigma}
isymi = 1
isymj = 1
isyma = 1
isymb = 1
!nocc = nfro(1)+nish(1)+nash(1)
nocc = nish(1)+nash(1)
call VVVO_Drv(nSym,nBas,nIsh,nFro,KEEP,iSym,iSymI,iSymA,iSymJ,iSymB,lT2AO,T2AO,vLag,nOcc,nBasT,nBMX,CMOPT2(1+nBasT*nFro(iSymA)), &
              DPT2AO,DPT2CAO,FPT2AO,FPT2CAO,DIA,DI,FIFA,FIMO)

!! Save the half transformed integrals on disk.
!! It will be used in drvg1 etc for gradient.
if (DoCholesky) then
  !! Do nothing
  close(LuGAMMA)
else
  !! This is only for conventional calculations!!
  call PrgmTranslate('CMOPT2',RealName,lRealName)
  LuCMOPT2 = isFreeUnit(LuCMOPT2)
  call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.false.,1,'OLD',is_error)
  !! First, CMOPT2 has to be saved. The MO coefficient matrix in
  !! grvg1 may be different from CMOPT2.
  do i=1,nBasT*nBasT
    write(LuCMOPT2) CMOPT2(i)
  end do
  do i=1,8
    write(LuCMOPT2) nIsh(i)+nAsh(i)
  end do
  do i=1,8
    write(LuCMOPT2) nFro(i)
  end do
  !! Number of state-specific density matrix
  !! It is needed, because the separable contribution to the
  !! electron repulsion integral is different
  nSSDM = 0
  if (.not. IFSADREF) then
    if (nState == 1) then
      nSSDM = 0
    else if (IFDW .or. IFRMS) then
      !! For (X)DW, use nState density matrix
      nSSDM = nState
    else if (IFXMS) then
      !! For XMS, use SA density matrix
      nSSDM = 0
    else if (IFMSCOUP) then
      !! For MS, use nState density matrix
      nSSDM = nState
    else
      !! Otherwise, use one SS density matrix
      nSSDM = 1
    end if
  end if
  write(LuCMOPT2) nSSDM

  close(LuCMOPT2)

  !write(u6,*) 'mo saved'
  !call sqprt(CMOPT2,nbast)

  call PrgmTranslate('GAMMA',RealName,lRealName)
  LuGAMMA = isFreeUnit(LuGAMMA)
  call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,nOcc*nOcc*8,'OLD',is_error)
  if (is_error) then
    write(u6,*) 'Something is wrong in opening LuGamma in olagvvvo'
    call abend()
  end if
  ! Setup for shell. Why do I have to call IniSew damn here?
  ! The number of shells should be able to be referred globally.
  nDiff = 1
  DoRys = .true.
  call IniSew(DoRys,nDiff)
  call Nr_Shells(nSkal)
  call Setup_iSD()
  !! see Include/info.fh
  call mma_allocate(iOffAO,nSkal+1,Label='iOffAO')
  MaxShlAO = 0
  iOffAO(1) = 0
  do iSh=1,nSkal
    nBasI = iSD(2,iSh)*iSD(3,iSh)
    if (nBasI > MaxShlAO) MaxShlAO = nBasI
    iOffAO(iSh+1) = iOffAO(iSh)+nBasI
  end do
  !nMax = MaxShlAO*MaxShlAO*nOcc*nOcc

  call mma_allocate(T_hbf,nOcc,nOcc,MaxShlAO,MaxShlAO,Label='T_hbf')
  do iSh=1,nSkal
    !! iSD(2,iSh): number of AOs of the shell
    !! iSD(3,iSh): number of cont. func. of the shell
    nBasI = iSD(2,iSh)*iSD(3,iSh)
    do jSh=1,nSkal
      nBasJ = iSD(2,jSh)*iSD(3,jSh)
      do iBas0=1,nBasI
        iBas = iOffAO(iSh)+iBas0
        do jBas0=1,nBasJ
          jBas = iOffAO(jSh)+jBas0
          do iOcc=1,nOcc
            do jOcc=1,nOcc
              loc1 = jOcc-1+(jBas-1)*nOcc+(iOcc-1)*nOcc*nBasT+(iBas-1)*nOcc*nBasT*nOcc
              loc2 = iOcc-1+(iBas-1)*nOcc+(jOcc-1)*nOcc*nBasT+(jBas-1)*nOcc*nBasT*nOcc
              T_hbf(iOcc,jOcc,iBas0,jBas0) = T2AO(loc1+1)+T2AO(loc2+1)
            end do
          end do
          iRec = iBas+nBasT*(jBas-1)
          if (ifmscoup .and. (jstate /= 1)) then
            read(lugamma,rec=irec) (wrk1(i),i=1,nocc*nocc)
            call daxpy_(nocc*nocc,One,wrk2,1,t_hbf(1,1,ibas0,jbas0),1)
          end if
          write(LuGamma,rec=iRec) ((T_hbf(i,j,iBas0,jBas0),i=1,nOcc),j=1,nOcc)
        end do
      end do
    end do
  end do
  call mma_deallocate(iOffAO)
  call mma_deallocate(T_hbf)
  close(LuGAMMA)
  call Free_iSD()
  call clssew()
end if

!! 5) L_{ai} = sum_{mu} C_{mu a} L_{mu i}
!call DGEMM_('T','N',nSsh(iSym),nOcc,nBasT,One,CMOPT2(1+nBasT*nOcc),nBasT,vLag,nBasT,One,OLAG(nOCC+1),nOrb(iSymA))
!write(u6,*) 'olag before vvvo'
!call sqprt(olag,nbast)
nOrbA = nFro(iSymA)+nIsh(iSymA)+nAsh(iSymA)+nSsh(iSymA)
if (DoCholesky) nOcc = nOrbA-nFro(iSymA)
!write(u6,*) 'vLag'
!call sqprt(vLag,nbast)
call DGEMM_('T','N',nOrbA,nOcc,nBasT,One,CMOPT2,nBasT,vLag,nBasT,One,OLAG(nOrbA*nFro(iSymA)+1),nOrbA)
!write(u6,*) 'olag after vvvo'
!call sqprt(olag,nbast)

#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) then
  if (DoCholesky) call GADGOP(OLag,nOLag,'+')
  if (nFroT == 0) then
    call GADGOP(FPT2AO,nBasT**2,'+')
    call GADGOP(FPT2CAO,nBasT**2,'+')
  end if
end if
#endif
call mma_deallocate(vLag)
call mma_deallocate(WRK1)

end subroutine OLagVVVO
