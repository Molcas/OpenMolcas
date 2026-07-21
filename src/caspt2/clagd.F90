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

subroutine CLagD(NASHT,NG3,NSTATE,G1,G2,G3,DG1,DG2,DG3,DF1,DF2,DF3,DEASUM,DEPSA,VECROT)

use Index_Functions, only: nTri_Elem
use EQSOLV, only: IDBMAT, IDSMAT, IRHS, IVECR, IVECW, IVECX
use fake_GA, only: GA_Arrays
use caspt2_global, only: imag_shift, ipea_shift, iVecL, LUSBT, LUSOLV, sigma_p_epsilon
use general_data, only: NASH
use caspt2_module, only: EASUM, EPSA, HZERO, IFMSCOUP, NAES, NASUP, NINDEP, NISUP, NSYM, NTGEUES, NTGTUES, NTUES
use BDerNEV, only: BDNA, BDNB, BDNC, BDND, BDNE, BDNF, BDNG, BDN_G3
use SC_NEVPT2, only: SC_prop, SC_NEVPT2_CLagD
#ifdef _MOLCAS_MPP_
use caspt2_global, only: do_lindep, idSDMat, LUSTD, real_shift
use Para_Info, only: Is_Real_Par, King
use GA_Wrapper, only: DBL_MB, GA_Destroy, GA_NodeId
use Definitions, only: u6
use EQSOLV, only: IDTMAT
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Four, Half
use Definitions, only: wp, iwp, byte

implicit none
integer(kind=iwp), intent(in) :: NASHT, NG3, NSTATE
real(kind=wp), intent(inout) :: G1(NASHT,NASHT), G2(NASHT,NASHT,NASHT,NASHT), G3(NG3), DG1(NASHT,NASHT), &
                                DG2(NASHT,NASHT,NASHT,NASHT), DG3(NG3), DF1(NASHT,NASHT), DF2(NASHT,NASHT,NASHT,NASHT), DF3(NG3), &
                                DEASUM, DEPSA(NASHT,NASHT)
real(kind=wp), intent(in) :: VECROT(NSTATE)
integer(kind=iwp) :: iCase, ID, idS, idum, iLUID, iSym, iTabs, iTgeUabs, iTgtUabs, iTU2, iTUabs, iUabs, iXabs, iXgeYabs, &
                     iXgtYabs, iXY2, iXYabs, iYabs, lg_V1, lg_V2, lg_V3, lg_V4, lg_V5, NAS, NIN, NIS, NS, NSEQ, NVEC
real(kind=wp) :: ATUX, ATUXY, ATUY, ATYU, ATYX, BDERval, bsBDER, ET, EU, EX, EY, ScalB1, ScalB2, ScalS1, ScalS2, SDERval
real(kind=wp), allocatable :: BDER(:), LBD(:), LID(:), SDER(:), SMat(:)
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: idT, IHI, ILO, JHI, JLO, LDV, lg_S, lg_T, mS, MYRANK
logical(kind=iwp) :: bStat
integer(kind=byte), allocatable :: idxG3(:,:)
real(kind=wp), allocatable :: VEC1(:), VEC2(:), VEC3(:), VEC4(:), VEC5(:)
#endif

if ((HZERO == 'DYALL') .and. SC_prop) then
  call SC_NEVPT2_CLagD(NASHT,NG3,NSTATE,G1,G2,G3,DG1,DG2,DG3,VECROT)
  return
end if

do iCase=1,11
  !cycle
  !if ((iCase /= 10) .and. (iCase /= 11)) cycle ! G
  !if (iCase /= 10) cycle ! GP
  !if ((iCase /= 6) .and. (iCase /= 7)) cycle ! E
  !if ((iCase /= 8) .and. (iCase /= 9)) cycle ! F
  !if (iCase /= 8) cycle ! FP
  !if ((iCase /= 2) .and. (iCase /= 3)) cycle ! B
  !if (icase /= 5) cycle ! D
  !if (iCase /= 4) cycle ! C
  !if (iCase /= 1) cycle ! A
  do iSym=1,nSym
    nIN = nINDEP(iSym,iCase)
    if (nIN == 0) cycle
    nIS = nISUP(iSym,iCase)
    NVEC = nIN*nIS
    nAS = nASUP(iSym,iCase)
    if (nVec == 0) cycle

#   ifdef _MOLCAS_MPP_
    if (is_real_par()) then
      if (((iCase /= 1) .and. (iCase /= 4)) .or. (HZERO == 'DYALL')) then
        call mma_allocate(BDER,NAS**2,Label='BDER')
        call mma_allocate(SDER,NAS**2,Label='SDER')
        BDER(:) = Zero
        SDER(:) = Zero
      end if
      if ((HZERO == 'DYALL') .and. ((iCase == 1) .or. (iCase == 4))) then
        call GA_CREATE_STRIPED('H',NAS,NIN,'TRANS',lg_T)
        call PSBMAT_READ('T',iCase,iSym,lg_T,NAS*NIN)
        if (King()) then
          call mma_allocate(VEC1,NAS*NIN,Label='VEC1')
          call GA_GET(lg_T,1,NAS,1,NIN,VEC1,NAS)
          idT = idTMAT(iSym,iCase)
          call DDAFILE(LUSBT,1,VEC1,NAS*NIN,idT)
          call mma_deallocate(VEC1)
        end if
        bStat = GA_destroy(lg_T)
      end if
    else
#   endif
      call mma_allocate(BDER,NAS**2,Label='BDER')
      call mma_allocate(SDER,NAS**2,Label='SDER')
      BDER(:) = Zero
      SDER(:) = Zero
#   ifdef _MOLCAS_MPP_
    end if
#   endif

#   if defined(_MOLCAS_MPP_) && defined(_GA_)
    if (is_real_par() .and. ((icase == 1) .or. (icase == 4)) .and. (HZERO /= 'DYALL')) then
      call CLagDX_MPP()
      cycle
    else
#   endif

      !write(u6,*) 'for icase = ',icase
      !write(u6,*) '# of independent vecs:',nin
      !write(u6,*) '# of non-active pairs:',nis
      !write(u6,*) '# of     active pairs:',nas
      !write(u6,*) 'dimension for Vec = ',nin*nis
      !! lg_V1 = T (solution; not quasi-variational)
      call RHS_ALLO(nIN,nIS,lg_V1)
      call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
      !! lg_V2 = lambda (shift correction)
      call RHS_ALLO(nIN,nIS,lg_V2)
      call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
      if (sigma_p_epsilon /= Zero) then
        call mma_allocate(LBD,nAS,Label='LBD')
        call mma_allocate(LID,nIS,Label='LID')
        iD = iDBMat(iSym,iCase)
        call dDaFile(LUSBT,2,LBD,nAS,iD)
        call dDaFile(LUSBT,2,LID,nIS,iD)
        call CASPT2_ResD(3,nIN,nIS,lg_V2,lg_V1,LBD,LID)
        call mma_deallocate(LBD)
        call mma_deallocate(LID)
      end if

      !! lg_V3 = RHS (in IC basis)
      call RHS_ALLO(nIN,nIS,lg_V3)
      call RHS_READ_SR(lg_V3,iCase,iSym,iRHS)
      !! lg_V4 = RHS (in MO basis)
      call RHS_ALLO(nAS,nIS,lg_V4)
      call RHS_READ(nAS,nIS,lg_V4,iCase,iSym,iVecW)
      !! lg_V5 = RHS2 (in IC basis)
      if (IFMSCOUP) then
        call RHS_ALLO(nIN,nIS,lg_V5)
        call RHS_READ_SR(lg_V5,iCase,iSym,iVecL) ! 7
      else
        lg_V5 = lg_V3
      end if

#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        if (KING()) then
          ! copy global array to local buffer
          call mma_allocate(VEC1,NVEC,Label='VEC1')
          call GA_GET(lg_V1,1,NIN,1,NIS,VEC1,NIN)
          call mma_allocate(VEC2,NVEC,Label='VEC2')
          call GA_GET(lg_V2,1,NIN,1,NIS,VEC2,NIN)
          call mma_allocate(VEC3,NIN*NIS,Label='VEC3')
          call GA_GET(lg_V3,1,NIN,1,NIS,VEC3,NIN)
          call mma_allocate(VEC4,NAS*NIS,Label='VEC4')
          call GA_GET(lg_V4,1,NAS,1,NIS,VEC4,NAS)
          !! is it possible to avoid allocating twice?
          !if (IFMSCOUP) then
          !  call mma_allocate(VEC5,NIN*NIS,Label='VEC5')
          !  call GA_GET(lg_V5,1,NIN,1,NIS,VEC5,NIN)
          !else
          !  LVEC5 = LVEC3
          !end if
          call mma_allocate(VEC5,NIN*NIS,Label='VEC5')
          call GA_GET(lg_V5,1,NIN,1,NIS,VEC5,NIN)

          call CLagDX(0,ISYM,ICASE,VEC1,VEC2,VEC3,VEC4,nIN,nIS,nAS,nState,VECROT,VEC5,lg_V2,BDER,SDER)

          ! free local buffer
          call mma_deallocate(VEC1)
          call mma_deallocate(VEC2)
          call mma_deallocate(VEC3)
          call mma_deallocate(VEC4)
          !if (IFMSCOUP) then
          call mma_deallocate(VEC5)
          !end if
        end if
        call GASYNC()
      else
#     endif
        call CLagDX(0,iSym,iCase,GA_Arrays(lg_V1)%A,GA_Arrays(lg_V2)%A,GA_Arrays(lg_V3)%A,GA_Arrays(lg_V4)%A,nIN,nIS,nAS,nState, &
                    VECROT,GA_Arrays(lg_V5)%A,lg_V2,BDER,SDER)
#     ifdef _MOLCAS_MPP_
      end if
#     endif

      if ((imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero)) then
        nAS = nASUP(iSym,iCase)
        call mma_allocate(LBD,nAS,Label='LBD')
        call mma_allocate(LID,nIS,Label='LID')
        iD = iDBMat(iSym,iCase)
        call dDaFile(LUSBT,2,LBD,nAS,iD)
        call dDaFile(LUSBT,2,LID,nIS,iD)

        call RHS_READ_SR(lg_V1,ICASE,ISYM,iVecX)
        call RHS_READ_SR(lg_V2,ICASE,ISYM,iVecR)
        call CASPT2_ResD(2,nIN,nIS,lg_V1,lg_V2,LBD,LID)
        call mma_deallocate(LBD)
        call mma_deallocate(LID)

#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          if (KING()) then
            ! copy global array to local buffer
            call mma_allocate(VEC1,NVEC,Label='VEC1')
            call GA_GET(lg_V1,1,NIN,1,NIS,VEC1,NIN)
            call mma_allocate(VEC2,NVEC,Label='VEC2')
            call GA_GET(lg_V2,1,NIN,1,NIS,VEC2,NIN)

            !! lvec3 to lvec5 are not used
            !! only lvec1 and lvec2?
            call mma_allocate(VEC3,1,Label='VEC3')
            call mma_allocate(VEC4,1,Label='VEC4')
            call mma_allocate(VEC5,1,Label='VEC5')
            call CLagDX(1,ISYM,ICASE,VEC1,VEC2,VEC3,VEC4,nIN,nIS,nAS,nState,VECROT,VEC5,lg_V2,BDER,SDER)

            ! free local buffer
            call mma_deallocate(VEC1)
            call mma_deallocate(VEC2)
            call mma_deallocate(VEC3)
            call mma_deallocate(VEC4)
            call mma_deallocate(VEC5)
          end if
          call GASYNC()
        else
#       endif
          call CLagDX(1,iSym,iCase,GA_Arrays(lg_V1)%A,GA_Arrays(lg_V2)%A,GA_Arrays(lg_V3)%A,GA_Arrays(lg_V4)%A,nIN,nIS,nAS,nState, &
                      VECROT,GA_Arrays(lg_V5)%A,lg_V2,BDER,SDER)
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      end if

      !! for non-separable density/derivative
      call RHS_READ_SR(lg_V1,ICASE,ISYM,iVecX)
      call RHS_READ_SR(lg_V2,ICASE,ISYM,iVecR)
#   if defined(_MOLCAS_MPP_) && defined(_GA_)
    end if
#   endif

    if (HZERO /= 'DYALL') then !! CASPT2
      select case (iCase)
        case (1)
          call CLagDXA(NAS,BDER,SDER)
        case (2)
          call CLagDXB(NAS,BDER,SDER)
        case (3)
          call CLagDXB(NAS,BDER,SDER)
        case (4)
          call CLagDXC(NAS,BDER,SDER)
        case (5)
          call CLagDXD(NAS,BDER,SDER)
        case (6)
          call CLagDXE(NAS,BDER,SDER)
        case (7)
          call CLagDXE(NAS,BDER,SDER)
        case (8)
          call CLagDXF(NAS,BDER,SDER)
        case (9)
          call CLagDXF(NAS,BDER,SDER)
        case (10)
          call CLagDXG(NAS,BDER,SDER)
        case (11)
          call CLagDXG(NAS,BDER,SDER)
      end select
    else !! NEVPT2
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        call GADGOP(BDER,NAS**2,'+')
        call GADGOP(SDER,NAS**2,'+')
      end if
#     endif
      select case (iCase)
        case (1)
          call BDNA(iSym,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)
        case (2,3)
          call BDNB(iSym,iCase,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)
        case (4)
          call BDNC(iSym,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)
        case (5)
          call BDND(iSym,NAS,NG3,BDER,SDER,G1,G2,G3,DG1,DG2,DG3)
        case (6,7)
          call BDNE(iSym,NAS,BDER,SDER,G1,G2,DG1,DG2)
        case (8,9)
          call BDNF(iSym,iCase,NAS,NG3,BDER,SDER,G2,G3,DG2,DG3)
        case (10,11)
          call BDNG(iSym,NAS,BDER,SDER,G1,G2,DG1,DG2)
      end select
    end if

    call RHS_FREE(lg_V1)
    call RHS_FREE(lg_V2)
    call RHS_FREE(lg_V3)
    call RHS_FREE(lg_V4)
    if (IFMSCOUP) call RHS_FREE(lg_V5)

#   ifdef _MOLCAS_MPP_
    if (is_real_par()) then
      if (((iCase /= 1) .and. (iCase /= 4)) .or. (HZERO == 'DYALL')) then
        call mma_deallocate(BDER)
        call mma_deallocate(SDER)
      end if
    else
#   endif
      call mma_deallocate(BDER)
      call mma_deallocate(SDER)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end do
end do

#ifdef _MOLCAS_MPP_
if (is_real_par() .and. (HZERO /= 'DYALL')) then
  iCase = 4
  MYRANK = GA_NODEID()
  do iSym=1,nSym
    nAS = nASUP(iSym,iCase)
    call PSBMAT_GETMEM('S',lg_S,NAS)
    call PSBMAT_READ('S',4,iSym,lg_S,NAS)
    call GA_Distribution(lg_S,myRank,ILO,IHI,JLO,JHI)
    call GA_Access(lg_S,ILO,IHI,JLO,JHI,mS,LDV)
    call mma_allocate(idxG3,6,NG3,label='idxG3')
    iLUID = 0
    call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
    !! DF3 is done after icase=4
    !! construct G3 matrix in lg_S
    call MKSC_G3_MPP(ISYM,DBL_MB(mS),ILO,NAS,LDV,NG3,G3,IDXG3)
    call GA_Release_Update(lg_S,ILO,IHI,JLO,JHI)
    call DF3_DEPSA_MPP(NG3,NASHT,DF3,DEPSA,lg_S,idxG3)

    call GA_Release(lg_S,ILO,IHI,JLO,JHI)
    call mma_deallocate(idxG3)
    call PSBMAT_FREEMEM(lg_S)
  end do
end if
#endif

if (HZERO == 'DYALL') call BDN_G3(DG1,DG2,DG3)

return

contains

subroutine CLagDXA(NAS,BDER,SDER)

  integer(kind=iwp), intent(in) :: NAS
  real(kind=wp), intent(in) :: BDER(NAS,NAS)
  real(kind=wp), intent(inout) :: SDER(NAS,NAS)
  integer(kind=byte), allocatable :: idxG3(:,:)

  NS = nTri_Elem(NAS)
  call mma_allocate(SMat,NS,Label='SMat')
  idS = idSMAT(iSym,1)
  call DDAFILE(LUSBT,2,SMat,NS,idS)

  idum = 0
  call CLagDXA_DP(iSym,nAS,nAshT,BDER,SDER,DG1,DG2,DF1,DF2,DEPSA,DEASUM,1,nAS,1,nAS,0,G1,G2,SMat,SMat,idum)

  !! G3 and F3 relevant
  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID = 0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
  !idS = idSMAT(iSym,4)
  !call DDAFILE(LUSBT,2,SMat,NS,idS)
  call MKSC_G3(iSym,SMat,NS,nG3,G3,idxG3)
  call CLagDXA_FG3(iSym,nAS,nAshT,NG3,NS,BDER,SDER,DF1,DF2,DF3,DG1,DG2,DG3,DEPSA,G2,SMat,idxG3)
  call mma_deallocate(idxG3)

  call mma_deallocate(SMat)

  return

end subroutine CLagDXA

subroutine CLagDXB(NAS,BDER,SDER)

  use SUPERINDEX, only: MTGEU, MTGTU

  integer(kind=iwp), intent(in) :: NAS
  real(kind=wp), intent(in) :: BDER(NAS,NAS)
  real(kind=wp), intent(inout) :: SDER(NAS,NAS)
  integer(kind=iwp) :: iT, iTU, iU, iX, iXY, iY
  real(kind=wp), allocatable :: WrkBbf(:,:,:,:), WrkSbf(:,:,:,:)

  call mma_allocate(WrkBbf,nAshT,nAshT,nAshT,nAshT,Label='WrkBbf')
  call mma_allocate(WrkSbf,nAshT,nAshT,nAshT,nAshT,Label='WrkSbf')
  WrkBbf(:,:,:,:) = Zero
  WrkSbf(:,:,:,:) = Zero

  if (ipea_shift /= Zero) then
    NS = nTri_Elem(NAS)
    call mma_allocate(SMat,NS,Label='SMat')
    idS = idSMAT(iSym,iCase)
    call DDAFILE(LUSBT,2,SMat,NS,idS)
  end if
  ScalB1 = Zero
  ScalB2 = Zero
  ScalS1 = Zero
  ScalS2 = Zero
  iTabs = 0
  iUabs = 0
  iXabs = 0
  iYabs = 0
  iTgeUabs = 0
  iTgtUabs = 0
  iXgeYabs = 0
  iXgtYabs = 0
  do iTU=1,nAS
    if (iCase == 2) then
      iTgeUabs = iTU+nTgeUes(iSym)
      iTabs = mTgeU(1,iTgeUabs)
      iUabs = mTgeU(2,iTgeUabs)
    else if (iCase == 3) then
      iTgtUabs = iTU+nTgtUes(iSym)
      iTabs = mTgtU(1,iTgtUabs)
      iUabs = mTgtU(2,iTgtUabs)
    end if
    ET = EPSA(iTabs)
    EU = EPSA(iUabs)
    do iXY=1,nAS
      if (iCase == 2) then
        iXgeYabs = iXY+nTgeUes(iSym)
        iXabs = mTgeU(1,iXgeYabs)
        iYabs = mTgeU(2,iXgeYabs)
      else if (iCase == 3) then
        iXgtYabs = iXY+nTgtUes(iSym)
        iXabs = mTgtU(1,iXgtYabs)
        iYabs = mTgtU(2,iXgtYabs)
      end if
      EX = EPSA(iXabs)
      EY = EPSA(iYabs)
      ATUXY = EASUM-ET-EU-EX-EY
      !iBadr = iTU + nAS*(iXY-1)
      BDERval = BDER(ITU,IXY)

      !! For IPEA shift
      if ((iTU == iXY) .and. (ipea_shift /= Zero)) then
        !idT = nTri_Elem(iTabs)
        !idU = nTri_Elem(iUabs)
        NSEQ = nTri_Elem(iTU)
        bsBDER = ipea_shift*Half*BDERval
        !! ipea_shift*Half*(DREF(IDT)+DREF(IDU))*SDP(ITGEU)
        DG1(iTabs,iTabs) = DG1(iTabs,iTabs)+bsBDER*SMat(NSEQ)
        DG1(iUabs,iUabs) = DG1(iUabs,iUabs)+bsBDER*SMat(NSEQ)
        SDER(iTU,iXY) = SDER(iTU,iXY)+(G1(iTabs,iTabs)+G1(iUabs,iUabs))*bsBDER
      end if
      SDERval = SDER(ITU,IXY)
      if (iTabs == iUabs) then
        BDERval = BDERval*Two
        SDERval = SDERval*Two
      end if

      if (iCase == 2) then
        ScalB1 = BDERval
        ScalB2 = BDERval
        ScalS1 = SDERval
        ScalS2 = SDERval
      else if (iCase == 3) then
        ScalB1 = BDERval
        ScalB2 = -BDERval
        ScalS1 = SDERval
        ScalS2 = -SDERval
      end if

      WRKBBF(iTabs,iUabs,iXabs,iYabs) = WRKBBF(iTabs,iUabs,iXabs,iYabs)+ScalB1
      WRKBBF(iTabs,iUabs,iYabs,iXabs) = WRKBBF(iTabs,iUabs,iYabs,iXabs)+ScalB2
      WRKSBF(iTabs,iUabs,iXabs,iYabs) = WRKSBF(iTabs,iUabs,iXabs,iYabs)+ScalS1
      WRKSBF(iTabs,iUabs,iYabs,iXabs) = WRKSBF(iTabs,iUabs,iYabs,iXabs)+ScalS2
      if (iTabs /= iUabs) then
        WRKBBF(iUabs,iTabs,iXabs,iYabs) = WRKBBF(iUabs,iTabs,iXabs,iYabs)+ScalB2
        WRKBBF(iUabs,iTabs,iYabs,iXabs) = WRKBBF(iUabs,iTabs,iYabs,iXabs)+ScalB1
        WRKSBF(iUabs,iTabs,iXabs,iYabs) = WRKSBF(iUabs,iTabs,iXabs,iYabs)+ScalS2
        WRKSBF(iUabs,iTabs,iYabs,iXabs) = WRKSBF(iUabs,iTabs,iYabs,iXabs)+ScalS1
      end if
    end do
  end do

  !! it,iu,... are actually itabs,iuabs,...
  do iT=1,nAshT
    ET = EPSA(iT)
    do iU=1,nAshT
      EU = EPSA(iU)
      do iX=1,nAshT
        EX = EPSA(iX)
        do iY=1,nAshT
          EY = EPSA(iY)
          BDERval = WRKBBF(iT,iU,iX,iY)
          SDERval = WRKSBF(iT,iU,iX,iY)

          !! term 1 (w/o delta)
          ATUXY = EASUM-ET-EU-EX-EY
          !! G1 and F1 derivative
          DF2(iX,iT,iY,iU) = DF2(iX,iT,iY,iU)+BDERval
          DG2(iX,iT,iY,iU) = DG2(iX,iT,iY,iU)-ATUXY*BDERval+SDERval
          !! EASUM derivative
          DEASUM = DEASUM-BDERval*G2(iX,iT,iY,iU)
          !! EPSA derivative
          DEPSA(iT,:) = DEPSA(iT,:)+BDERval*G2(iX,:,iY,iU)
          DEPSA(iU,:) = DEPSA(iU,:)+BDERval*G2(iX,iT,iY,:)
          DEPSA(iX,:) = DEPSA(iX,:)+BDERval*G2(:,iT,iY,iU)
          DEPSA(iY,:) = DEPSA(iY,:)+BDERval*G2(iX,iT,:,iU)

          BDERval = BDERval*Two
          SDERval = SDERval*Two

          !! term 2 (dxt)
          if (iX == iT) then
            ATYU = EASUM-ET-EY-EU
            !! G1 and F1 derivative
            DF1(iY,iU) = DF1(iY,iU)-BDERval
            DG1(iY,iU) = DG1(iY,iU)+ATYU*BDERval-SDERval
            !! EASUM derivative
            DEASUM = DEASUM+BDERval*G1(iY,iU)
            !! EPSA derivative
            DEPSA(iY,:) = DEPSA(iY,:)-BDERval*G1(:,iU)
            DEPSA(iU,:) = DEPSA(iU,:)-BDERval*G1(iY,:)
          end if
          !! Additional EPSA derivative
          DEPSA(iX,iT) = DEPSA(iX,iT)-BDERval*G1(iY,iU)
          !! dxt*dyu term
          if (iY == iU) DEPSA(iX,iT) = DEPSA(iX,iT)+Two*BDERval
          if (iX == iT) DEPSA(iY,iU) = DEPSA(iY,iU)+Two*BDERval

          !! term 3 (dyu)
          if (iY == iU) then
            ATYX = EASUM-ET-EY-EX
            !! G1 and F1 derivative
            DF1(iX,iT) = DF1(iX,iT)-BDERval
            DG1(iX,iT) = DG1(iX,iT)+ATYX*BDERval-SDERval
            !! EASUM derivative
            DEASUM = DEASUM+BDERval*G1(iX,iT)
            !! EPSA derivative
            DEPSA(iX,:) = DEPSA(iX,:)-BDERval*G1(:,iT)
            DEPSA(iT,:) = DEPSA(iT,:)-BDERval*G1(iX,:)
          end if
          !! Additional EPSA derivative
          DEPSA(iY,iU) = DEPSA(iY,iU)-BDERval*G1(iX,iT)

          BDERval = BDERval*Half
          SDERval = SDERval*Half

          !! term 4 (dyt)
          if (iY == iT) then
            ATUX = EASUM-ET-EU-EX
            !! G1 and F1 derivative
            DF1(iX,iU) = DF1(iX,iU)+BDERval
            DG1(iX,iU) = DG1(iX,iU)-ATUX*BDERval+SDERval
            !! EASUM derivative
            DEASUM = DEASUM-BDERval*G1(iX,iU)
            !! EPSA derivative
            DEPSA(iX,:) = DEPSA(iX,:)+BDERval*G1(:,iU)
            DEPSA(iU,:) = DEPSA(iU,:)+BDERval*G1(iX,:)
          end if
          !! Additional EPSA derivative
          DEPSA(iY,iT) = DEPSA(iY,iT)+BDERval*G1(iX,iU)
          !! dxu*dyt term
          if (iY == iT) DEPSA(iX,iU) = DEPSA(iX,iU)-Two*BDERval
          if (iX == iU) DEPSA(iY,iT) = DEPSA(iY,iT)-Two*BDERval

          !! term 5 (dxu)
          if (iX == iU) then
            ATUY = EASUM-ET-EU-EY
            !! G1 and F1 derivative
            DF1(iY,iT) = DF1(iY,iT)+BDERval
            DG1(iY,iT) = DG1(iY,iT)-ATUY*BDERval+SDERval
            !! EASUM derivative
            DEASUM = DEASUM-BDERval*G1(iY,iT)
            !! EPSA derivative
            DEPSA(iY,:) = DEPSA(iY,:)+BDERval*G1(:,iT)
            DEPSA(iT,:) = DEPSA(iT,:)+BDERval*G1(iY,:)
          end if
          !! Additional EPSA derivative
          DEPSA(iX,iU) = DEPSA(iX,iU)+BDERval*G1(iY,iT)
        end do
      end do
    end do
  end do
  if (ipea_shift /= Zero) call mma_deallocate(SMat)

  call mma_deallocate(WrkBbf)
  call mma_deallocate(WrkSbf)

  return

end subroutine CLagDXB

subroutine CLagDXC(NAS,BDER,SDER)

  integer(kind=iwp), intent(in) :: NAS
  real(kind=wp), intent(in) :: BDER(NAS,NAS)
  real(kind=wp), intent(inout) :: SDER(NAS,NAS)
  integer(kind=byte), allocatable :: idxG3(:,:)

  NS = nTri_Elem(NAS)
  call mma_allocate(SMat,NS,Label='SMat')
  idS = idSMAT(iSym,4)
  call DDAFILE(LUSBT,2,SMat,NS,idS)

  idum = 0
  call CLagDXC_DP(iSym,nAS,nAshT,BDER,SDER,DG1,DG2,DF1,DF2,DEPSA,DEASUM,1,nAS,1,nAS,0,G1,G2,SMat,SMat,idum)

  !! G3 and F3 relevant
  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID = 0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)
  !idS = idSMAT(iSym,4)
  !call DDAFILE(LUSBT,2,SMat,NS,idS)
  call MKSC_G3(iSym,SMat,NS,nG3,G3,idxG3)
  call CLagDXC_FG3(iSym,nAS,nAshT,NG3,NS,BDER,SDER,DF1,DF2,DF3,DG1,DG2,DG3,DEPSA,G2,SMat,idxG3)
  call mma_deallocate(idxG3)

  call mma_deallocate(SMat)

  return

end subroutine CLagDXC

subroutine CLagDXD(NAS,BDER,SDER)

  use SUPERINDEX, only: MTU

  integer(kind=iwp), intent(in) :: NAS
  real(kind=wp), intent(in) :: BDER(NAS,NAS)
  real(kind=wp), intent(inout) :: SDER(NAS,NAS)
  integer(kind=iwp) :: iTU, iXY
  real(kind=wp) :: BDER1, BDER2, ETX, SDER1, SDER2

  if (ipea_shift /= Zero) then
    NS = nTri_Elem(NAS)
    call mma_allocate(SMat,NS,Label='SMat')
    idS = idSMAT(iSym,iCase)
    call DDAFILE(LUSBT,2,SMat,NS,idS)
  end if

  do iTU=1,nAS/2
    iTU2 = iTU+nAS/2
    iTUabs = iTU+nTUes(iSym)
    iTabs = mTU(1,iTUabs)
    iUabs = mTU(2,iTUabs)
    ET = EPSA(iTabs)
    do iXY=1,nAS/2
      iXY2 = iXY+nAS/2
      iXYabs = iXY+nTUes(iSym)
      iXabs = mTU(1,iXYabs)
      iYabs = mTU(2,iXYabs)
      EX = EPSA(iXabs)
      ETX = ET+EX

      BDER1 = BDER(iTU,iXY)-BDER(iTU,iXY2)*Half-BDER(iTU2,iXY)*Half
      BDER2 = BDER(iTU2,iXY2)

      !! Derivative of B11
      DF2(iUabs,iTabs,iXabs,iYabs) = DF2(iUabs,iTabs,iXabs,iYabs)+Two*BDER1
      DG2(iUabs,iTabs,iXabs,iYabs) = DG2(iUabs,iTabs,iXabs,iYabs)+Two*(ETX-EASUM)*BDER1
      DEASUM = DEASUM-Two*G2(iUabs,iTabs,iXabs,iYabs)*BDER1
      if (iXabs == iTabs) then
        DF1(iUabs,iYabs) = DF1(iUabs,iYabs)+Two*BDER1
        DG1(iUabs,iYabs) = DG1(iUabs,iYabs)+Two*(ET-EASUM)*BDER1
        DEASUM = DEASUM-Two*G1(iUabs,iYabs)*BDER1
      end if
      DEPSA(iTabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM)) = DEPSA(iTabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM))+ &
                                                        Two*BDER1*G2(iUabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM),iXabs,iYabs)
      DEPSA(iXabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM)) = DEPSA(iXabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM))+ &
                                                        Two*BDER1*G2(iUabs,iTabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM),iYabs)
      DEPSA(iTabs,iXabs) = DEPSA(iTabs,iXabs)+Two*G1(iUabs,iYabs)*BDER1

      !! Derivative of B22
      DF2(iXabs,iTabs,iUabs,iYabs) = DF2(iXabs,iTabs,iUabs,iYabs)-BDER2
      DG2(iXabs,iTabs,iUabs,iYabs) = DG2(iXabs,iTabs,iUabs,iYabs)-(ETX-EASUM)*BDER2
      DEASUM = DEASUM+G2(iXabs,iTabs,iUabs,iYabs)*BDER2
      if (iXabs == iTabs) then
        DF1(iUabs,iYabs) = DF1(iUabs,iYabs)+Two*BDER2
        DG1(iUabs,iYabs) = DG1(iUabs,iYabs)+Two*(EX-EASUM)*BDER2
        DEASUM = DEASUM-Two*G1(iUabs,iYabs)*BDER2
      end if
      DEPSA(iTabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM)) = DEPSA(iTabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM))- &
                                                        BDER2*G2(iXabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM),iUabs,iYabs)
      DEPSA(iXabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM)) = DEPSA(iXabs,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM))- &
                                                        BDER2*G2(NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM),iTabs,iUabs,iYabs)
      DEPSA(iXabs,iTabs) = DEPSA(iXabs,iTabs)+Two*G1(iUabs,iYabs)*BDER2

      if ((iTU == iXY) .and. (ipea_shift /= Zero)) then
        !! ipea_shift*Half*(Two-DREF(IDU)+DREF(IDT))*SD(ITU)
        bsBDER = ipea_shift*Half*BDER(iTU,iXY)
        NSEQ = nTri_Elem(iTU)
        DG1(iTabs,iTabs) = DG1(iTabs,iTabs)+bsBDER*SMat(NSEQ)
        DG1(iUabs,iUabs) = DG1(iUabs,iUabs)-bsBDER*SMat(NSEQ)
        SDER(iTU,iXY) = SDER(iTU,iXY)+bsBDER*(Two+G1(iTabs,iTabs)-G1(iUabs,iUabs))
        !! ipea_shift*Half*(Two-DREF(IDU)+DREF(IDT))*SD(ITU+NAS)
        bsBDER = ipea_shift*Half*BDER(iTU2,iXY2)
        NSEQ = nTri_Elem(iTU2)
        DG1(iTabs,iTabs) = DG1(iTabs,iTabs)+bsBDER*SMat(NSEQ)
        DG1(iUabs,iUabs) = DG1(iUabs,iUabs)-bsBDER*SMat(NSEQ)
        SDER(iTU2,iXY2) = SDER(iTU2,iXY2)+bsBDER*(Two+G1(iTabs,iTabs)-G1(iUabs,iUabs))
      end if

      SDER1 = SDER(iTU,iXY)-SDER(iTU,iXY2)*Half-SDER(iTU2,iXY)*Half
      SDER2 = SDER(iTU2,iXY2)

      !! Derivative of S11
      DG2(iUabs,iTabs,iXabs,iYabs) = DG2(iUabs,iTabs,iXabs,iYabs)+Two*SDER1
      if (iXabs == iTabs) DG1(iUabs,iYabs) = DG1(iUabs,iYabs)+Two*SDER1
      !! Derivative of S22
      DG2(iXabs,iTabs,iUabs,iYabs) = DG2(iXabs,iTabs,iUabs,iYabs)-SDER2
      if (iXabs == iTabs) DG1(iUabs,iYabs) = DG1(iUabs,iYabs)+Two*SDER2
    end do
  end do
  if (ipea_shift /= Zero) call mma_deallocate(SMat)

  return

end subroutine CLagDXD

subroutine CLagDXE(NAS,BDER,SDER)

  integer(kind=iwp), intent(in) :: NAS
  real(kind=wp), intent(in) :: BDER(NAS,NAS)
  real(kind=wp), intent(inout) :: SDER(NAS,NAS)
  integer(kind=iwp) :: IT, IU
  real(kind=wp) :: VAL

  if (ipea_shift /= Zero) then
    NS = nTri_Elem(NAS)
    call mma_allocate(SMat,NS,Label='SMat')
    idS = idSMAT(iSym,6)
    call DDAFILE(LUSBT,2,SMat,NS,idS)
    !! ipea_shift*Half*DREF(IDT)*SD(IT)
    do IT=1,NAS
      ITABS = IT+NAES(ISYM)
      VAL = ipea_shift*Half*BDER(IT,IT)
      SDER(IT,IT) = SDER(IT,IT)+G1(ITABS,ITABS)*VAL
      NSEQ = nTri_Elem(IT)
      DG1(ITABS,ITABS) = DG1(ITABS,ITABS)+SMat(NSEQ)*VAL
    end do
    call mma_deallocate(SMat)
  end if

  do IT=1,NAS
    ITABS = IT+NAES(ISYM)
    ET = EPSA(ITABS)
    do IU=1,NAS
      IUABS = IU+NAES(ISYM)
      EU = EPSA(IUABS)
      !! Derivative of the B matrix
      !! B_{tu} = -F1_{tu} + (Esum-e_t-e_u)*G1(tu)
      DG1(ITABS,IUABS) = DG1(ITABS,IUABS)+(EASUM-ET-EU)*BDER(IT,IU)
      DEASUM = DEASUM+G1(ITABS,IUABS)*BDER(IT,IU)
      DF1(ITABS,IUABS) = DF1(ITABS,IUABS)-BDER(IT,IU)
      DEPSA(ITABS,IUABS) = DEPSA(ITABS,IUABS)-sum(G1(ITABS,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM))*BDER(1:NASH(ISYM),IU)+ &
                                                  G1(IUABS,NAES(ISYM)+1:NAES(ISYM)+NASH(ISYM))*BDER(1:NASH(ISYM),IT))
      DEPSA(ITABS,IUABS) = DEPSA(ITABS,IUABS)+Two*BDER(IT,IU)
      !! Derivative of the S matrix
      DG1(ITABS,IUABS) = DG1(ITABS,IUABS)-SDER(IT,IU)
    end do
  end do

  return

end subroutine CLagDXE

subroutine CLagDXF(NAS,BDER,SDER)

  use SUPERINDEX, only: MTGTU, MTGEU

  integer(kind=iwp), intent(in) :: NAS
  real(kind=wp), intent(in) :: BDER(NAS,NAS)
  real(kind=wp), intent(inout) :: SDER(NAS,NAS)
  integer(kind=iwp) :: iTU, iXY

  if (ipea_shift /= Zero) then
    NS = nTri_Elem(NAS)
    call mma_allocate(SMat,NS,Label='SMat')
    idS = idSMAT(iSym,iCase)
    call DDAFILE(LUSBT,2,SMat,NS,idS)
  end if
  ScalB1 = Zero
  ScalB2 = Zero
  ScalS1 = Zero
  ScalS2 = Zero
  iXabs = 0
  iYabs = 0
  iTabs = 0
  iUabs = 0
  iTgeUabs = 0
  iTgtUabs = 0
  iXgeYabs = 0
  iXgtYabs = 0
  do iTU=1,nAS
    if (iCase == 8) then
      iTgeUabs = iTU+nTgeUes(iSym)
      iTabs = mTgeU(1,iTgeUabs)
      iUabs = mTgeU(2,iTgeUabs)
    else if (iCase == 9) then
      iTgtUabs = iTU+nTgtUes(iSym)
      iTabs = mTgtU(1,iTgtUabs)
      iUabs = mTgtU(2,iTgtUabs)
    end if
    do iXY=1,nAS !! iTU
      if (iCase == 8) then
        iXgeYabs = iXY+nTgeUes(iSym)
        iXabs = mTgeU(1,iXgeYabs)
        iYabs = mTgeU(2,iXgeYabs)
      else if (iCase == 9) then
        iXgtYabs = iXY+nTgtUes(iSym)
        iXabs = mTgtU(1,iXgtYabs)
        iYabs = mTgtU(2,iXgtYabs)
      end if

      BDERval = BDER(ITU,IXY)
      if ((iTU == iXY) .and. (ipea_shift /= Zero)) then
        !idT = nTri_Elem(iTabs)
        !idU = nTri_Elem(iUabs)
        NSEQ = nTri_Elem(iTU)
        bsBDER = ipea_shift*Half*BDERval
        !! ipea_shift*Half*(Four-DREF(IDT)-DREF(IDU))*SDP(ITGEU)
        DG1(iTabs,iTabs) = DG1(iTabs,iTabs)-SMat(NSEQ)*bsBDER
        DG1(iUabs,iUabs) = DG1(iUabs,iUabs)-SMat(NSEQ)*bsBDER
        SDER(ITU,IXY) = SDER(ITU,IXY)+(Four-G1(iTabs,iTabs)-G1(iUabs,iUabs))*bsBDER
      end if
      SDERval = SDER(ITU,IXY)
      if (iTabs == iUabs) then
        BDERval = Two*BDERval
        SDERval = Two*SDERval
      end if

      if (iCase == 8) then
        ScalB1 = BDERval
        ScalB2 = BDERval
        ScalS1 = SDERval
        ScalS2 = SDERval
      else if (iCase == 9) then
        ScalB1 = BDERval
        ScalB2 = -BDERval
        ScalS1 = SDERval
        ScalS2 = -SDERval
      end if

      !! Derivative of the B matrix
      !! B(tuxy) -> PREF(tx,uy)
      DEASUM = DEASUM-ScalB1*G2(iTabs,iXabs,iUabs,iYabs)-ScalB2*G2(iTabs,iYabs,iUabs,iXabs)
      if (iTabs /= iUabs) DEASUM = DEASUM-ScalB2*G2(iUabs,iXabs,iTabs,iYabs)-ScalB1*G2(iUabs,iYabs,iTabs,iXabs)

      !iTX = iTabs+nAshT*(iXabs-1)
      !iUY = iUabs+nAshT*(iYabs-1)
      !iTY = iTabs+nAshT*(iYabs-1)
      !iUX = iUabs+nAshT*(iXabs-1)

      DF2(iTabs,iXabs,iUabs,iYabs) = DF2(iTabs,iXabs,iUabs,iYabs)+ScalB1
      DF2(iTabs,iYabs,iUabs,iXabs) = DF2(iTabs,iYabs,iUabs,iXabs)+ScalB2
      if (iTabs /= iUabs) then
        DF2(iUabs,iXabs,iTabs,iYabs) = DF2(iUabs,iXabs,iTabs,iYabs)+ScalB2
        DF2(iUabs,iYabs,iTabs,iXabs) = DF2(iUabs,iYabs,iTabs,iXabs)+ScalB1
      end if
      DG2(iTabs,iXabs,iUabs,iYabs) = DG2(iTabs,iXabs,iUabs,iYabs)+ScalS1-EASUM*ScalB1
      DG2(iTabs,iYabs,iUabs,iXabs) = DG2(iTabs,iYabs,iUabs,iXabs)+ScalS2-EASUM*ScalB2
      if (iTabs /= iUabs) then
        DG2(iUabs,iXabs,iTabs,iYabs) = DG2(iUabs,iXabs,iTabs,iYabs)+ScalS2-EASUM*ScalB2
        DG2(iUabs,iYabs,iTabs,iXabs) = DG2(iUabs,iYabs,iTabs,iXabs)+ScalS1-EASUM*ScalB1
      end if
    end do
  end do
  if (ipea_shift /= Zero) call mma_deallocate(SMat)

  return

end subroutine CLagDXF

subroutine CLagDXG(NAS,BDER,SDER)

  integer(kind=iwp), intent(in) :: NAS
  real(kind=wp), intent(in) :: BDER(NAS,NAS)
  real(kind=wp), intent(inout) :: SDER(NAS,NAS)
  integer(kind=iwp) :: IT, IU
  real(kind=wp) :: VAL

  if (ipea_shift /= Zero) then
    NS = nTri_Elem(NAS)
    call mma_allocate(SMat,NS,Label='SMat')
    idS = idSMAT(iSym,10)
    call DDAFILE(LUSBT,2,SMat,NS,idS)
    !! ipea_shift*Half*(Two-DREF(IDT))*SD(IT)
    do IT=1,NAS
      ITABS = IT+NAES(ISYM)
      VAL = ipea_shift*Half*BDER(IT,IT)
      SDER(IT,IT) = SDER(IT,IT)+(Two-G1(ITABS,ITABS))*VAL
      NSEQ = nTri_Elem(IT)
      DG1(ITABS,ITABS) = DG1(ITABS,ITABS)-SMat(NSEQ)*VAL
    end do
    call mma_deallocate(SMat)
  end if

  do IT=1,NAS
    ITABS = IT+NAES(ISYM)
    do IU=1,NAS
      IUABS = IU+NAES(ISYM)
      !! Derivative of the B matrix
      DG1(ITABS,IUABS) = DG1(ITABS,IUABS)-EASUM*BDER(IT,IU)
      DEASUM = DEASUM-G1(ITABS,IUABS)*BDER(IT,IU)
      DF1(ITABS,IUABS) = DF1(ITABS,IUABS)+BDER(IT,IU)
      !! Derivative of the S matrix
      DG1(ITABS,IUABS) = DG1(ITABS,IUABS)+SDER(IT,IU)
    end do
  end do

  return

end subroutine CLagDXG

#if defined(_MOLCAS_MPP_) && defined(_GA_)
subroutine CLagDX_MPP()

  use caspt2_global, only: iVecL
  use caspt2_module, only: JSTATE, MAXIT
  use Constants, only: One

  integer(kind=iwp) :: i, idB, idSD, iHiV1, iICB, iLoV1, j, jHiV1, jICB, jLoV1, LDV1, lg_BDER, lg_SDER, lg_T, lg_WRK, lg_WRK2, &
                       mBDER, mSDER, mV1, myrank, NCOL, NROW
  real(kind=wp) :: EigI, EigJ, SCAL, tmp
  logical(kind=iwp) :: bStat, invar_act
  integer(kind=byte), allocatable :: idxG3(:,:)
  real(kind=wp), allocatable :: EIG(:), WRK(:,:)

  ! Construct active density in NAS basis
  ! Although non-GA version is also implemented, I noticed that
  ! scatter operations require GA, so I should just use GA_DGEMM

  SCAL = One
  if (IFMSCOUP) SCAL = VECROT(jState)
  MYRANK = GA_NODEID()

  invar_act = .true.
  if (sigma_p_epsilon /= Zero) invar_act = .false.
  !! First, distribute the transformation matrix
  call GA_CREATE_STRIPED('H',NAS,NIN,'TRANS',lg_T)
  call PSBMAT_READ('T',iCase,iSym,lg_T,NAS*NIN)
  !call PSBMAT_READ('T',iCase,iSym,lg_T,NAS)
  !! only the master node has written it
  !if ((.not. ifmscoup) .and. King()) then
  !  call mma_allocate(TRANS,NAS*NIN,Label='TRANS')
  !  IDT = IDTMAT(ISYM,ICASE)
  !  call DDAFILE(LUSBT,2,TRANS,NAS*NIN,IDT)
  !  call GA_PUT(lg_T,1,NAS,1,NIN,TRANS,NAS)
  !  call mma_deallocate(TRANS)
  !endif
  call GA_SYNC()

  !! Allocate BDER in NIN
  call GA_CREATE_STRIPED('V',NIN,NIN,'WRK',lg_WRK)
  call GA_ZERO(lg_WRK)
  call GA_SYNC()

  ! mode = 0 operations for B derivative

  !! First, construct the density with NIN basis
  !! lg_V1 = T (solution; not quasi-variational)
  call RHS_ALLO(nIN,nIS,lg_V1)
  call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
  call GA_DGEMM('N','T',NIN,NIN,NIS,SCAL,lg_V1,lg_V1,Zero,lg_WRK)

  if ((real_shift /= Zero) .or. (imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero) .or. IFMSCOUP) then
    if (sigma_p_epsilon /= Zero) then
      !! lg_V2 = lambda (shift correction)
      call RHS_ALLO(nIN,nIS,lg_V2)
      call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
      call mma_allocate(LBD,nAS,Label='LBD')
      call mma_allocate(LID,nIS,Label='LID')
      iD = iDBMat(iSym,iCase)
      call dDaFile(LUSBT,2,LBD,nAS,iD)
      call dDaFile(LUSBT,2,LID,nIS,iD)
      !! this scaling is needed, so GA_DGEMM cannot be used
      call CASPT2_ResD(3,nIN,nIS,lg_V2,lg_V1,LBD,LID)
      call mma_deallocate(LBD)
      call mma_deallocate(LID)
    else
      !! lg_V2 = lambda (shift correction)
      call RHS_ALLO(nIN,nIS,lg_V2)
      call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
    end if
    call GA_DGEMM('N','T',NIN,NIN,NIS,Half,lg_V1,lg_V2,One,lg_WRK)
    call GA_DGEMM('N','T',NIN,NIN,NIS,Half,lg_V2,lg_V1,One,lg_WRK)
    call RHS_FREE(lg_V2)
  end if

  !if (sigma_p_epsilon /= Zero) call RHS_READ_SR(lg_V2,ICASE,ISYM,iVecR)

  call GA_SYNC()

  ! mode = 1 operations for B derivative
  ! lg_V1 is still loaded

  if ((imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero)) then
    call GA_Scale(lg_WRK,-One)

    !! T*T is skipped

    !! lg_V2 = lambda (shift correction)
    call RHS_ALLO(nIN,nIS,lg_V2)
    call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
    !call GA_Distribution (lg_V1,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
    !call GA_Distribution (lg_V2,myRank,iLoV2,iHiV2,jLoV2,jHiV2)
    !NROW1 = iHiV1-iLoV1+1
    !NCOL1 = jHiV1-jLoV1+1
    !NROW2 = iHiV2-iLoV2+1
    !NCOL2 = jHiV2-jLoV2+1
    !if (nrow1*ncol1*nrow2*ncol2 > 0) then
    call mma_allocate(LBD,nAS,Label='LBD')
    call mma_allocate(LID,nIS,Label='LID')
    iD = iDBMat(iSym,iCase)
    call dDaFile(LUSBT,2,LBD,nAS,iD)
    call dDaFile(LUSBT,2,LID,nIS,iD)
    !! this scaling is needed, so GA_DGEMM cannot be used
    call CASPT2_ResD(2,nIN,nIS,lg_V1,lg_V2,LBD,LID)
    call mma_deallocate(LBD)
    call mma_deallocate(LID)
    !call GA_Access(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
    !call GA_Access(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2,mV2,LDV2)
    !call CLag_MPP_NT(Half,DBL_MB(mV1),NROW1,DBL_MB(mV2),NROW2,lg_WRK,NIN,NCOL1)
    !call CLag_MPP_NT(Half,DBL_MB(mV2),NROW2,DBL_MB(mV1),NROW1,lg_WRK,NIN,NCOL2)
    !call GA_Release(lg_V1,iLoV1,iHiV1,jLoV1,jHiV1)
    !call GA_Release(lg_V2,iLoV2,iHiV2,jLoV2,jHiV2)
    call GA_DGEMM('N','T',NIN,NIN,NIS,Half,lg_V1,lg_V2,One,lg_WRK)
    call GA_DGEMM('N','T',NIN,NIN,NIS,Half,lg_V2,lg_V1,One,lg_WRK)
    !end if
    call RHS_FREE(lg_V2)

    !! Restore the original T
    call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
    call GA_Scale(lg_WRK,-One)
    call GA_SYNC()
  end if

  call RHS_FREE(lg_V1)

  ! B derivative in NIN completed

  if (invar_act) then
    call GA_CREATE_STRIPED('H',NAS,NIN,'WRK2',lg_WRK2)
    !! Use the same stripe as PSBMAT_GEMEM?
    call GA_CREATE_STRIPED('H',NAS,NAS,'BDER',lg_BDER)

    !! NIN -> NAS transformation of B derivative
    !! Need 4 GAs; is it possible to reduce?
    !! lg_WRK is used later, so probably not
    call GA_DGEMM('N','N',NAS,NIN,NIN,One,lg_T,lg_WRK,Zero,lg_WRK2)
    call GA_DGEMM('N','T',NAS,NAS,NIN,One,lg_WRK2,lg_T,Zero,lg_BDER)
    !if (King()) then
    !  call mma_allocate(VEC1,NAS*NAS,Label='WRK1')
    !  call GA_GET(lg_bder,1,NAS,1,NAS,VEC1,NAS)
    !  write(u6,*) 'B DERIVATIVE IN NAS'
    !  call SQPRT(VEC1,NAS)
    !  call mma_deallocate(VEC1)
    !end if

    !! cannot destroy lg_WRK; it is used for overlap derivative
    bStat = GA_destroy(lg_WRK2)
  else
    !! allocate BDER temporarily with the same dimension of WRK
    call GA_CREATE_STRIPED('V',NIN,NIN,'BDER',lg_BDER)
    call GA_COPY(lg_WRK,lg_BDER)
    call GA_ZERO(lg_WRK)
  end if

  !! At present, lg_T, lg_WRK, and lg_BDER are in GA.
  !! It is possible to destroy lg_T here and can reduce the max
  !! allocated GA by one, but we anyway have to restore the whole
  !! transformation matrix in memory in the master node, so we can
  !! probably assume that there is no problem in the availability
  !! of memory.

  !! It seems that it is not good to decompose the B derivative
  !! matrix to G and F derivative contributions here, because some
  !! contributions are computed as overlap derivatives, so it is
  !! best to decompose later simultaneously. Fortunately, the max
  !! number of allocated GAs remains unchanged (4 at most).

  !! at this point, lg_T, lg_WRK, and lg_BDER are distributed
  !if (king()) then
  !  call mma_allocate(VEC1,NAS*NAS,Label='WRK1')
  !  call GA_GET(lg_wrk,1,NIN,1,NIN,VEC1,NIN)
  !  write(u6,*) 'B DERIVATIVE IN NIN'
  !  call SQPRT(VEC1,NIN)
  !  call mma_deallocate(VEC1)
  !end if

  !mode = 0 operations for S derivative

  !! Scale with the eigenvalue
  if (invar_act) then
    call GA_Distribution(lg_WRK,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
    NROW = iHiV1-iLoV1+1
    NCOL = jHiV1-jLoV1+1
    if ((NROW > 0) .and. (NCOL > 0)) then
      call mma_allocate(EIG,NIN,Label='EIG')
      idB = idBMAT(iSym,iCase)
      call DDAFILE(LUSBT,2,EIG,NIN,IDB)
      call GA_Access(lg_WRK,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)

      do j=1,NCOL
        jICB = j+jLoV1-1
        EigJ = EIG(jICB)
        ! can't use array statement because DBL_MB is out of bounds!
        !DBL_MB(mV1+NROW*(j-1):mV1+NROW*j-1) = -DBL_MB(mV1+NROW*(j-1):mV1+NROW*j-1)*(EIG(iLoV1:iLoV1+NROW-1)+EigJ)*Half
        do i=1,NROW
          iICB = i+iLoV1-1
          EigI = EIG(iICB)
          DBL_MB(mV1+i-1+NROW*(j-1)) = -DBL_MB(mV1+i-1+NROW*(j-1))*(EigI+EigJ)*Half
        end do
      end do

      call GA_Release_Update(lg_WRK,iLoV1,iHiV1,jLoV1,jHiV1)
      call mma_deallocate(EIG)
    end if
  end if
  !if (King()) then
  !  call mma_allocate(VEC1,NIN*NIN,Label='VEC1'))
  !  call GA_GET(lg_wrk,1,NIN,1,NIN,VEC1,NIN)
  !  write(u6,*) 'SCALED B DERIVATIVE IN NIN'
  !  call SQPRT(VEC1,NIN)
  !  call mma_deallocate(VEC1)
  !end if

  !! Max memory allocation of 5 GAs
  ! 1) Implicit overlap derivative of the 2<1|H|0> part
  !! lg_V1 = T (solution; not quasi-variational)
  call RHS_ALLO(nIN,nIS,lg_V1)
  call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
  !! lg_V2 = RHS2 (in IC basis)
  call RHS_ALLO(nIN,nIS,lg_V2)
  if (ifmscoup) then
    call RHS_READ_SR(lg_V2,iCase,iSym,iVecL)
  else
    call RHS_READ_SR(lg_V2,iCase,iSym,iRHS)
  end if
  call GA_DGEMM('N','T',NIN,NIN,NIS,-One,lg_V2,lg_V1,One,lg_WRK)
  call RHS_FREE(lg_V1)
  call RHS_FREE(lg_V2)

  if ((real_shift /= Zero) .or. (imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero) .or. IFMSCOUP) then
    !! WRK1 = -RHS*(T+lambda/2)
    !! lg_V1 = RHS (in IC basis)
    call RHS_ALLO(nIN,nIS,lg_V1)
    call RHS_READ_SR(lg_V1,iCase,iSym,iRHS)
    !! lg_V2 = lambda (shift correction)
    call RHS_ALLO(nIN,nIS,lg_V2)
    call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
    call GA_DGEMM('N','T',NIN,NIN,NIS,-Half,lg_V1,lg_V2,One,lg_WRK)
    call RHS_FREE(lg_V1)
    call RHS_FREE(lg_V2)
  end if

  ! S derivative in NIN completed (some NAS operations remain)

  !! For sigma_p, non-canonical condition
  if (.not. invar_act) then
    call mma_allocate(WRK,NIN,NIN,Label='WRK1')
    call GA_GET(lg_WRK,1,NIN,1,NIN,WRK,NIN)
    call GA_Distribution(lg_BDER,myRank,iLoV1,iHiV1,jLoV1,jHiV1)
    NROW = iHiV1-iLoV1+1
    NCOL = jHiV1-jLoV1+1
    if ((NROW > 0) .and. (NCOL > 0)) then
      call mma_allocate(EIG,NIN,Label='EIG')
      idB = idBMAT(iSym,iCase)
      call DDAFILE(LUSBT,2,EIG,NIN,IDB)
      call GA_Access(lg_BDER,iLoV1,iHiV1,jLoV1,jHiV1,mBDER,LDV1)
      !! construct the off-diagonal
      do i=1,NROW
        iICB = i+iLoV1-1
        EigI = EIG(iICB)
        do j=1,NCOL
          jICB = j+jLoV1-1
          EigJ = EIG(jICB)
          if (iICB == jICB) cycle
          tmp = WRK(iICB,jICB)-WRK(jICB,iICB)
          tmp = tmp/(EigI-EigJ)
          DBL_MB(mBDER+i-1+NROW*(j-1)) = tmp
        end do
      end do
      !! -(e_o + e_p)*dS/da
      call GA_Access(lg_WRK,iLoV1,iHiV1,jLoV1,jHiV1,mV1,LDV1)
      do j=1,NCOL
        jICB = j+jLoV1-1
        EigJ = EIG(jICB)
        ! can't use array statement because DBL_MB is out of bounds!
        !DBL_MB(mV1+NROW*(j-1):mV1+NROW*j-1) = DBL_MB(mV1+NROW*(j-1):mV1+NROW*j-1)- &
        !                                      DBL_MB(mBDER+NROW*(j-1):mBDER+NROW*j-1)*(EIG(iLoV1:iLoV1+NROW-1)+EigJ)*Half
        do i=1,NROW
          iICB = i+iLoV1-1
          EigI = EIG(iICB)
          DBL_MB(mV1+i-1+NROW*(j-1)) = DBL_MB(mV1+i-1+NROW*(j-1))-DBL_MB(mBDER+i-1+NROW*(j-1))*(EigI+EigJ)*Half
        end do
      end do
      call GA_Release_Update(lg_BDER,iLoV1,iHiV1,jLoV1,jHiV1)
      call GA_Release_Update(lg_WRK,iLoV1,iHiV1,jLoV1,jHiV1)
      call mma_deallocate(EIG)
    end if
    call mma_deallocate(WRK)

    !! IC -> MO (B matrix)
    call GA_CREATE_STRIPED('H',NAS,NIN,'WRK2',lg_WRK2)
    call GA_DGEMM('N','N',NAS,NIN,NIN,One,lg_T,lg_BDER,Zero,lg_WRK2)
    bStat = GA_destroy(lg_BDER)
    call GA_CREATE_STRIPED('H',NAS,NAS,'BDER',lg_BDER)
    call GA_DGEMM('N','T',NAS,NAS,NIN,One,lg_WRK2,lg_T,Zero,lg_BDER)
    bStat = GA_destroy(lg_WRK2)
  end if

  call GA_CREATE_STRIPED('H',NAS,NIN,'WRK2',lg_WRK2)
  !! Use the same stripe as PSBMAT_GEMEM?
  call GA_CREATE_STRIPED('H',NAS,NAS,'SDER',lg_SDER)

  !! NIN -> NAS transformation of S derivative
  call GA_DGEMM('N','N',NAS,NIN,NIN,One,lg_T,lg_WRK,Zero,lg_WRK2)
  call GA_DGEMM('N','T',NAS,NAS,NIN,One,lg_WRK2,lg_T,Zero,lg_SDER)
  !if (King()) then
  !  call mma_allocate(VEC1,NAS*NAS,Label='VEC1')
  !  call GA_GET(lg_sder,1,NAS,1,NAS,VEC1,NAS)
  !  write(u6,*) 'S DERIVATIVE IN NAS (1)'
  !  call SQPRT(VEC1,NAS)
  !  call mma_deallocate(VEC1)
  !end if

  bStat = GA_destroy(lg_WRK)
  bStat = GA_destroy(lg_WRK2)

  !! Add some trivial contributions due to the dependence
  !! on the linearly independent space
  if (do_lindep .and. (nAS /= nIN)) call LinDepLag_MPP(lg_BDER,lg_SDER,nAS,nIN,iSym,iCase)

  ! 2) Explicit overlap derivative of the 2<1|H|0> part
  !    Again, not for imaginary shift-specific terms
  !! E = 2<1|H|0> + <1|H0-E0|1>
  !! lg_V1 = VEC1 = T (solution; not quasi-variational)
  call RHS_ALLO(nIN,nIS,lg_V1)
  call RHS_READ_SR(lg_V1,iCase,iSym,iVecX)
  call RHS_SCAL(nIN,nIS,lg_V1,SCAL)
  if ((real_shift /= Zero) .or. (imag_shift /= Zero) .or. (sigma_p_epsilon /= Zero) .or. IFMSCOUP) then
    !! lg_V2 = VEC2 = lambda (shift correction)
    call RHS_ALLO(nIN,nIS,lg_V2)
    call RHS_READ_SR(lg_V2,iCase,iSym,iVecR)
    call RHS_DAXPY(NIN,NIS,Half,lg_V2,lg_V1)
    call RHS_FREE(lg_V2)
  end if

  call GA_CREATE_STRIPED('V',NAS,NIS,'WRK',lg_WRK)
  call GA_DGEMM('N','N',NAS,NIS,NIN,One,lg_T,lg_V1,Zero,lg_WRK)

  call RHS_FREE(lg_V1)
  bStat = GA_destroy(lg_T)
  !! lg_V1 = VEC4 = RHS (in MO basis)
  call RHS_ALLO(nAS,nIS,lg_V1)
  call RHS_READ(nAS,nIS,lg_V1,iCase,iSym,iVecW)

  call GA_DGEMM('N','T',NAS,NAS,NIS,Two,lg_WRK,lg_V1,One,lg_SDER)

  call RHS_FREE(lg_V1)
  bStat = GA_destroy(lg_WRK)
  call GA_SYNC()

  !! Add the contributions from the off-diagonal coupling
  !! (i.e., CASPT2-N). Of course, this is not for imaginary shift-
  !! specific terms.
  if (MAXIT /= 0) then
    call mma_allocate(WRK,NAS,NAS,Label='WRK')
    idSD = idSDMat(iSym,iCase)
    call DDAFILE(LuSTD,2,WRK,nAS*nAS,idSD)
    call GA_Distribution(lg_SDER,myRank,ILO,IHI,JLO,JHI)
    call GA_Access(lg_SDER,ILO,IHI,JLO,JHI,mSDER,LDV)
    NROW = IHI-ILO+1
    NCOL = JHI-JLO+1
    do J=1,NCOL
      ! can't use array statement because DBL_MB is out of bounds!
      !DBL_MB(mSDER+NROW*(J-1):mSDER+NROW*J-1) = DBL_MB(mSDER+NROW*(J-1):mSDER+NROW*J-1)+WRK(ILO:ILO+NROW-1,J+JLO-1)*Half
      do I=1,NROW
        DBL_MB(mSDER+I-1+NROW*(J-1)) = DBL_MB(mSDER+I-1+NROW*(J-1))+WRK(I+ILO-1,J+JLO-1)*Half
      end do
    end do
    call GA_Release_Update(lg_SDER,ILO,IHI,JLO,JHI)
    call mma_deallocate(WRK)
  end if
  call GA_SYNC()
  !if (King()) then
  !  call mma_allocate(VEC1,NAS*NAS,Label='VEC1')
  !  call GA_GET(lg_sder,1,NAS,1,NAS,VEC1,NAS)
  !  write(u6,*) 'S DERIVATIVE IN NAS'
  !  call SQPRT(VEC1,NAS)
  !  call mma_deallocate(VEC1)
  !end if

  !! Compute G and F derivative contributions for B and S der
  call PSBMAT_GETMEM('S',lg_S,NAS)
  call PSBMAT_READ('S',iCase,iSym,lg_S,NAS)
  call GA_Distribution(lg_BDER,myRank,ILO,IHI,JLO,JHI)
  call GA_Distribution(lg_SDER,myRank,ILO,IHI,JLO,JHI)
  call GA_Distribution(lg_S,myRank,ILO,IHI,JLO,JHI)
  call GA_Access(lg_BDER,ILO,IHI,JLO,JHI,mBDER,LDV)
  call GA_Access(lg_SDER,ILO,IHI,JLO,JHI,mSDER,LDV)
  call GA_Access(lg_S,ILO,IHI,JLO,JHI,mS,LDV)
  NROW = iHi-iLo+1+1 !! 1 is added so that the work space is used
  NCOL = jHi-jLo+1+1 !! for all procs
  call mma_allocate(WRK,NROW,NCOL,Label='WRK')
  if (iCase == 1) then
    call CLagDXA_DP(iSym,nAS,nAshT,DBL_MB(mBDER),DBL_MB(mSDER),DG1,DG2,DF1,DF2,DEPSA,DEASUM,ILO,IHI,JLO,JHI,LDV,G1,G2,DBL_MB(mS), &
                    WRK,lg_S)
  else if (iCase == 4) then
    call CLagDXC_DP(iSym,nAS,nAshT,DBL_MB(mBDER),DBL_MB(mSDER),DG1,DG2,DF1,DF2,DEPSA,DEASUM,ILO,IHI,JLO,JHI,LDV,G1,G2,DBL_MB(mS), &
                    WRK,lg_S)
  else
    write(u6,*) 'Invalid iCase in ...'
    call abend()
  end if
  call mma_deallocate(WRK)
  call GA_Release_Update(lg_BDER,ILO,IHI,JLO,JHI)
  call GA_Release_Update(lg_SDER,ILO,IHI,JLO,JHI)

  call ga_sync()

  call mma_allocate(idxG3,6,NG3,label='idxG3')
  iLUID = 0
  call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

  if (iCase == 1) then
    call CLagDXA_FG3_MPP(iSym,NASHT,NG3,lg_BDER,lg_SDER,DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G2,idxG3)
  else if (iCase == 4) then
    call CLagDXC_FG3_MPP(iSym,NASHT,NG3,lg_BDER,lg_SDER,DG1,DG2,DG3,DF1,DF2,DF3,DEPSA,G2,idxG3)

    !! DF3 is done after icase=4
    !! construct G3 matrix in lg_S
    !call MKSC_G3_MPP(ISYM,DBL_MB(mS),ILO,NAS,LDV,NG3,G3,IDXG3)
    !call GA_Release_Update(lg_S,ILO,IHI,JLO,JHI)
    !call DF3_DEPSA_MPP(DF3,DEPSA,lg_S,idxG3)
  else
    write(u6,*) 'Invalid iCase in ...'
    call abend()
  end if
  call GA_Release(lg_S,ILO,IHI,JLO,JHI)
  call mma_deallocate(idxG3)
  call PSBMAT_FREEMEM(lg_S)

  bStat = GA_destroy(lg_BDER)
  bStat = GA_destroy(lg_SDER)

# include "macros.fh"
  unused_var(bStat)

end subroutine CLagDX_MPP
#endif

end subroutine CLagD
