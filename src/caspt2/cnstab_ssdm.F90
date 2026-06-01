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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine CnstAB_SSDM(NBSQT,DPT2AO,SSDM)

use ChoVec_io, only: NVLOC_CHOBATCH
use Cholesky, only: InfVec, nDimRS, nnBstR
use ChoCASPT2, only: MaxVec_PT2, MXNVC, NCHSPC, NumCho_PT2
use caspt2_global, only: LuAPT2, LuGAMMA
use caspt2_module, only: NBAS, NBAST, NBTCHES, NSYM
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par, myRank, nProcs
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
#include "warnings.h"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
integer(kind=iwp), intent(in) :: NBSQT
real(kind=wp), intent(in) :: DPT2AO(NBSQT), SSDM(NBSQT)
integer(kind=iwp) :: IBATCH, IBATCH_TOT, id, ILOC, iost, ipV1, ipV2, ipVecL, ipWRK(8), IRC, iSkip(8), iSym, iVec, JBATCH, &
                     JBATCH_TOT, JNUM, JRED, JRED1, JRED2, JREDC, JREDL, JSTART, jSym, JV1, JV2, KNUM, KV1, KV2, lRealName, lscr, &
                     MUSED, nBasI, NBATCH, NUMV, NUMVI, NUMVJ, NVECS_RED
logical(kind=iwp) :: is_error
character(len=4096) :: RealName
real(kind=wp), allocatable :: A_PT2(:), B_SSDM(:), CHSPC(:), HTVec(:), V1(:), V2(:), WRK(:)
integer(kind=iwp), external :: isFreeUnit
real(kind=wp), external :: ddot_
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: i, IHIV1, ILOV1, iRank, JHIV1, JLOV1, JVG, lg_V1, MV1, NDIM1
logical(kind=iwp) :: bStat
integer(kind=iwp), allocatable :: map2(:)
#endif

iSym = 1 !! iSym0

do jSym=1,nSym
  iSkip(jSym) = 1
# ifdef _MOLCAS_MPP_
  !! I do not know why this is necessary
  if (is_real_par()) then
    ipWRK(jSym) = 1
  else
# endif
    ipWRK(jSym) = 1
# ifdef _MOLCAS_MPP_
  end if
# endif
end do

nBasI = nBas(iSym)

call mma_allocate(A_PT2,MaxVec_PT2**2,Label='A_PT2')

! Read A_PT2 from LUAPT2
id = 0
call ddafile(LUAPT2,2,A_PT2,MaxVec_PT2**2,id)

!! Open B_PT2
call PrgmTranslate('GAMMA',RealName,lRealName)
LuGAMMA = isFreeUnit(LuGAMMA)
call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),'DIRECT','UNFORMATTED',iost,.true.,nBas(iSym)**2*8,'OLD',is_error)

call mma_allocate(CHSPC,NCHSPC,Label='CHSPC')
call mma_allocate(HTVec,nBasT*nBasT,Label='HTVec')
call mma_allocate(WRK,nBasT**2,Label='WRK')
!! V(P) = (mu nu|P)*D_{mu nu}
call mma_allocate(V1,MaxVec_PT2,Label='V1')
call mma_allocate(V2,MaxVec_PT2,Label='V2')
!! B_SSDM(mu,nu,P) = D_{mu rho}*D_{nu sigma}*(rho sigma|P)
!! Add one more vector
call mma_allocate(B_SSDM,NCHSPC+NBSQT,Label='B_SSDM')

if (NUMCHO_PT2(iSym) == 0) return

!ipnt = ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(iSym-1))
!JRED1 = iWork(ipnt)
!JRED2 = iWork(ipnt-1+NumCho_PT2(iSym))
JRED1 = InfVec(1,2,iSym)
JRED2 = InfVec(NumCho_PT2(iSym),2,iSym)

#ifdef _MOLCAS_MPP_
if (is_real_par()) then
  call mma_allocate(MAP2,nProcs,Label='MAP2')
  MAP2(:) = 0
  MAP2(myRank+1) = sum(NumCho_PT2(1:nSym)) ! MJRED2-JRED1+1
  call GAIGOP(MAP2,NPROCS,'+')
  !ndim2 = sum(map2)

  do i=nprocs,2,-1
    map2(i) = sum(map2(1:i-1))+1
  end do
  map2(1) = 1
  ipV1 = map2(myRank+1)
  ipV2 = map2(myRank+1)
  V1(:) = Zero
  V2(:) = Zero

  bStat = GA_CREATE_IRREG(MT_DBL,nBasT**2,MaxVec_PT2,'WRK',1,1,MAP2,NPROCS,lg_V1)
  call GA_DISTRIBUTION(LG_V1,MYRANK,ILOV1,IHIV1,JLOV1,JHIV1)
  call GA_ACCESS(LG_V1,ILOV1,IHIV1,JLOV1,JHIV1,MV1,NDIM1)
else
#endif
  ipV1 = 1
  ipV2 = 1
#ifdef _MOLCAS_MPP_
end if
#endif

!! Prepare density matrix
!! subtract the state-averaged density matrix

IBATCH_TOT = NBTCHES(iSym)

! Loop over JRED
do JRED=JRED1,JRED2

  call Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
  if (NVECS_RED == 0) cycle

  ILOC = 3
  call CHO_X_SETRED(IRC,ILOC,JRED)
  ! For a reduced set, the structure is known, including
  ! the mapping between reduced index and basis set pairs.
  ! The reduced set is divided into suitable batches.
  ! First vector is JSTART. Nr of vectors in r.s. is NVECS_RED.

  ! Determine batch length for this reduced set.
  ! Make sure to use the same formula as in the creation of disk
  ! address tables, etc, above:
  NBATCH = 1+(NVECS_RED-1)/MXNVC

  ! Loop over IBATCH
  JV1 = JSTART
  do IBATCH=1,NBATCH
    !write(u6,*) 'ibatch,nbatch = ',ibatch,nbatch
    IBATCH_TOT = IBATCH_TOT+1

    JNUM = NVLOC_CHOBATCH(IBATCH_TOT)
    JV2 = JV1+JNUM-1

    JREDC = JRED
    ! Read a batch of reduced vectors
    call CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,iSym,NUMV,JREDC,MUSED)
    if (NUMV /= JNUM) then
      write(u6,*) ' Rats! CHO_VECRD was called, assuming it to'
      write(u6,*) ' read JNUM vectors. Instead it returned NUMV'
      write(u6,*) ' vectors: JNUM, NUMV=',JNUM,NUMV
      write(u6,*) ' Back to the drawing board?'
      call QUIT(_RC_INTERNAL_ERROR_)
    end if
    if (JREDC /= JRED) then
      write(u6,*) ' Rats! It was assumed that the Cholesky vectors'
      write(u6,*) ' in HALFTRNSF all belonged to a given reduced'
      write(u6,*) ' set, but they don''t!'
      write(u6,*) ' JRED, JREDC:',JRED,JREDC
      write(u6,*) ' Back to the drawing board?'
      write(u6,*) ' Let the program continue and see what happens.'
    end if

    ipVecL = 1
    do iVec=1,NUMV

      !! reduced form -> squared AO vector (mu nu|iVec)
      !if (l_NDIMRS < 1) then
      if (size(nDimRS) < 1) then
        lscr = NNBSTR(iSym,3)
      else
        JREDL = INFVEC(iVec,2,iSym)
        !lscr = iWork(ip_nDimRS+iSym-1+nSym*(JREDL-1)) !! JRED?
        lscr = nDimRS(iSym,JREDL)
      end if
      WRK(:) = Zero
      call Cho_ReOrdr(irc,CHSPC(ipVecL),lscr,1,1,1,1,iSym,JREDC,2,ipWRK,WRK,iSkip)
      ipVecL = ipVecL+lscr

      V1(ipV1+JV1+iVec-2) = DDot_(nBasI**2,DPT2AO,1,WRK,1)
      V2(ipV2+JV1+iVec-2) = DDot_(nBasI**2,SSDM,1,WRK,1)

      call DGemm_('N','N',nBasI,nBasI,nBasI,One,DPT2AO,nBasI,WRK,nBasI,Zero,HTVec,nBasI)
      call DGemm_('N','N',nBasI,nBasI,nBasI,One,HTVec,nBasI,SSDM,nBasI,Zero,B_SSDM(1+nBasT**2*(iVec-1)),nBasI)
    end do
    NUMVI = NUMV

#   ifdef _MOLCAS_MPP_
    if (is_real_par()) then
      ! can't use array statement because DBL_MB is out of bounds!
      !DBL_MB(mV1+NDIM1*(JV1-1):mV1+NDIM1*(JV1-1)+nBasT**2*NUMV-1) = B_SSDM(1:nBasT**2*NUMV)
      call DCopy_(nBasT**2*NUMV,B_SSDM,1,DBL_MB(mV1+NDIM1*(JV1-1)),1)
    else
#   endif
      KV1 = JSTART
      JBATCH_TOT = NBTCHES(iSym)
      do JBATCH=1,NBATCH
        JBATCH_TOT = JBATCH_TOT+1

        KNUM = NVLOC_CHOBATCH(JBATCH_TOT)
        KV2 = KV1+KNUM-1

        JREDC = JRED
        call CHO_VECRD(CHSPC,NCHSPC,KV1,KV2,iSym,NUMV,JREDC,MUSED)
        call R2FIP(CHSPC,size(CHSPC),WRK,ipWRK,NUMV,nBasT,iSym,iSkip,irc,JREDC)

        !! Exchange part of A_PT2
        NUMVJ = NUMV
        call DGEMM_('T','N',NUMVI,NUMVJ,nBasT**2,-One,B_SSDM,nBasT**2,CHSPC,nBasT**2,One,A_PT2(JV1+MaxVec_PT2*(KV1-1)),MaxVec_PT2)
        KV1 = KV1+KNUM
      end do
#   ifdef _MOLCAS_MPP_
    end if
#   endif

    !! Read, add, and save the B_PT2 contribution
    do iVec=1,NUMVI
      read(LuGAMMA,rec=JV1+iVec-1) WRK(1:nBasT**2)
      !! The contributions are doubled,
      !! because halved in PGet1_RI3?
      !! Coulomb
      WRK(1:nBasT**2) = WRK(1:nBasT**2)+V2(ipV2+iVec-1)*DPT2AO(1:nBasT**2)+V1(ipV1+iVec-1)*SSDM(1:nBasT**2)
      !! Exchange
      WRK(1:nBasT**2) = WRK(1:nBasT**2)-B_SSDM(1+nBasT**2*(iVec-1):nBasT**2*iVec)
      write(LuGAMMA,rec=JV1+iVec-1) WRK(1:nBasT**2)
    end do
    JV1 = JV1+JNUM
  end do
end do

#ifdef _MOLCAS_MPP_
if (is_real_par()) then
  !! Parallel for the exchange part of A_PT2
  A_PT2(:) = A_PT2/real(NPROCS,kind=wp)
  IBATCH_TOT = NBTCHES(iSym)
  do JRED=JRED1,JRED2
    call Cho_X_nVecRS(JRED,iSym,JSTART,NVECS_RED)
    if (NVECS_RED == 0) cycle
    ILOC = 3
    call CHO_X_SETRED(IRC,ILOC,JRED)
    NBATCH = 1+(NVECS_RED-1)/MXNVC
    ! Loop over IBATCH
    JV1 = JSTART
    do IBATCH=1,NBATCH
      !write(u6,*) 'ibatch,nbatch = ',ibatch,nbatch
      IBATCH_TOT = IBATCH_TOT+1

      JNUM = NVLOC_CHOBATCH(IBATCH_TOT)
      JV2 = JV1+JNUM-1

      JREDC = JRED
      ! Read a batch of reduced vectors
      call CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,iSym,NUMV,JREDC,MUSED)
      call R2FIP(CHSPC,size(CHSPC),WRK,ipWRK,NUMV,nBasT,iSym,iSkip,irc,JREDC)
      NUMVJ = NUMV

      !! Exchange part of A_PT2
      JVG = JV1+MAP2(myRank+1)-1
      do iRank=0,NPROCS-1
        call GA_DISTRIBUTION(LG_V1,iRank,ILOV1,IHIV1,JLOV1,JHIV1)
        call GA_GET(LG_V1,ILOV1,IHIV1,JLOV1,JHIV1,B_SSDM,NDIM1)
        NUMVI = JHIV1-JLOV1+1
        call DGEMM_('T','N',NUMVI,NUMVJ,nBasT**2,-One,B_SSDM,nBasT**2,CHSPC,nBasT**2,One,A_PT2(JLOV1+MaxVec_PT2*(JVG-1)),MaxVec_PT2)
      end do
      JV1 = JV1+JNUM
    end do
  end do
  call GADGOP(A_PT2,MaxVec_PT2**2,'+')
  bStat = GA_Destroy(lg_V1)
# include "macros.fh"
  unused_var(bStat)

  call GADGOP(V1,MaxVec_PT2,'+')
  call GADGOP(V2,MaxVec_PT2,'+')
  call mma_deallocate(MAP2)
end if
#endif

!! Coulomb for A_PT2
!! Consider using DGER?
call DGEMM_('N','T',MaxVec_PT2,MaxVec_PT2,1,Two,V1,MaxVec_PT2,V2,MaxVec_PT2,One,A_PT2,MaxVec_PT2)

! write to A_PT2 in LUAPT2
id = 0
call ddafile(LUAPT2,1,A_PT2,MaxVec_PT2**2,id)

!! close B_PT2
close(LuGAMMA)

call mma_deallocate(A_PT2)

call mma_deallocate(CHSPC)
call mma_deallocate(HTVec)
call mma_deallocate(WRK)
call mma_deallocate(V1)
call mma_deallocate(V2)
call mma_deallocate(B_SSDM)

end subroutine CnstAB_SSDM
