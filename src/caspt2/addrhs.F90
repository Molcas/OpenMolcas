!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module ADDRHS

#ifdef _MOLCAS_MPP_
use GA_Wrapper, only: DBL_MB, GA_NodeId
#endif
use caspt2_global, only: iParRHS
use Symmetry_Info, only: Mul
use Constants, only: Zero, One, Two, Three, Half, Quart, OneHalf
use Definitions, only: wp, iwp, u6
use SC_NEVPT2, only: Do_SC

implicit none
private

public :: ADDRHSA, ADDRHSB, ADDRHSC, ADDRHSD1, ADDRHSD2, ADDRHSE, ADDRHSF, ADDRHSG, ADDRHSH

contains

subroutine ADDRHSA(IVEC,JSYM,ISYJ,ISYX,NT,NJ,NV,NX,TJVX,nBuff,Buff,idxBuf,Cho_Bra,Cho_Ket,NCHO)

  use SUPERINDEX, only: KTUV
  use caspt2_module, only: NINDEP, NTUV, NISH, NAES, NTUVES

  integer(kind=iwp), intent(in) :: IVEC, JSYM, ISYJ, ISYX, NT, NJ, NV, NX, nBuff, NCHO
  real(kind=wp), intent(out) :: TJVX(NT,NJ,NV,NX), Buff(nBuff)
  integer(kind=iwp), intent(out) :: idxBuf(nBuff)
  real(kind=wp), intent(in) :: Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NX,NCHO)
  integer(kind=iwp) :: IBUF, ICASE, IJ, ISYM, ISYT, ISYV, IT, ITABS, IV, IVABS, IW, IW1, IW2, IX, IXABS, LDA, lg_A, NAS, NIS, NWA
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHIV, ILOV, JHIV, JLOV, LDV, MV, myRank
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *

  ISYT = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYX)
  ISYM = ISYJ
  NAS = NTUV(ISYM)
  if (NINDEP(ISYM,1) == 0 .and. .not.Do_SC) return
  NIS = NISH(ISYM)
  NWA = NAS*NIS
  if (NWA == 0) return
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Incore or partitioned option.

  !Incore = nBuff >= NWA+NT*NJ*NV*NX
  !if (.not. Incore) then
  !  write(u6,*) 'Sort out of memory in ADDRHSA'
  !  call Abend()
  !end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call DGEMM_('N','T',NT*NJ,NV*NX,NCHO,One,Cho_Bra,NT*NJ,Cho_Ket,NV*NX,Zero,TJVX,NT*NJ)

  ! Compute W(tvx,j)=(tj,vx) + FIMO(t,j)*delta(v,x)/NACTEL
  ICASE = 1
  !LWA = 1+NT*NJ*NV*NX
  LDA = NAS
  ! Read W:
  call RHS_ALLO(NAS,NIS,lg_A)
  call RHS_READ(NAS,NIS,lg_A,iCASE,iSYM,iVEC)

  if (iParRHS == 1) then
    IBUF = 0
    do IT=1,NT
      ITABS = IT+NAES(ISYT)
      do IJ=1,NJ
        do IV=1,NV
          IVABS = IV+NAES(ISYV)
          do IX=1,NX
            IXABS = IX+NAES(ISYX)
            IW1 = KTUV(ITABS,IVABS,IXABS)-NTUVES(ISYM)
            IW2 = IJ
            IW = IW1+NAS*(IW2-1)
            !SVC Buff(LWA-1+IW) = Buff(LWA-1+IW)+TJVX(IT,IJ,IV,IX)
            IBUF = IBUF+1
            idxBuf(IBUF) = IW
            Buff(IBUF) = TJVX(IT,IJ,IV,IX)
            if (IBUF == NBUFF) then
              call RHS_SCATTER(LDA,lg_A,Buff,idxBuf,IBUF)
              IBUF = 0
            end if
          end do
        end do
      end do
    end do
    if (IBUF /= 0) call RHS_SCATTER(LDA,lg_A,Buff,idxBuf,IBUF)
# ifdef _MOLCAS_MPP_
  else if (iParRHS == 2) then
    ! Note that iParRHS = 2 only when is_real_par() is true and
    ! PRHS = 2 is specified in the input file (see procinp_caspt2)
    call GADSUM_ADDRHS(TJVX,NT*NJ*NV*NX)
    myRank = GA_NodeID()
    call GA_Distribution(lg_A,myRank,ILOV,IHIV,JLOV,JHIV)
    if (JLOV > 0) then
      call GA_Access(lg_A,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      do IT=1,NT
        ITABS = IT+NAES(ISYT)
        do IJ=JLOV,JHIV
          IW2 = IJ
          do IV=1,NV
            IVABS = IV+NAES(ISYV)
            do IX=1,NX
              IXABS = IX+NAES(ISYX)
              IW1 = KTUV(ITABS,IVABS,IXABS)-NTUVES(ISYM)
              if ((IW1 >= ILOV) .and. (IW1 <= IHIV)) &
                DBL_MB(MV+IW1-ILOV+LDA*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDA*(IW2-JLOV))+TJVX(IT,IJ,IV,IX)
            end do
          end do
        end do
      end do
      call GA_Release_Update(lg_A,ILOV,IHIV,JLOV,JHIV,MV,LDV)
    end if
# endif
  end if

  ! Put W on disk:
  call RHS_SAVE(NAS,NIS,lg_A,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_A)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine ADDRHSA

subroutine ADDRHSB(IVEC,JSYM,ISYJ,ISYL,NT,NJ,NV,NL,TJVL,nBuff,Buff,idxBuf,Cho_Bra,Cho_Ket,NCHO)
  ! Case B:
  ! t>v j>l WP(tv,jl)=add ((tj,vl))*(1/2)
  ! t>v j=l WP(tv,jl)=add ((tj,vl))*(1/2)*(SQRT(2))
  ! t>v j<l WP(tv,lj)=add ((tj,vl))*(1/2)

  ! t=v j>l WP(tv,jl)=add ((tj,vl))*(1/4)
  ! t=v j=l WP(tv,jl)=add ((tj,vl))*(1/4)*(SQRT(2))
  ! t=v j<l WP(tv,lj)=add ((tj,vl))*(1/4)

  use SUPERINDEX, only: KIGEJ, KIGTJ, KTGEU, KTGTU
  use caspt2_module, only: NAES, NIES, NIGEJ, NIGEJES, NIGTJ, NIGTJES, NINDEP, NTGEU, NTGEUES, NTGTU, NTGTUES

  integer(kind=iwp), intent(in) :: IVEC, JSYM, ISYJ, ISYL, NT, NJ, NV, NL, nBuff, NCHO
  real(kind=wp), intent(out) :: TJVL(NT,NJ,NV,NL), Buff(nBuff)
  integer(kind=iwp), intent(out) :: idxBuf(nBuff)
  real(kind=wp), intent(in) :: Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NL,NCHO)
  integer(kind=iwp) :: IBUF, ICASE, IJ, IJABS, IL, ILABS, ISYM, ISYT, ISYV, IT, ITABS, IV, IVABS, IVMAX, IW, IW1, IW2, LDBM, LDBP, &
                       lg_BM, lg_BP, NASM, NASP, NISM, NISP, NWBM, NWBP
  real(kind=wp) :: SCL, SCL1, SQ2
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHIV, ILOV, JHIV, JLOV, LDV, MV, myRank
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *

  ISYT = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYL)
  if (ISYT < ISYV) return
  SQ2 = sqrt(Two)
  ISYM = Mul(ISYJ,ISYL)

  if (NINDEP(ISYM,2) > 0 .or. Do_SC) then
    ! The plus combination:
    ICASE = 2
    NASP = NTGEU(ISYM)
    NISP = NIGEJ(ISYM)
    NWBP = NASP*NISP
  else
    NWBP = 0
  end if
  if (NINDEP(ISYM,3) > 0 .or. Do_SC) then
    ! The minus combination:
    ICASE = 3
    NASM = NTGTU(ISYM)
    NISM = NIGTJ(ISYM)
    NWBM = NASM*NISM
  else
    NWBM = 0
  end if
  if (max(NWBP,NWBM) <= 0) return
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Incore or partitioned option.

  !Incore = nBuff >= max(NWBP,NWBM)+NT*NJ*NV*NL
  !if (.not. Incore) then
  !  write(u6,*) 'Sort out of memory in ADDRHSB'
  !  call Abend()
  !end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call DGEMM_('N','T',NT*NJ,NV*NL,NCHO,One,Cho_Bra,NT*NJ,Cho_Ket,NV*NL,Zero,TJVL,NT*NJ)
# ifdef _MOLCAS_MPP_
  if (iParRHS == 2) call GADSUM_ADDRHS(TJVL,NT*NJ*NV*NL)
# endif
  if (NWBP > 0) then

    if (NINDEP(ISYM,2) > 0 .or. Do_SC) then
      ! The plus combination:
      ICASE = 2
      NASP = NTGEU(ISYM)
      NISP = NIGEJ(ISYM)
      NWBP = NASP*NISP
      !LWBP = 1+NT*NJ*NV*NL
      LDBP = NASP
      ! Read WP:
      call RHS_ALLO(NASP,NISP,lg_BP)
      call RHS_READ(NASP,NISP,lg_BP,iCASE,iSYM,iVEC)

      if (iParRHS == 1) then
        IBUF = 0
        do IT=1,NT
          ITABS = IT+NAES(ISYT)
          IVMAX = NV
          if (ISYV == ISYT) IVMAX = IT
          do IV=1,IVMAX
            IVABS = IV+NAES(ISYV)
            SCL1 = Half
            IW1 = KTGEU(ITABS,IVABS)-NTGEUES(ISYM)
            if (ITABS == IVABS) SCL1 = Quart
            do IJ=1,NJ
              IJABS = IJ+NIES(ISYJ)
              do IL=1,NL
                ILABS = IL+NIES(ISYL)
                SCL = SCL1
                if (IJABS >= ILABS) then
                  IW2 = KIGEJ(IJABS,ILABS)-NIGEJES(ISYM)
                  if (IJABS == ILABS) SCL = SQ2*SCL1
                else
                  IW2 = KIGEJ(ILABS,IJABS)-NIGEJES(ISYM)
                end if
                IW = IW1+NASP*(IW2-1)
                !Buff(LWBP-1+IW) = Buff(LWBP-1+IW)+SCL*TJVL(IT,IJ,IV,IL)
                IBUF = IBUF+1
                idxBuf(IBUF) = IW
                Buff(IBUF) = SCL*TJVL(IT,IJ,IV,IL)
                if (IBUF == NBUFF) then
                  call RHS_SCATTER(LDBP,lg_BP,Buff,idxBuf,IBUF)
                  IBUF = 0
                end if
              end do
            end do
          end do
        end do
        if (IBUF /= 0) call RHS_SCATTER(LDBP,lg_BP,Buff,idxBuf,IBUF)
#     ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
        myRank = GA_NodeID()
        call GA_Distribution(lg_BP,myRank,ILOV,IHIV,JLOV,JHIV)
        if (JLOV > 0) then
          call GA_Access(lg_BP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
          do IT=1,NT
            ITABS = IT+NAES(ISYT)
            IVMAX = NV
            if (ISYV == ISYT) IVMAX = IT
            do IV=1,IVMAX
              IVABS = IV+NAES(ISYV)
              SCL1 = Half
              IW1 = KTGEU(ITABS,IVABS)-NTGEUES(ISYM)
              if ((IW1 < ILOV) .or. (IW1 > IHIV)) cycle
              if (ITABS == IVABS) SCL1 = Quart
              do IJ=1,NJ
                IJABS = IJ+NIES(ISYJ)
                do IL=1,NL
                  ILABS = IL+NIES(ISYL)
                  SCL = SCL1
                  if (IJABS >= ILABS) then
                    IW2 = KIGEJ(IJABS,ILABS)-NIGEJES(ISYM)
                    if (IJABS == ILABS) SCL = SQ2*SCL1
                  else
                    IW2 = KIGEJ(ILABS,IJABS)-NIGEJES(ISYM)
                  end if
                  if ((IW2 >= JLOV) .and. (IW2 <= JHIV)) &
                    DBL_MB(MV+IW1-ILOV+LDBP*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDBP*(IW2-JLOV))+SCL*TJVL(IT,IJ,IV,IL)
                end do
              end do
            end do
          end do
          call GA_Release_Update(lg_BP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
        end if
#     endif
      end if

      ! Put WBP on disk:
      call RHS_SAVE(NASP,NISP,lg_BP,iCASE,iSYM,iVEC)
      call RHS_FREE(lg_BP)
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (NINDEP(ISYM,3) > 0 .or. Do_SC) then
    ! The minus combination:
    ICASE = 3
    NASM = NTGTU(ISYM)
    NISM = NIGTJ(ISYM)
    NWBM = NASM*NISM
    !LWBM = 1+NT*NJ*NV*NL
    LDBM = NASM
    ! Read WM:
    call RHS_ALLO(NASM,NISM,lg_BM)
    call RHS_READ(NASM,NISM,lg_BM,iCASE,iSYM,iVEC)

    if (iParRHS == 1) then
      IBUF = 0
      do IT=1,NT
        ITABS = IT+NAES(ISYT)
        IVMAX = NV
        if (ISYV == ISYT) IVMAX = IT-1
        do IV=1,IVMAX
          IVABS = IV+NAES(ISYV)
          IW1 = KTGTU(ITABS,IVABS)-NTGTUES(ISYM)
          do IJ=1,NJ
            IJABS = IJ+NIES(ISYJ)
            do IL=1,NL
              ILABS = IL+NIES(ISYL)
              if (IJABS > ILABS) then
                IW2 = KIGTJ(IJABS,ILABS)-NIGTJES(ISYM)
                IW = IW1+NASM*(IW2-1)
                !Buff(LWBM-1+IW) = Buff(LWBM-1+IW)+Half*TJVL(IT,IJ,IV,IL)
                IBUF = IBUF+1
                idxBuf(IBUF) = IW
                Buff(IBUF) = Half*TJVL(IT,IJ,IV,IL)
              else if (IJABS < ILABS) then
                IW2 = KIGTJ(ILABS,IJABS)-NIGTJES(ISYM)
                IW = IW1+NASM*(IW2-1)
                !Buff(LWBM-1+IW) = Buff(LWBM-1+IW)-Half*TJVL(IT,IJ,IV,IL)
                IBUF = IBUF+1
                idxBuf(IBUF) = IW
                Buff(IBUF) = -Half*TJVL(IT,IJ,IV,IL)
              end if
              if (IBUF == NBUFF) then
                call RHS_SCATTER(LDBM,lg_BM,Buff,idxBuf,IBUF)
                IBUF = 0
              end if
            end do
          end do
        end do
      end do
      if (IBUF /= 0) call RHS_SCATTER(LDBM,lg_BM,Buff,idxBuf,IBUF)
#   ifdef _MOLCAS_MPP_
    else if (iParRHS == 2) then
      myRank = GA_NodeID()
      call GA_Distribution(lg_BM,myRank,ILOV,IHIV,JLOV,JHIV)
      if (JLOV > 0) then
        call GA_Access(lg_BM,ILOV,IHIV,JLOV,JHIV,MV,LDV)

        do IT=1,NT
          ITABS = IT+NAES(ISYT)
          IVMAX = NV
          if (ISYV == ISYT) IVMAX = IT-1
          do IV=1,IVMAX
            IVABS = IV+NAES(ISYV)
            IW1 = KTGTU(ITABS,IVABS)-NTGTUES(ISYM)
            if ((IW1 < ILOV) .or. (IW1 > IHIV)) cycle
            do IJ=1,NJ
              IJABS = IJ+NIES(ISYJ)
              do IL=1,NL
                ILABS = IL+NIES(ISYL)
                if (IJABS > ILABS) then
                  IW2 = KIGTJ(IJABS,ILABS)-NIGTJES(ISYM)
                  if ((IW2 >= JLOV) .and. (IW2 <= JHIV)) &
                    DBL_MB(MV+IW1-ILOV+LDBM*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDBM*(IW2-JLOV))+Half*TJVL(IT,IJ,IV,IL)
                else if (IJABS < ILABS) then
                  IW2 = KIGTJ(ILABS,IJABS)-NIGTJES(ISYM)
                  if ((IW2 >= JLOV) .and. (IW2 <= JHIV)) &
                    DBL_MB(MV+IW1-ILOV+LDBM*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDBM*(IW2-JLOV))-Half*TJVL(IT,IJ,IV,IL)
                end if
              end do
            end do
          end do
        end do
        call GA_Release_Update(lg_BM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      end if
#   endif
    end if

    ! Put WBM on disk:
    call RHS_SAVE(NASM,NISM,lg_BM,iCASE,iSYM,iVEC)
    call RHS_FREE(lg_BM)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine ADDRHSB

subroutine ADDRHSC(IVEC,JSYM,ISYU,ISYX,NA,NU,NV,NX,AUVX,nBuff,Buff,idxBuf,Cho_Bra,Cho_Ket,NCHO)
  ! Case C:
  !   Allocate W. Put in W(uvx,a)=(au,vx) +
  !             (FIMO(a,t)-sum(y)(ay,yt))*delta(u,v)/NACTEL.

  use SUPERINDEX, only: KTUV
  use caspt2_module, only: NAES, NINDEP, NSSH, NTUV, NTUVES

  integer(kind=iwp), intent(in) :: IVEC, JSYM, ISYU, ISYX, NA, NU, NV, NX, nBuff, NCHO
  real(kind=wp), intent(out) :: AUVX(NA,NU,NV,NX), Buff(nBuff)
  integer(kind=iwp), intent(out) :: idxBuf(nBuff)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NX,NCHO)
  integer(kind=iwp) :: IA, IBUF, ICASE, ISYA, ISYM, ISYV, IU, IUABS, IV, IVABS, IW, IW1, IW2, IX, IXABS, LDC, lg_C, NAS, NIS, NWC
#ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHIV, ILOV, JHIV, JLOV, LDV, MV, myRank
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *

  ISYA = Mul(JSYM,ISYU)
  ISYV = Mul(JSYM,ISYX)
  ISYM = ISYA
  if (NINDEP(ISYM,4) == 0 .and. .not.Do_SC) return
  NAS = NTUV(ISYM)
  NIS = NSSH(ISYM)
  NWC = NAS*NIS
  if (NWC == 0) return
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Incore or partitioned option.

  !Incore = nBuff >= NWC+NA*NU*NV*NX
  !if (.not. Incore) then
  !  write(u6,*) 'Sort out of memory in ADDRHSC'
  !  call Abend()
  !end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call DGEMM_('N','T',NA*NU,NV*NX,NCHO,One,Cho_Bra,NA*NU,Cho_Ket,NV*NX,Zero,AUVX,NA*NU)

  ICASE = 4
  !LWC = 1+NA*NU*NV*NX
  LDC = NAS
  ! Read W:
  call RHS_ALLO(NAS,NIS,lg_C)
  call RHS_READ(NAS,NIS,lg_C,iCASE,iSYM,iVEC)
  if (iParRHS == 1) then
    IBUF = 0
    do IA=1,NA
      do IU=1,NU
        IUABS = IU+NAES(ISYU)
        do IV=1,NV
          IVABS = IV+NAES(ISYV)
          do IX=1,NX
            IXABS = IX+NAES(ISYX)
            IW1 = KTUV(IUABS,IVABS,IXABS)-NTUVES(ISYM)
            IW2 = IA
            IW = IW1+NAS*(IW2-1)
            !Buff(LWC-1+IW) = Buff(LWC-1+IW)+AUVX(IA,IU,IV,IX)
            IBUF = IBUF+1
            idxBuf(IBUF) = IW
            Buff(IBUF) = AUVX(IA,IU,IV,IX)
            if (IBUF == NBUFF) then
              call RHS_SCATTER(LDC,lg_C,Buff,idxBuf,IBUF)
              IBUF = 0
            end if
          end do
        end do
      end do
    end do
    if (IBUF /= 0) call RHS_SCATTER(LDC,lg_C,Buff,idxBuf,IBUF)
# ifdef _MOLCAS_MPP_
  else if (iParRHS == 2) then
    call GADSUM_ADDRHS(AUVX,NA*NU*NV*NX)
    myRank = GA_NodeID()
    call GA_Distribution(lg_C,myRank,ILOV,IHIV,JLOV,JHIV)
    if (JLOV > 0) then
      call GA_Access(lg_C,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      do IA=JLOV,JHIV
        IW2 = IA
        do IU=1,NU
          IUABS = IU+NAES(ISYU)
          do IV=1,NV
            IVABS = IV+NAES(ISYV)
            do IX=1,NX
              IXABS = IX+NAES(ISYX)
              IW1 = KTUV(IUABS,IVABS,IXABS)-NTUVES(ISYM)
              if ((IW1 >= ILOV) .and. (IW1 <= IHIV)) &
                DBL_MB(MV+IW1-ILOV+LDC*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDC*(IW2-JLOV))+AUVX(IA,IU,IV,IX)
            end do
          end do
        end do
      end do
      call GA_Release_Update(lg_C,ILOV,IHIV,JLOV,JHIV,MV,LDV)
    end if
# endif
  end if
  ! Put W on disk:
  call RHS_SAVE(NAS,NIS,lg_C,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_C)

  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine ADDRHSC

subroutine ADDRHSD1(IVEC,JSYM,ISYJ,ISYX,NA,NJ,NV,NX,AJVX,nBuff,Buff,idxBuf,Cho_Bra,Cho_Ket,NCHO)
  ! Case D:
  ! Compute W1(vx,aj)=(aj,vx) + FIMO(a,j)*delta(v,x)/NACTEL
  ! Compute W2(vu,al)=(au,vl)

  use SUPERINDEX, only: KTU
  use caspt2_module, only: NAES, NINABX, NINDEP, NISH, NISUP, NSECBX, NSSH, NSYM, NTU, NTUES

  integer(kind=iwp), intent(in) :: IVEC, JSYM, ISYJ, ISYX, NA, NJ, NV, NX, nBuff, NCHO
  real(kind=wp), intent(out) :: AJVX(NV,NX,*), Buff(nBuff)
  integer(kind=iwp), intent(out) :: idxBuf(nBuff)
  real(kind=wp), intent(in) :: Cho_Bra(NA*NJ,NCHO), Cho_Ket(NV,NX,NCHO)
  integer(kind=iwp) :: IA, IAEND, IAJ, IAJSTA, IASTA, IBUF, ICASE, IJ, IJEND, IJSTA, IO, IOFFD(8,8), ISA, ISI, ISW, ISYA, ISYM, &
                       ISYV, IV, IVABS, IW, IW1, IW2, IX, IXABS, LDD, lg_D, NAS, NAS1, NASZ, NBXSZA, NBXSZJ, NIS, NJSZ, NWD
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHIV, ILOV, JHIV, JLOV, LDV, MV, myRank
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *

  do ISW=1,NSYM
    IO = 0
    do ISA=1,NSYM
      IOFFD(ISA,ISW) = IO
      ISI = Mul(ISA,ISW)
      IO = IO+NSSH(ISA)*NISH(ISI)
    end do
  end do

  ISYA = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYX)
  ISYM = JSYM
  if (NINDEP(ISYM,5) == 0 .and. .not.Do_SC) return
  NAS1 = NTU(ISYM)
  NAS = 2*NAS1
  NIS = NISUP(ISYM,5)
  NWD = NAS*NIS
  if (NWD == 0) return
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Incore or partitioned option?

  !Incore = nBuff >= NWD+NA*NJ*NV*NX
  !if (.not. Incore) then
  !  write(u6,*) 'Sort out of memory in ADDRHSD1'
  !  call Abend()
  !end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !call DGEMM_('N','T',NA*NJ,NV*NX,NCHO,One,Cho_Bra,NA*NJ,Cho_Ket,NV*NX,Zero,AJVX,NA*NJ)

  ! Compute W1(vx,aj)=(aj,vx) + FIMO(a,j)*delta(v,x)/NACTEL
  ICASE = 5
  !LWD = 1+NA*NJ*NV*NX
  LDD = NAS
  ! Read W:
  call RHS_ALLO(NAS,NIS,lg_D)
  call RHS_READ(NAS,NIS,lg_D,iCASE,iSYM,iVEC)
# ifdef _MOLCAS_MPP_
  if (iParRHS == 2) then
    myRank = GA_NodeID()
    call GA_Distribution(lg_D,myRank,ILOV,IHIV,JLOV,JHIV)
    if (JLOV > 0) call GA_Access(lg_D,ILOV,IHIV,JLOV,JHIV,MV,LDV)
  end if
# endif

  NBXSZA = NSECBX
  NBXSZJ = NINABX

  do IASTA=1,NA,NBXSZA
    IAEND = min(IASTA-1+NBXSZA,NA)
    NASZ = IAEND-IASTA+1
    do IJSTA=1,NJ,NBXSZJ
      IJEND = min(IJSTA-1+NBXSZJ,NJ)
      NJSZ = IJEND-IJSTA+1

      IAJSTA = 1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
      call DGEMM_('N','T',NV*NX,NASZ*NJSZ,NCHO,One,Cho_Ket,NV*NX,Cho_Bra(IAJSTA,1),NA*NJ,Zero,AJVX,NV*NX)

      if (iParRHS == 1) then
        IAJ = 0
        IBUF = 0
        do IJ=IJSTA,IJEND
          do IA=IASTA,IAEND
            IAJ = IAJ+1

            do IX=1,NX
              IXABS = IX+NAES(ISYX)
              do IV=1,NV
                IVABS = IV+NAES(ISYV)

                IW1 = KTU(IVABS,IXABS)-NTUES(ISYM)
                IW2 = IOFFD(ISYA,ISYM)+IJ+NJ*(IA-1)
                IW = IW1+NAS*(IW2-1)
                IBUF = IBUF+1
                !Buff(LWD-1+IW) = Buff(LWD-1+IW)+AJVX(IV,IX,IAJ)
                idxBuf(IBUF) = IW
                Buff(IBUF) = AJVX(IV,IX,IAJ)
                if (IBUF == NBUFF) then
                  call RHS_SCATTER(LDD,lg_D,Buff,idxBuf,IBUF)
                  IBUF = 0
                end if
              end do
            end do
          end do
        end do
        if (IBUF /= 0) call RHS_SCATTER(LDD,lg_D,Buff,idxBuf,IBUF)
#     ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
        call GADSUM_ADDRHS(AJVX,NV*NX*NASZ*NJSZ)
        if (JLOV > 0) then
          IAJ = 0
          do IJ=IJSTA,IJEND
            do IA=IASTA,IAEND
              IAJ = IAJ+1
              IW2 = IOFFD(ISYA,ISYM)+IJ+NJ*(IA-1)
              if ((IW2 < JLOV) .or. (IW2 > JHIV)) cycle

              do IX=1,NX
                IXABS = IX+NAES(ISYX)
                do IV=1,NV
                  IVABS = IV+NAES(ISYV)

                  IW1 = KTU(IVABS,IXABS)-NTUES(ISYM)
                  if ((IW1 >= ILOV) .and. (IW1 <= IHIV)) &
                    DBL_MB(MV+IW1-ILOV+LDD*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDD*(IW2-JLOV))+AJVX(IV,IX,IAJ)
                end do
              end do
            end do
          end do
        end if
#     endif
      end if

    end do
  end do

# ifdef _MOLCAS_MPP_
  if ((iParRHS == 2) .and. (JLOV > 0)) call GA_Release_Update(lg_D,ILOV,IHIV,JLOV,JHIV,MV,LDV)
# endif
  ! Put W on disk:
  call RHS_SAVE(NAS,NIS,lg_D,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_D)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine ADDRHSD1

subroutine ADDRHSD2(IVEC,JSYM,ISYU,ISYL,NA,NU,NV,NL,AUVL,nBuff,Buff,idxBuf,Cho_Bra,Cho_Ket,NCHO)
  ! Case D:
  ! Compute W1(vx,aj)=(aj,vx) + FIMO(a,j)*delta(v,x)/NACTEL
  ! Compute W2(vu,al)=(au,vl)

  use SUPERINDEX, only: KTU
  use caspt2_module, only: NAES, NINDEP, NISH, NISUP, NSSH, NSYM, NTU, NTUES

  integer(kind=iwp), intent(in) :: IVEC, JSYM, ISYU, ISYL, NA, NU, NV, NL, nBuff, NCHO
  real(kind=wp), intent(out) :: AUVL(NA,NU,NV,NL), Buff(nBuff)
  integer(kind=iwp), intent(out) :: idxBuf(nBuff)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NL,NCHO)
  integer(kind=iwp) :: IA, IBUF, ICASE, IL, IO, IOFFD(8,8), ISYA, ISYI, ISYM, ISYV, ISYW, IU, IUABS, IV, IVABS, IW, IW1, IW2, LDD, &
                       lg_D, NAS, NAS1, NIS, NWD
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHIV, ILOV, JHIV, JLOV, LDV, MV, myRank
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *

  do ISYW=1,NSYM
    IO = 0
    do ISYA=1,NSYM
      IOFFD(ISYA,ISYW) = IO
      ISYI = Mul(ISYA,ISYW)
      IO = IO+NSSH(ISYA)*NISH(ISYI)
    end do
  end do

  ISYA = Mul(JSYM,ISYU)
  ISYV = Mul(JSYM,ISYL)
  ISYM = Mul(ISYU,ISYV)
  if (NINDEP(ISYM,5) == 0 .and. .not.Do_SC) return
  NAS1 = NTU(ISYM)
  NAS = 2*NAS1
  NIS = NISUP(ISYM,5)
  NWD = NAS*NIS
  if (NWD == 0) return
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Incore or partitioned option?

  !Incore = nBuff >= NWD+NA*NU*NV*NL
  !if (.not. Incore) then
  !  write(u6,*) 'Sort out of memory in ADDRHSD2'
  !  call Abend()
  !end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *

  call DGEMM_('N','T',NA*NU,NV*NL,NCHO,One,Cho_Bra,NA*NU,Cho_Ket,NV*NL,Zero,AUVL,NA*NU)

  ! Compute W2(vu,al)=(au,vl)
  ICASE = 5
  !LWD = 1+NA*NU*NV*NL
  LDD = NAS
  ! Read W:
  call RHS_ALLO(NAS,NIS,lg_D)
  call RHS_READ(NAS,NIS,lg_D,iCASE,iSYM,iVEC)

  if (iParRHS == 1) then
    IBUF = 0
    do IA=1,NA
      do IU=1,NU
        IUABS = IU+NAES(ISYU)
        do IV=1,NV
          IVABS = IV+NAES(ISYV)
          do IL=1,NL
            IW1 = NAS1+KTU(IVABS,IUABS)-NTUES(ISYM)
            IW2 = IOFFD(ISYA,ISYM)+IL+NL*(IA-1)
            IW = IW1+NAS*(IW2-1)
            !Buff(LWD-1+IW) = Buff(LWD-1+IW)+AUVL(IA,IU,IV,IL)
            IBUF = IBUF+1
            idxBuf(IBUF) = IW
            Buff(IBUF) = AUVL(IA,IU,IV,IL)
            if (IBUF == NBUFF) then
              call RHS_SCATTER(LDD,lg_D,Buff,idxBuf,IBUF)
              IBUF = 0
            end if
          end do
        end do
      end do
    end do
    if (IBUF /= 0) call RHS_SCATTER(LDD,lg_D,Buff,idxBuf,IBUF)
# ifdef _MOLCAS_MPP_
  else if (iParRHS == 2) then
    call GADSUM_ADDRHS(AUVL,NA*NU*NV*NL)
    myRank = GA_NodeID()
    call GA_Distribution(lg_D,myRank,ILOV,IHIV,JLOV,JHIV)
    if (JLOV > 0) then
      call GA_Access(lg_D,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      do IA=1,NA
        do IU=1,NU
          IUABS = IU+NAES(ISYU)
          do IV=1,NV
            IVABS = IV+NAES(ISYV)
            IW1 = NAS1+KTU(IVABS,IUABS)-NTUES(ISYM)
            if ((IW1 < ILOV) .and. (IW1 > IHIV)) cycle
            do IL=1,NL
              IW2 = IOFFD(ISYA,ISYM)+IL+NL*(IA-1)
              if ((IW2 >= JLOV) .and. (IW2 <= JHIV)) &
                DBL_MB(MV+IW1-ILOV+LDD*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDD*(IW2-JLOV))+AUVL(IA,IU,IV,IL)
            end do
          end do
        end do
      end do
      call GA_Release_Update(lg_D,ILOV,IHIV,JLOV,JHIV,MV,LDV)
    end if
# endif
  end if

  ! Put W on disk:
  call RHS_SAVE(NAS,NIS,lg_D,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_D)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine ADDRHSD2

subroutine ADDRHSE(IVEC,JSYM,ISYJ,ISYL,NA,NJ,NV,NL,AJVL,nBuff,Buff,idxBuf,Cho_Bra,Cho_Ket,NCHO)
  ! Case E:

  use SUPERINDEX, only: KIGEJ, KIGTJ
  use caspt2_module, only: NASH, NIES, NIGEJ, NIGEJES, NIGTJ, NIGTJES, NINABX, NISUP, NSECBX, NSSH, NSYM

  integer(kind=iwp), intent(in) :: IVEC, JSYM, ISYJ, ISYL, NA, NJ, NV, NL, nBuff, NCHO
  real(kind=wp), intent(out) :: AJVL(NV,NL,*), Buff(nBuff)
  integer(kind=iwp), intent(out) :: idxBuf(nBuff)
  real(kind=wp), intent(in) :: Cho_Bra(NA*NJ,NCHO), Cho_Ket(NV,NL,NCHO)
  integer(kind=iwp) :: IA, IAEND, IAJ, IAJSTA, IASTA, IBUF, ICASE, IJ, IJABS, IJEND, IJSTA, IL, ILABS, IO1, IO2, IOFF1(8), &
                       IOFF2(8), ISA, ISIJ, ISYA, ISYJL, ISYM, ISYV, IV, IW, IW1, IW2, JGEL, JGTL, LDEM, LDEP, lg_EM, lg_EP, NAS, &
                       NASZ, NBXSZA, NBXSZJ, NISM, NISP, NJSZ, NW, NWM, NWP
  real(kind=wp) :: SCL, SQ32
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHIV, ILOV, JHIV, JLOV, LDV, MV, myRank
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *

  SQ32 = sqrt(OneHalf)
  ISYA = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYL)
  ISYM = ISYV
  ISYJL = Mul(ISYJ,ISYL)

  ! Set up offset table:
  IO1 = 0
  IO2 = 0
  do ISA=1,NSYM
    IOFF1(ISA) = IO1
    IOFF2(ISA) = IO2
    ISIJ = Mul(ISA,ISYM)
    IO1 = IO1+NSSH(ISA)*NIGEJ(ISIJ)
    IO2 = IO2+NSSH(ISA)*NIGTJ(ISIJ)
  end do

  NAS = NASH(ISYM)
  NISP = NISUP(ISYM,6)
  NISM = NISUP(ISYM,7)
  NWP = NAS*NISP
  NWM = NAS*NISM
  NW = NWP+NWM
  if (NW == 0) return
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Incore or partitioned option?

  !Incore = nBuff >= Max(NWP,NWM)+NA*NJ*NV*NL
  !if (.not. Incore) then
  !  write(u6,*) 'Sort out of memory in ADDRHSE'
  !  call Abend()
  !end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !call DGEMM_('N','T',NA*NJ,NV*NL,NCHO,One,Cho_Bra,NA*NJ,Cho_Ket,NV*NL,Zero,AJVL,NA*NJ)
  !
  !LWE = 1+NA*NJ*NV*NL
  !LWEP = LWE
  !LWEM = LWE
  LDEP = NAS
  LDEM = NAS
  ! The plus combination:
  if (NWP > 0) then
    ICASE = 6
    ! Read WP:
    call RHS_ALLO(NAS,NISP,lg_EP)
    call RHS_READ(NAS,NISP,lg_EP,iCASE,iSYM,iVEC)
#   ifdef _MOLCAS_MPP_
    if (iParRHS == 2) then
      myRank = GA_NodeID()
      call GA_Distribution(lg_EP,myRank,ILOV,IHIV,JLOV,JHIV)
      if (JLOV > 0) call GA_Access(lg_EP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
    end if
#   endif

    NBXSZA = NSECBX
    NBXSZJ = NINABX

    do IASTA=1,NA,NBXSZA
      IAEND = min(IASTA-1+NBXSZA,NA)
      NASZ = IAEND-IASTA+1
      do IJSTA=1,NJ,NBXSZJ
        IJEND = min(IJSTA-1+NBXSZJ,NJ)
        NJSZ = IJEND-IJSTA+1

        IAJSTA = 1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
        call DGEMM_('N','T',NV*NL,NASZ*NJSZ,NCHO,One,Cho_Ket,NV*NL,Cho_Bra(IAJSTA,1),NA*NJ,Zero,AJVL,NV*NL)
        if (iParRHS == 1) then
          IAJ = 0
          IBUF = 0
          do IJ=IJSTA,IJEND
            IJABS = IJ+NIES(ISYJ)
            do IA=IASTA,IAEND
              IAJ = IAJ+1

              do IV=1,NV
                do IL=1,NL
                  ILABS = IL+NIES(ISYL)
                  SCL = sqrt(Half)
                  if (IJABS >= ILABS) then
                    JGEL = KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
                    if (IJABS == ILABS) SCL = One
                  else
                    JGEL = KIGEJ(ILABS,IJABS)-NIGEJES(ISYJL)
                  end if
                  IW1 = IV
                  IW2 = IA+NA*(JGEL-1)+IOFF1(ISYA)
                  IW = IW1+NAS*(IW2-1)
                  !WP(v,a,jl) = ((ajvl)+(alvj))/sqrt(2+2*Kron(jl))
                  !Buff(LWEP-1+IWEP) = Buff(LWEP-1+IWEP)+SCL*AJVL(IV,IL,IAJ)
                  IBUF = IBUF+1
                  idxBuf(IBUF) = IW
                  Buff(IBUF) = SCL*AJVL(IV,IL,IAJ)
                  if (IBUF == NBUFF) then
                    call RHS_SCATTER(LDEP,lg_EP,Buff,idxBuf,IBUF)
                    IBUF = 0
                  end if
                end do
              end do
            end do
          end do
          if (IBUF /= 0) call RHS_SCATTER(LDEP,lg_EP,Buff,idxBuf,IBUF)
#       ifdef _MOLCAS_MPP_
        else
          call GADSUM_ADDRHS(AJVL,NV*NL*NASZ*NJSZ)
          if (JLOV > 0) then
            IAJ = 0
            do IJ=IJSTA,IJEND
              IJABS = IJ+NIES(ISYJ)
              do IA=IASTA,IAEND
                IAJ = IAJ+1

                do IV=ILOV,IHIV
                  IW1 = IV
                  do IL=1,NL
                    ILABS = IL+NIES(ISYL)
                    SCL = sqrt(Half)
                    if (IJABS >= ILABS) then
                      JGEL = KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
                      if (IJABS == ILABS) SCL = One
                    else
                      JGEL = KIGEJ(ILABS,IJABS)-NIGEJES(ISYJL)
                    end if
                    IW2 = IA+NA*(JGEL-1)+IOFF1(ISYA)
                    if ((IW2 >= JLOV) .and. (IW2 <= JHIV)) &
                      DBL_MB(MV+IW1-ILOV+LDEP*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDEP*(IW2-JLOV))+SCL*AJVL(IV,IL,IAJ)
                  end do
                end do
              end do
            end do
          end if
#       endif
        end if

      end do
    end do

#   ifdef _MOLCAS_MPP_
    if ((iParRHS == 2) .and. (JLOV > 0)) call GA_Release_Update(lg_EP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
#   endif
    call RHS_SAVE(NAS,NISP,lg_EP,iCASE,iSYM,iVEC)
    call RHS_FREE(lg_EP)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! The minus combination:
  if (NWM > 0) then
    ICASE = 7
    ! Read WM:
    call RHS_ALLO(NAS,NISM,lg_EM)
    call RHS_READ(NAS,NISM,lg_EM,iCASE,iSYM,iVEC)
#   ifdef _MOLCAS_MPP_
    if (iParRHS == 2) then
      myRank = GA_NodeID()
      call GA_Distribution(lg_EM,myRank,ILOV,IHIV,JLOV,JHIV)
      if (JLOV > 0) call GA_Access(lg_EM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
    end if
#   endif

    NBXSZA = NSECBX
    NBXSZJ = NINABX

    do IASTA=1,NA,NBXSZA
      IAEND = min(IASTA-1+NBXSZA,NA)
      NASZ = IAEND-IASTA+1
      do IJSTA=1,NJ,NBXSZJ
        IJEND = min(IJSTA-1+NBXSZJ,NJ)
        NJSZ = IJEND-IJSTA+1

        IAJSTA = 1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
        call DGEMM_('N','T',NV*NL,NASZ*NJSZ,NCHO,One,Cho_Ket,NV*NL,Cho_Bra(IAJSTA,1),NA*NJ,Zero,AJVL,NV*NL)

        if (iParRHS == 1) then
          IAJ = 0
          IBUF = 0
          do IJ=IJSTA,IJEND
            IJABS = IJ+NIES(ISYJ)
            do IA=IASTA,IAEND
              IAJ = IAJ+1

              do IV=1,NV
                do IL=1,NL
                  ILABS = IL+NIES(ISYL)
                  if (IJABS /= ILABS) then
                    if (IJABS > ILABS) then
                      SCL = SQ32
                      JGTL = KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
                    else
                      SCL = -SQ32
                      JGTL = KIGTJ(ILABS,IJABS)-NIGTJES(ISYJL)
                    end if
                    IW1 = IV
                    IW2 = IA+NA*(JGTL-1)+IOFF2(ISYA)
                    IW = IW1+NAS*(IW2-1)
                    !Buff(LWEM-1+IWEM) = Buff(LWEM-1+IWEM)+SCL*AJVL(IV,IL,IAJ)
                    IBUF = IBUF+1
                    idxBuf(IBUF) = IW
                    Buff(IBUF) = SCL*AJVL(IV,IL,IAJ)
                    if (IBUF == NBUFF) then
                      call RHS_SCATTER(LDEM,lg_EM,Buff,idxBuf,IBUF)
                      IBUF = 0
                    end if
                  end if
                end do
              end do
            end do
          end do
          if (IBUF /= 0) call RHS_SCATTER(LDEM,lg_EM,Buff,idxBuf,IBUF)
#       ifdef _MOLCAS_MPP_
        else if (iParRHS == 2) then
          call GADSUM_ADDRHS(AJVL,NV*NL*NASZ*NJSZ)
          if (JLOV > 0) then
            IAJ = 0
            do IJ=IJSTA,IJEND
              IJABS = IJ+NIES(ISYJ)
              do IA=IASTA,IAEND
                IAJ = IAJ+1

                do IV=ILOV,IHIV
                  IW1 = IV
                  do IL=1,NL
                    ILABS = IL+NIES(ISYL)
                    if (IJABS /= ILABS) then
                      if (IJABS > ILABS) then
                        SCL = SQ32
                        JGTL = KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
                      else
                        SCL = -SQ32
                        JGTL = KIGTJ(ILABS,IJABS)-NIGTJES(ISYJL)
                      end if
                      IW2 = IA+NA*(JGTL-1)+IOFF2(ISYA)
                      if ((IW2 >= JLOV) .and. (IW2 <= JHIV)) &
                        DBL_MB(MV+IW1-ILOV+LDEM*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDEM*(IW2-JLOV))+SCL*AJVL(IV,IL,IAJ)
                    end if
                  end do
                end do
              end do
            end do
          end if
#       endif
        end if

      end do
    end do

#   ifdef _MOLCAS_MPP_
    if ((iParRHS == 2) .and. (JLOV > 0)) call GA_Release_Update(lg_EM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
#   endif
    call RHS_SAVE(NAS,NISM,lg_EM,iCASE,iSYM,iVEC)
    call RHS_FREE(lg_EM)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine ADDRHSE

subroutine ADDRHSF(IVEC,JSYM,ISYU,ISYX,NA,NU,NC,NX,AUCX,nBuff,Buff,idxBuf,Cho_Bra,Cho_Ket,NCHO)
  ! Case F:
  !   WP(ux,ac)=((aucx)+(axcu))*(1-Kron(x,u)/2) /2
  ! With new normalisation, replace /2 with /(2*SQRT(1+Kron(ac))
  !   WM(ux,ac)= -((aucx)-(axcu))/2

  use SUPERINDEX, only: KAGEB, KAGTB, KTGEU, KTGTU
  use caspt2_module, only: NAES, NAGEB, NAGEBES, NAGTB, NAGTBES, NINDEP, NSES, NTGEU, NTGEUES, NTGTU, NTGTUES
# ifdef _MOLCAS_MPP_
  use Para_Info, only: Is_Real_Par
# endif

  integer(kind=iwp), intent(in) :: IVEC, JSYM, ISYU, ISYX, NA, NU, NC, NX, nBuff, NCHO
  real(kind=wp), intent(out) :: AUCX(NA,NU,NC,NX), Buff(nBuff)
  integer(kind=iwp), intent(out) :: idxBuf(nBuff)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NC,NX,NCHO)
  integer(kind=iwp) :: IA, IAABS, IBUF, IC, ICABS, ICASE, ISYA, ISYC, ISYM, IU, IUABS, IW, IW1, IW2, IX, IXABS, IXMAX, LDFM, LDFP, &
                       lg_FM, lg_FP, NASM, NASP, NISM, NISP, NWFM, NWFP
  real(kind=wp) :: SCL, SCL1
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHIV, ILOV, JHIV, JLOV, LDV, MV, myRank
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *

  if (ISYU < ISYX) return

  ISYA = Mul(JSYM,ISYU)
  ISYC = Mul(JSYM,ISYX)
  ISYM = Mul(ISYU,ISYX)

  if (NINDEP(ISYM,8) > 0 .or. Do_SC) then
    ! The plus combination:
    NASP = NTGEU(ISYM)
    NISP = NAGEB(ISYM)
    NWFP = NASP*NISP
  else
    NWFP = 0
  end if
  if (NINDEP(ISYM,9) > 0 .or. Do_SC) then
    ICASE = 9
    ! The minus combination:
    NASM = NTGTU(ISYM)
    NISM = NAGTB(ISYM)
    NWFM = NASM*NISM
  else
    NWFM = 0
  end if
  if (NWFP+NWFM <= 0) return
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Incore or partitioned option?

  !Incore = nBuff >= max(NWFP,NWFM)+NA*NU*NC*NX
  !if (.not. Incore) then
  !  write(u6,*) 'Sort out of memory in ADDRHSF'
  !  call Abend()
  !end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call DGEMM_('N','T',NA*NU,NC*NX,NCHO,One,Cho_Bra,NA*NU,Cho_Ket,NC*NX,Zero,AUCX,NA*NU)
# ifdef _MOLCAS_MPP_
  if (iParRHS == 2) call GADSUM_ADDRHS(AUCX,NA*NU*NC*NX)
# endif

  if (NWFP > 0) then
    if (NINDEP(ISYM,8) > 0 .or. Do_SC) then
      ! The plus combination:
      ICASE = 8
      NASP = NTGEU(ISYM)
      NISP = NAGEB(ISYM)
      NWFP = NASP*NISP
      !LWFP = 1+NA*NU*NC*NX
      LDFP = NASP
      ! Read WP:
      call RHS_ALLO(NASP,NISP,lg_FP)
      call RHS_READ(NASP,NISP,lg_FP,iCASE,iSYM,iVEC)

      if (iParRHS == 1) then
        IBUF = 0
        do IU=1,NU
          IUABS = IU+NAES(ISYU)
          IXMAX = NX
          if (ISYU == ISYX) IXMAX = IU
          do IX=1,IXMAX
            IXABS = IX+NAES(ISYX)
            SCL1 = Half
            if (IUABS == IXABS) SCL1 = Quart
            IW1 = KTGEU(IUABS,IXABS)-NTGEUES(ISYM)
            do IA=1,NA
              IAABS = IA+NSES(ISYA)
              do IC=1,NC
                ICABS = IC+NSES(ISYC)
                SCL = SCL1
                if (IAABS >= ICABS) then
                  IW2 = KAGEB(IAABS,ICABS)-NAGEBES(ISYM)
                  if (IAABS == ICABS) SCL = sqrt(Two)*SCL1
                else
                  IW2 = KAGEB(ICABS,IAABS)-NAGEBES(ISYM)
                end if
                IW = IW1+NASP*(IW2-1)
                !WFP = Buff(LWFP-1+IW)+SCL*AUCX(IA,IU,IC,IX)
                !Buff(LWFP-1+IW) = WFP
                IBUF = IBUF+1
                idxBuf(IBUF) = IW
                Buff(IBUF) = SCL*AUCX(IA,IU,IC,IX)
                if (IBUF == NBUFF) then
                  call RHS_SCATTER(LDFP,lg_FP,Buff,idxBuf,IBUF)
                  IBUF = 0
                end if
              end do
            end do
          end do
        end do
        if (IBUF /= 0) call RHS_SCATTER(LDFP,lg_FP,Buff,idxBuf,IBUF)
#     ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
        myRank = GA_NodeID()
        call GA_Distribution(lg_FP,myRank,ILOV,IHIV,JLOV,JHIV)
        if (JLOV > 0) then
          call GA_Access(lg_FP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
          IBUF = 0
          do IU=1,NU
            IUABS = IU+NAES(ISYU)
            IXMAX = NX
            if (ISYU == ISYX) IXMAX = IU
            do IX=1,IXMAX
              IXABS = IX+NAES(ISYX)
              SCL1 = Half
              if (IUABS == IXABS) SCL1 = Quart
              IW1 = KTGEU(IUABS,IXABS)-NTGEUES(ISYM)
              if ((IW1 < ILOV) .or. (IW1 > IHIV)) cycle
              do IA=1,NA
                IAABS = IA+NSES(ISYA)
                do IC=1,NC
                  ICABS = IC+NSES(ISYC)
                  SCL = SCL1
                  if (IAABS >= ICABS) then
                    IW2 = KAGEB(IAABS,ICABS)-NAGEBES(ISYM)
                    if (IAABS == ICABS) SCL = sqrt(Two)*SCL1
                  else
                    IW2 = KAGEB(ICABS,IAABS)-NAGEBES(ISYM)
                  end if
                  if ((IW2 >= JLOV) .and. (IW2 <= JHIV)) &
                    DBL_MB(MV+IW1-ILOV+LDFP*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDFP*(IW2-JLOV))+SCL*AUCX(IA,IU,IC,IX)
                end do
              end do
            end do
          end do
          call GA_Release_Update(lg_FP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
        end if
#     endif
      end if

      call RHS_SAVE(NASP,NISP,lg_FP,iCASE,iSYM,iVEC)
      call RHS_FREE(lg_FP)
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (NWFM <= 0) return
  if (NINDEP(ISYM,9) > 0 .or. Do_SC) then
    ICASE = 9
    ! The minus combination:
    NASM = NTGTU(ISYM)
    NISM = NAGTB(ISYM)
    NWFM = NASM*NISM
    !LWFM = 1+NA*NU*NC*NX
    LDFM = NASM
    ! Read WM:
    call RHS_ALLO(NASM,NISM,lg_FM)
    call RHS_READ(NASM,NISM,lg_FM,iCASE,iSYM,iVEC)
#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      myRank = GA_NodeID()
      call GA_Distribution(lg_FM,myRank,ILOV,IHIV,JLOV,JHIV)
      if (JLOV > 0) call GA_Access(lg_FM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
    end if
#   endif

    if (iParRHS == 1) then
      IBUF = 0
      do IU=1,NU
        IUABS = IU+NAES(ISYU)
        IXMAX = NX
        if (ISYU == ISYX) IXMAX = IU-1
        do IX=1,IXMAX
          IXABS = IX+NAES(ISYX)
          IW1 = KTGTU(IUABS,IXABS)-NTGTUES(ISYM)
          do IA=1,NA
            IAABS = IA+NSES(ISYA)
            do IC=1,NC
              ICABS = IC+NSES(ISYC)
              if (IAABS > ICABS) then
                IW2 = KAGTB(IAABS,ICABS)-NAGTBES(ISYM)
                IW = IW1+NASM*(IW2-1)
                !Buff(LWFM-1+IW) = Buff(LWFM-1+IW)-Half*AUCX(IA,IU,IC,IX)
                IBUF = IBUF+1
                idxBuf(IBUF) = IW
                Buff(IBUF) = -Half*AUCX(IA,IU,IC,IX)
              else if (IAABS < ICABS) then
                IW2 = KAGTB(ICABS,IAABS)-NAGTBES(ISYM)
                IW = IW1+NASM*(IW2-1)
                !Buff(LWFM-1+IW) = Buff(LWFM-1+IW)+Half*AUCX(IA,IU,IC,IX)
                IBUF = IBUF+1
                idxBuf(IBUF) = IW
                Buff(IBUF) = Half*AUCX(IA,IU,IC,IX)
              end if
              if (IBUF == NBUFF) then
                call RHS_SCATTER(LDFM,lg_FM,Buff,idxBuf,IBUF)
                IBUF = 0
              end if
            end do
          end do
        end do
      end do
      if (IBUF /= 0) call RHS_SCATTER(LDFM,lg_FM,Buff,idxBuf,IBUF)
#   ifdef _MOLCAS_MPP_
    else if (iParRHS == 2) then
      if (JLOV > 0) then
        do IU=1,NU
          IUABS = IU+NAES(ISYU)
          IXMAX = NX
          if (ISYU == ISYX) IXMAX = IU-1
          do IX=1,IXMAX
            IXABS = IX+NAES(ISYX)
            IW1 = KTGTU(IUABS,IXABS)-NTGTUES(ISYM)
            if ((IW1 < ILOV) .or. (IW1 > IHIV)) cycle
            do IA=1,NA
              IAABS = IA+NSES(ISYA)
              do IC=1,NC
                ICABS = IC+NSES(ISYC)
                if (IAABS > ICABS) then
                  IW2 = KAGTB(IAABS,ICABS)-NAGTBES(ISYM)
                  if ((IW2 >= JLOV) .and. (IW2 <= JHIV)) &
                    DBL_MB(MV+IW1-ILOV+LDFM*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDFM*(IW2-JLOV))-Half*AUCX(IA,IU,IC,IX)
                else if (IAABS < ICABS) then
                  IW2 = KAGTB(ICABS,IAABS)-NAGTBES(ISYM)
                  if ((IW2 >= JLOV) .and. (IW2 <= JHIV)) &
                    DBL_MB(MV+IW1-ILOV+LDFM*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDFM*(IW2-JLOV))+Half*AUCX(IA,IU,IC,IX)
                end if
              end do
            end do
          end do
        end do
        call GA_Release_Update(lg_FM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
      end if
#   endif
    end if

    ! Put WFM on disk:
    call RHS_SAVE(NASM,NISM,lg_FM,iCASE,iSYM,iVEC)
    call RHS_FREE(lg_FM)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine ADDRHSF

subroutine ADDRHSG(IVEC,JSYM,ISYU,ISYL,NA,NU,NC,NL,AUCL,NAUCL,nBuff,Buff,idxBuf,Cho_Bra,Cho_Ket,NCHO)
  ! Case G:
  !   WP(u,l,ac)=  ((aucl)+cual))/SQRT(2+2*Kron(ab))
  !   WM(u,l,ac)=  ((aucl)-cual))*SQRT(OneHalf)

  use SUPERINDEX, only: KAGEB, KAGTB
  use caspt2_module, only: NAGEB, NAGEBES, NAGTB, NAGTBES, NASH, NISH, NISUP, NSECBX, NSES, NSYM

  integer(kind=iwp), intent(in) :: IVEC, JSYM, ISYU, ISYL, NA, NU, NC, NL, NAUCL, nBuff, NCHO
  real(kind=wp), intent(out) :: AUCL(NA,NU,*), Buff(nBuff)
  integer(kind=iwp), intent(out) :: idxBuf(nBuff)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NC*NL,NCHO)
  integer(kind=iwp) :: IA, IAABS, IAGEC, IAGTC, IBUF, IC, ICABS, ICASE, ICEND, ICL, ICLSTA, ICSTA, IL, ILEND, ILSTA, IO1, IO2, &
                       IOFF1(8), IOFF2(8), ISAB, ISI, ISYA, ISYAC, ISYC, ISYM, IU, IW, IW1, IW2, KCL, LDGM, LDGP, lg_GM, lg_GP, &
                       NAS, NBXSZC, NBXSZL, NCSZ, NISM, NISP, NLSZ, NWGM, NWGP
  real(kind=wp) :: SCL
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHIV, ILOV, ITMP1, ITMP2, JHIV, JLOV, LDV, MV, myRank
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *

  ISYA = Mul(JSYM,ISYU)
  ISYC = Mul(JSYM,ISYL)
  ISYM = ISYU
  ISYAC = Mul(ISYA,ISYC)
  ! Set up offset table:
  IO1 = 0
  IO2 = 0
  do ISI=1,NSYM
    IOFF1(ISI) = IO1
    IOFF2(ISI) = IO2
    ISAB = Mul(ISI,ISYM)
    IO1 = IO1+NISH(ISI)*NAGEB(ISAB)
    IO2 = IO2+NISH(ISI)*NAGTB(ISAB)
  end do

  ! Allocate W with parts WP,WM
  NAS = NASH(ISYM)
  NISP = NISUP(ISYM,10)
  NISM = NISUP(ISYM,11)
  NWGP = NAS*NISP
  NWGM = NAS*NISM
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Incore or partitioned option?

  !Incore = nBuff >= max(NWGP,NWGM)+NA*NU*NC*NL
  !if (.not. Incore) then
  !  write(u6,*) 'Sort out of memory in ADDRHSG'
  !  call Abend()
  !end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !call DGEMM_('N','T',NA*NU,NC*NL,NCHO,One,Cho_Bra,NA*NU,Cho_Ket,NC*NL,Zero,AUCL,NA*NU)

  !LWG = 1+NA*NU*NC*NL
  !LWGP = LWG
  !LWGM = LWG
  LDGP = NAS
  LDGM = NAS

  ! The plus combination:
  if (NWGP > 0) then
    ICASE = 10
    ! Read WP:
    call RHS_ALLO(NAS,NISP,lg_GP)
    call RHS_READ(NAS,NISP,lg_GP,iCASE,iSYM,iVEC)
#   ifdef _MOLCAS_MPP_
    if (iParRHS == 2) then
      myRank = GA_NodeID()
      call GA_Distribution(lg_GP,myRank,ILOV,IHIV,JLOV,JHIV)
      if (JLOV > 0) call GA_Access(lg_GP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
    end if
#   endif

    ! to keep memory advantage, scale NL so that NC*NL fits in KCL
    ! (scaling NL does not affect ordering of integrals = safe)
    NBXSZC = NSECBX
    !NBXSZJ = NINABX
    KCL = NAUCL/(NA*NU)
    NBXSZL = KCL/NC
    if (NBXSZL <= 0) then
      write(u6,*) 'Not enough memory in ADDRHSG, I give up'
      call Abend()
    end if

    do ICSTA=1,NC,NBXSZC
      ICEND = min(ICSTA-1+NBXSZC,NC)
      NCSZ = ICEND-ICSTA+1
      do ILSTA=1,NL,NBXSZL
        ILEND = min(ILSTA-1+NBXSZL,NL)
        NLSZ = ILEND-ILSTA+1

        ICLSTA = 1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
        call DGEMM_('N','T',NA*NU,NCSZ*NLSZ,NCHO,One,Cho_Bra,NA*NU,Cho_Ket(ICLSTA,1),NC*NL,Zero,AUCL,NA*NU)

        if (iParRHS == 1) then
          ICL = 0
          IBUF = 0
          do IL=ILSTA,ILEND
            do IC=ICSTA,ICEND
              ICABS = IC+NSES(ISYC)
              ICL = ICL+1

              do IA=1,NA
                IAABS = IA+NSES(ISYA)
                SCL = sqrt(Half)
                if (IAABS >= ICABS) then
                  IAGEC = KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                  if (IAABS == ICABS) SCL = One
                else
                  IAGEC = KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
                end if
                do IU=1,NU
                  IW1 = IU
                  IW2 = IL+NL*(IAGEC-1)+IOFF1(ISYL)
                  IW = IW1+NAS*(IW2-1)
                  IBUF = IBUF+1
                  !Buff(LWGP-1+IW) = Buff(LWGP-1+IW)+SCL*AUCL(IA,IU,ICL)
                  idxBuf(IBUF) = IW
                  Buff(IBUF) = SCL*AUCL(IA,IU,ICL)
                  if (IBUF == NBUFF) then
                    call RHS_SCATTER(LDGP,lg_GP,Buff,idxBuf,IBUF)
                    IBUF = 0
                  end if
                end do
              end do
            end do
          end do
          if (IBUF /= 0) call RHS_SCATTER(LDGP,lg_GP,Buff,idxBuf,IBUF)
#       ifdef _MOLCAS_MPP_
        else if (iParRHS == 2) then
          call GADSUM_ADDRHS(AUCL,NA*NU*NCSZ*NLSZ)
          if (JLOV > 0) then
            ICL = 0
            do IL=ILSTA,ILEND
              do IC=ICSTA,ICEND
                ICABS = IC+NSES(ISYC)
                ICL = ICL+1

                do IA=1,NA
                  IAABS = IA+NSES(ISYA)
                  SCL = sqrt(Half)
                  if (IAABS >= ICABS) then
                    IAGEC = KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                    if (IAABS == ICABS) SCL = One
                  else
                    IAGEC = KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
                  end if
                  IW2 = IL+NL*(IAGEC-1)+IOFF1(ISYL)
                  if ((IW2 < JLOV) .or. (IW2 > JHIV)) cycle
                  do IU=1,NU
                    IW1 = IU
                    if ((IW1 >= ILOV) .and. (IW1 <= IHIV)) &
                      DBL_MB(MV+IW1-ILOV+LDGP*(IW2-JLOV)) = DBL_MB(MV+IW1-ILOV+LDGP*(IW2-JLOV))+SCL*AUCL(IA,IU,ICL)
                  end do
                end do
              end do
            end do
          end if
#       endif
        end if

      end do
    end do

#   ifdef _MOLCAS_MPP_
    if ((iParRHS == 2) .and. (JLOV > 0)) call GA_Release_Update(lg_GP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
#   endif
    call RHS_SAVE(NAS,NISP,lg_GP,iCASE,iSYM,iVEC)
    call RHS_FREE(lg_GP)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! The minus combination:
  if (NWGM > 0) then
    ICASE = 11
    ! Read WGM:
    call RHS_ALLO(NAS,NISM,lg_GM)
    call RHS_READ(NAS,NISM,lg_GM,iCASE,iSYM,iVEC)
#   ifdef _MOLCAS_MPP_
    if (iParRHS == 2) then
      myRank = GA_NodeID()
      call GA_Distribution(lg_GM,myRank,ILOV,IHIV,JLOV,JHIV)
      if (JLOV > 0) call GA_Access(lg_GM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
    end if
#   endif

    ! to keep memory advantage, scale NL so that NC*NL fits in KCL
    ! (scaling NL does not affect ordering of integrals = safe)
    NBXSZC = NSECBX
    !NBXSZJ = NINABX
    KCL = NAUCL/(NA*NU)
    NBXSZL = KCL/NC
    if (NBXSZL <= 0) then
      write(u6,*) 'Not enough memory in ADDRHSG, I give up'
      call Abend()
    end if

    do ICSTA=1,NC,NBXSZC
      ICEND = min(ICSTA-1+NBXSZC,NC)
      NCSZ = ICEND-ICSTA+1
      do ILSTA=1,NL,NBXSZL
        ILEND = min(ILSTA-1+NBXSZL,NL)
        NLSZ = ILEND-ILSTA+1

        ICLSTA = 1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
        call DGEMM_('N','T',NA*NU,NCSZ*NLSZ,NCHO,One,Cho_Bra,NA*NU,Cho_Ket(ICLSTA,1),NC*NL,Zero,AUCL,NA*NU)

        if (iParRHS == 1) then
          ICL = 0
          IBUF = 0
          do IL=ILSTA,ILEND
            do IC=ICSTA,ICEND
              ICABS = IC+NSES(ISYC)
              ICL = ICL+1

              do IA=1,NA
                IAABS = IA+NSES(ISYA)
                if (IAABS > ICABS) then
                  IAGTC = KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                  SCL = sqrt(OneHalf)
                  do IU=1,NU
                    IW1 = IU
                    IW2 = IL+NL*(IAGTC-1)+IOFF2(ISYL)
                    IW = IW1+NAS*(IW2-1)
                    !Buff(LWGM-1+IW) = Buff(LWGM-1+IW)+SCL*AUCL(IA,IU,ICL)
                    IBUF = IBUF+1
                    idxBuf(IBUF) = IW
                    Buff(IBUF) = SCL*AUCL(IA,IU,ICL)
                    if (IBUF == NBUFF) then
                      call RHS_SCATTER(LDGM,lg_GM,Buff,idxBuf,IBUF)
                      IBUF = 0
                    end if
                  end do
                else if (IAABS < ICABS) then
                  IAGTC = KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                  SCL = -sqrt(OneHalf)
                  do IU=1,NU
                    IW1 = IU
                    IW2 = IL+NL*(IAGTC-1)+IOFF2(ISYL)
                    IW = IW1+NAS*(IW2-1)
                    !Buff(LWGM-1+IW) = Buff(LWGM-1+IW)+SCL*AUCL(IA,IU,ICL)
                    IBUF = IBUF+1
                    idxBuf(IBUF) = IW
                    Buff(IBUF) = SCL*AUCL(IA,IU,ICL)
                    if (IBUF == NBUFF) then
                      call RHS_SCATTER(LDGM,lg_GM,Buff,idxBuf,IBUF)
                      IBUF = 0
                    end if
                  end do
                end if
              end do
            end do
          end do
          if (IBUF /= 0) call RHS_SCATTER(LDGM,lg_GM,Buff,idxBuf,IBUF)
#       ifdef _MOLCAS_MPP_
        else if (iParRHS == 2) then
          call GADSUM_ADDRHS(AUCL,NA*NU*NCSZ*NLSZ)
          if (JLOV > 0) then
            ICL = 0
            do IL=ILSTA,ILEND
              do IC=ICSTA,ICEND
                ICABS = IC+NSES(ISYC)
                ICL = ICL+1

                do IA=1,NA
                  IAABS = IA+NSES(ISYA)
                  if (IAABS > ICABS) then
                    IAGTC = KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                    ITMP2 = IL+NL*(IAGTC-1)+IOFF2(ISYL)
                    if ((ITMP2 < JLOV) .or. (ITMP2 > JHIV)) cycle
                    SCL = sqrt(OneHalf)
                    do IU=1,NU
                      ITMP1 = IU
                      if ((ITMP1 >= ILOV) .and. (ITMP1 <= IHIV)) &
                        DBL_MB(MV+ITMP1-ILOV+LDGP*(ITMP2-JLOV)) = DBL_MB(MV+ITMP1-ILOV+LDGP*(ITMP2-JLOV))+SCL*AUCL(IA,IU,ICL)
                    end do
                  else if (IAABS < ICABS) then
                    IAGTC = KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                    ITMP2 = IL+NL*(IAGTC-1)+IOFF2(ISYL)
                    if ((ITMP2 < JLOV) .or. (ITMP2 > JHIV)) cycle
                    SCL = -sqrt(OneHalf)
                    do IU=1,NU
                      ITMP1 = IU
                      if ((ITMP1 >= ILOV) .and. (ITMP1 <= IHIV)) &
                        DBL_MB(MV+ITMP1-ILOV+LDGP*(ITMP2-JLOV)) = DBL_MB(MV+ITMP1-ILOV+LDGP*(ITMP2-JLOV))+SCL*AUCL(IA,IU,ICL)
                    end do
                  end if
                end do
              end do
            end do
          end if
#       endif
        end if

      end do
    end do

#   ifdef _MOLCAS_MPP_
    if ((iParRHS == 2) .and. (JLOV > 0)) call GA_Release_Update(lg_GM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
#   endif
    call RHS_SAVE(NAS,NISM,lg_GM,iCASE,iSYM,iVEC)
    call RHS_FREE(lg_GM)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine ADDRHSG

subroutine ADDRHSH(IVEC,JSYM,ISYJ,ISYL,NA,NJ,NC,NL,AJCL,NAJCL,nBuff,Buff,idxBuf,Cho_Bra,Cho_Ket,NCHO)
  ! Case H:
  !   WP(jl,ac)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
  !   WM(jl,ac)=((ajcl)-(alcj))*SQRT(Three)

  use SUPERINDEX, only: KAGEB, KAGTB, KIGEJ, KIGTJ
  use caspt2_module, only: NAGEB, NAGEBES, NAGTB, NAGTBES, NIES, NIGEJ, NIGEJES, NIGTJ, NIGTJES, NINABX, NSECBX, NSES

  integer(kind=iwp), intent(in) :: IVEC, JSYM, ISYJ, ISYL, NA, NJ, NC, NL, NAJCL, nBuff, NCHO
  real(kind=wp), intent(out) :: AJCL(NC*NL,*), Buff(nBuff)
  integer(kind=iwp), intent(out) :: idxBuf(nBuff)
  real(kind=wp), intent(in) :: Cho_Bra(NA*NJ,NCHO), Cho_Ket(NC*NL,NCHO)
  integer(kind=iwp) :: IA, IAABS, IAEND, IAGEC, IAGTC, IAJ, IAJSTA, IASTA, IBUF, IC, ICABS, ICASE, ICEND, ICL, ICLSTA, ICSTA, IJ, &
                       IJABS, IJEND, IJGEL, IJGTL, IJSTA, IL, ILABS, ILEND, ILMAX, ILSTA, ISYA, ISYAC, ISYC, ISYJL, ISYM, IW, KAJ, &
                       LDHM, LDHP, lg_HM, lg_HP, NASM, NASP, NASZ, NBXSZA, NBXSZC, NBXSZJ, NBXSZL, NCSZ, NISM, NISP, NJSZ, NWHM, &
                       NWHP
  real(kind=wp) :: SCL, SCL1
# ifdef _MOLCAS_MPP_
  integer(kind=iwp) :: IHIV, ILOV, ITMP1, ITMP2, JHIV, JLOV, LDV, MV, myRank
# endif

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  if (ISYJ < ISYL) return
  ISYA = Mul(JSYM,ISYJ)
  ISYC = Mul(JSYM,ISYL)
  ISYM = Mul(ISYA,ISYC)
  ISYAC = ISYM
  ISYJL = ISYM

  ! Allocate WHP,WHM
  NASP = NAGEB(ISYM)
  NISP = NIGEJ(ISYM)
  NWHP = NASP*NISP
  if (NWHP == 0) return
  NASM = NAGTB(ISYM)
  NISM = NIGTJ(ISYM)
  NWHM = NASM*NISM
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Incore or partitioned option?

  !Incore = nBuff >= max(NWHP,NWHM)+NC*NL
  !if (.not. Incore) then
  !   write(u6,*) 'Sort out of memory in ADDRHSH'
  !   call Abend()
  !end if

  !nBatch = min((nBuff-max(NWHP,NWHM))/(NC*NL),NA*NJ)
  !if ((iPrGlb >= DEBUG) .and. (nBatch < NA*NJ)) then
  !  write(u6,'(1X,A)') 'less memory than ideal for ADDRHSH:'
  !  write(u6,'(1X,A12,I12)') 'needed    = ',NA*NJ*NC*NL
  !  write(u6,'(1X,A12,I12)') 'available = ',nBatch*NC*NL
  !end if

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !call DGEMM_('N','T',NA*NJ,NC*NL,NCHO,One,Cho_Bra,NA*NJ,Cho_Ket,NC*NL,Zero,AJCL,NA*NJ)
  !LWHP = 1+nBatch*NC*NL
  !LWHM = LWHP
  LDHP = NASP
  LDHM = NASM

  ! The plus combination:
  if (NWHP > 0) then
    ICASE = 12
    ! Read WP:
    call RHS_ALLO(NASP,NISP,lg_HP)
    call RHS_READ(NASP,NISP,lg_HP,iCASE,iSYM,iVEC)

#   ifdef _MOLCAS_MPP_
    if (iParRHS == 2) then
      myRank = GA_NodeID()
      call GA_Distribution(lg_HP,myRank,ILOV,IHIV,JLOV,JHIV)
      if (JLOV > 0) call GA_Access(lg_HP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
    end if
#   endif

    ! WP(jl,ac)=((ajcl)+(alcj))/SQRT((1+Kron(jl))*(1+Kron(ac))
    !-SVC20100313: use only part of aj but all of cl, this introduces
    ! a reordered loop over parts of the RHS array.

    ! to keep memory advantage, scale NJ so that NA*NJ fits in KAJ
    ! (scaling NJ or NL does not affect ordering of integrals = safe)
    NBXSZA = NSECBX
    !NBXSZJ = NINABX
    KAJ = NAJCL/(NC*NL)
    NBXSZJ = KAJ/NA
    if (NBXSZJ <= 0) then
      write(u6,*) 'Not enough memory in ADDRHSH, I give up'
      call Abend()
    end if

    NBXSZC = NSECBX
    NBXSZL = NINABX

    do IASTA=1,NA,NBXSZA
      IAEND = min(IASTA-1+NBXSZA,NA)
      NASZ = IAEND-IASTA+1
      do IJSTA=1,NJ,NBXSZJ
        IJEND = min(IJSTA-1+NBXSZJ,NJ)
        NJSZ = IJEND-IJSTA+1

        IAJSTA = 1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
        call DGEMM_('N','T',NC*NL,NASZ*NJSZ,NCHO,One,Cho_Ket,NC*NL,Cho_Bra(IAJSTA,1),NA*NJ,Zero,AJCL,NC*NL)

        if (iParRHS == 1) then
          do ICSTA=1,NC,NBXSZC
            ICEND = min(ICSTA-1+NBXSZC,NC)
            NCSZ = ICEND-ICSTA+1
            do ILSTA=1,NL,NBXSZL
              ILEND = min(ILSTA-1+NBXSZL,NL)

              ICLSTA = 1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

              IAJ = 0
              IBUF = 0
              do IJ=IJSTA,IJEND
                IJABS = IJ+NIES(ISYJ)
                ILMAX = NL
                if (ISYJ == ISYL) ILMAX = IJ
                do IA=IASTA,IAEND
                  IAABS = IA+NSES(ISYA)
                  IAJ = IAJ+1

                  ICL = 0
                  do IL=ILSTA,min(ILEND,ILMAX)
                    ILABS = IL+NIES(ISYL)
                    SCL1 = One
                    IJGEL = KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
                    if (IJABS == ILABS) SCL1 = sqrt(Half)
                    do IC=ICSTA,ICEND
                      ICABS = IC+NSES(ISYC)
                      ICL = ICL+1

                      SCL = SCL1
                      if (IAABS >= ICABS) then
                        IAGEC = KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                        if (IAABS == ICABS) SCL = sqrt(Two)*SCL1
                      else
                        IAGEC = KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
                      end if
                      IW = IAGEC+NAGEB(ISYM)*(IJGEL-1)
                      !Buff(LWHP-1+IW) = Buff(LWHP-1+IW)+SCL*AJCL(ICLSTA+ICL-1,IAJ)
                      IBUF = IBUF+1
                      idxBuf(IBUF) = IW
                      Buff(IBUF) = SCL*AJCL(ICLSTA+ICL-1,IAJ)
                      if (IBUF == NBUFF) then
                        call RHS_SCATTER(LDHP,lg_HP,Buff,idxBuf,IBUF)
                        IBUF = 0
                      end if
                    end do
                  end do
                end do
              end do
              if (IBUF /= 0) call RHS_SCATTER(LDHP,lg_HP,Buff,idxBuf,IBUF)

            end do
          end do
#       ifdef _MOLCAS_MPP_
        else if (iParRHS == 2) then
          call GADSUM_ADDRHS(AJCL,NC*NL*NASZ*NJSZ)
          if (JLOV > 0) then

            do ICSTA=1,NC,NBXSZC
              ICEND = min(ICSTA-1+NBXSZC,NC)
              NCSZ = ICEND-ICSTA+1
              do ILSTA=1,NL,NBXSZL
                ILEND = min(ILSTA-1+NBXSZL,NL)

                ICLSTA = 1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

                IAJ = 0
                do IJ=IJSTA,IJEND
                  IJABS = IJ+NIES(ISYJ)
                  ILMAX = NL
                  if (ISYJ == ISYL) ILMAX = IJ
                  do IA=IASTA,IAEND
                    IAABS = IA+NSES(ISYA)
                    IAJ = IAJ+1

                    ICL = 0
                    do IL=ILSTA,min(ILEND,ILMAX)
                      ILABS = IL+NIES(ISYL)
                      SCL1 = One
                      IJGEL = KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
                      if ((IJGEL < JLOV) .or. (IJGEL > JHIV)) then
                        ICL = ICL+ICEND-ICSTA+1
                        cycle
                      end if
                      if (IJABS == ILABS) SCL1 = sqrt(Half)
                      do IC=ICSTA,ICEND
                        ICABS = IC+NSES(ISYC)
                        ICL = ICL+1

                        SCL = SCL1
                        if (IAABS >= ICABS) then
                          IAGEC = KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                          if (IAABS == ICABS) SCL = sqrt(Two)*SCL1
                        else
                          IAGEC = KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
                        end if
                        if ((IAGEC >= ILOV) .and. (IAGEC <= IHIV)) then
                          ITMP1 = IAGEC
                          ITMP2 = IJGEL
                          DBL_MB(MV+ITMP1-ILOV+LDHP*(ITMP2-JLOV)) = &
                            DBL_MB(MV+ITMP1-ILOV+LDHP*(ITMP2-JLOV))+SCL*AJCL(ICLSTA+ICL-1,IAJ)
                        end if
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end if
#       endif
        end if
      end do
    end do

#   ifdef _MOLCAS_MPP_
    if ((iParRHS == 2) .and. (JLOV > 0)) call GA_Release_Update(lg_HP,ILOV,IHIV,JLOV,JHIV,MV,LDV)
#   endif
    call RHS_SAVE(NASP,NISP,lg_HP,iCASE,iSYM,iVEC)
    call RHS_FREE(lg_HP)
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! The minus combination:
  if (NWHM == 0) return
  ICASE = 13
  ! Read WM:
  call RHS_ALLO(NASM,NISM,lg_HM)
  call RHS_READ(NASM,NISM,lg_HM,iCASE,iSYM,iVEC)

# ifdef _MOLCAS_MPP_
  if (iParRHS == 2) then
    myRank = GA_NodeID()
    call GA_Distribution(lg_HM,myRank,ILOV,IHIV,JLOV,JHIV)
    if (JLOV > 0) call GA_Access(lg_HM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
  end if
# endif

  !VM(jl,ac) = ((ajcl)-(alcj))*sqrt(Three)
  NBXSZA = NSECBX
  !NBXSZJ = NINABX
  KAJ = NAJCL/(NC*NL)
  NBXSZJ = KAJ/NA
  if (NBXSZJ <= 0) then
    write(u6,*) 'Not enough memory in ADDRHSH, I give up'
    call Abend()
  end if
  NBXSZC = NSECBX
  NBXSZL = NINABX

  do IASTA=1,NA,NBXSZA
    IAEND = min(IASTA-1+NBXSZA,NA)
    NASZ = IAEND-IASTA+1
    do IJSTA=1,NJ,NBXSZJ
      IJEND = min(IJSTA-1+NBXSZJ,NJ)
      NJSZ = IJEND-IJSTA+1

      IAJSTA = 1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
      call DGEMM_('N','T',NC*NL,NASZ*NJSZ,NCHO,One,Cho_Ket,NC*NL,Cho_Bra(IAJSTA,1),NA*NJ,Zero,AJCL,NC*NL)

      if (iParRHS == 1) then
        do ICSTA=1,NC,NBXSZC
          ICEND = min(ICSTA-1+NBXSZC,NC)
          NCSZ = ICEND-ICSTA+1
          do ILSTA=1,NL,NBXSZL
            ILEND = min(ILSTA-1+NBXSZL,NL)

            ICLSTA = 1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

            IAJ = 0
            IBUF = 0
            do IJ=IJSTA,IJEND
              IJABS = IJ+NIES(ISYJ)
              ILMAX = NL
              if (ISYJ == ISYL) ILMAX = IJ-1
              do IA=IASTA,IAEND
                IAABS = IA+NSES(ISYA)
                IAJ = IAJ+1

                ICL = 0
                do IL=ILSTA,min(ILMAX,ILEND)
                  ILABS = IL+NIES(ISYL)
                  IJGTL = KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
                  do IC=ICSTA,ICEND
                    ICABS = IC+NSES(ISYC)
                    ICL = ICL+1

                    if (IAABS > ICABS) then
                      IAGTC = KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                      SCL = sqrt(Three)
                    else if (IAABS < ICABS) then
                      IAGTC = KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                      SCL = -sqrt(Three)
                    else
                      cycle
                    end if
                    IW = IAGTC+NAGTB(ISYM)*(IJGTL-1)
                    !Buff(LWHM-1+IW) = Buff(LWHM-1+IW)+SCL*AJCL(ICLSTA+ICL-1,IAJ)
                    IBUF = IBUF+1
                    idxBuf(IBUF) = IW
                    Buff(IBUF) = SCL*AJCL(ICLSTA+ICL-1,IAJ)
                    if (IBUF == NBUFF) then
                      call RHS_SCATTER(LDHM,lg_HM,Buff,idxBuf,IBUF)
                      IBUF = 0
                    end if
                  end do
                end do
              end do
            end do
            if (IBUF /= 0) call RHS_SCATTER(LDHM,lg_HM,Buff,idxBuf,IBUF)

          end do
        end do
#     ifdef _MOLCAS_MPP_
      else if (iParRHS == 2) then
        call GADSUM_ADDRHS(AJCL,NC*NL*NASZ*NJSZ)
        if (JLOV > 0) then

          do ICSTA=1,NC,NBXSZC
            ICEND = min(ICSTA-1+NBXSZC,NC)
            NCSZ = ICEND-ICSTA+1
            do ILSTA=1,NL,NBXSZL
              ILEND = min(ILSTA-1+NBXSZL,NL)

              ICLSTA = 1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

              IAJ = 0
              do IJ=IJSTA,IJEND
                IJABS = IJ+NIES(ISYJ)
                ILMAX = NL
                if (ISYJ == ISYL) ILMAX = IJ-1
                do IA=IASTA,IAEND
                  IAABS = IA+NSES(ISYA)
                  IAJ = IAJ+1

                  ICL = 0
                  do IL=ILSTA,min(ILMAX,ILEND)
                    ILABS = IL+NIES(ISYL)
                    IJGTL = KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
                    if ((IJGTL < JLOV) .or. (IJGTL > JHIV)) then
                      ICL = ICL+ICEND-ICSTA+1
                      cycle
                    end if
                    do IC=ICSTA,ICEND
                      ICABS = IC+NSES(ISYC)
                      ICL = ICL+1

                      if (IAABS > ICABS) then
                        IAGTC = KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                        SCL = sqrt(Three)
                      else if (IAABS < ICABS) then
                        IAGTC = KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                        SCL = -sqrt(Three)
                      else
                        cycle
                      end if
                      if ((IAGTC >= ILOV) .and. (IAGTC <= IHIV)) then
                        ITMP1 = IAGTC
                        ITMP2 = IJGTL
                        DBL_MB(MV+ITMP1-ILOV+LDHM*(ITMP2-JLOV)) = DBL_MB(MV+ITMP1-ILOV+LDHM*(ITMP2-JLOV))+SCL*AJCL(ICLSTA+ICL-1,IAJ)
                      end if
                    end do
                  end do
                end do
              end do
            end do
          end do
        end if
#     endif
      end if
    end do
  end do

# ifdef _MOLCAS_MPP_
  if ((iParRHS == 2) .and. (JLOV > 0)) call GA_Release_Update(lg_HM,ILOV,IHIV,JLOV,JHIV,MV,LDV)
# endif
  call RHS_SAVE(NASM,NISM,lg_HM,iCASE,iSYM,iVEC)
  call RHS_FREE(lg_HM)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end subroutine ADDRHSH

#ifdef _MOLCAS_MPP_
subroutine GADSUM_ADDRHS(buff,nbuff)

  use caspt2_global, only: MAXBUF

  integer(kind=iwp), intent(in) :: nbuff
  real(kind=wp), intent(inout) :: buff(nbuff)
  integer(kind=iwp) :: istart

  ! GADGOP wrapper: avoid the 2 GB limit of 32-bit MPI
  do istart=1,nbuff,MAXBUF
    call GADGOP(buff(istart),min(nbuff-istart+1,MAXBUF),'+')
  end do

end subroutine GADSUM_ADDRHS
#endif

end module ADDRHS
