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
! Copyright (C) 2026, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine nevint(nAshT,INT1,INT2,Hbar,Htilde)

  use CHOVEC_IO, only: NVLOC_CHOBATCH
! use ChoCASPT2
  use caspt2_global, only: FIMO
  use caspt2_module, only: NAES, NASH, NBSQT, NBTCH, NBTCHES, NISH, NORB, NSSH, NORB, NSYM, IfChol
  use stdalloc, only: mma_allocate, mma_deallocate
  use Definitions, only: iwp,wp,u6
  use Constants, only: Zero, One, Half
  use Symmetry_Info, only: Mul
#ifdef _MOLCAS_MPP_
  use Para_Info, only: Is_Real_Par
#endif

  implicit none

  integer(kind=iwp), intent(in) :: nAshT
  real(kind=wp), intent(out) :: INT1(nAshT,nAshT), INT2(nAshT,nAshT,nAshT,nAshT), Hbar(nAshT,nAshT), Htilde(nAshT,nAshT)

  integer(kind=iwp), allocatable :: BGRP(:,:)
  real(kind=wp), allocatable :: PIQK(:), BRA(:), KET(:), WRK1(:), WRK2(:)

  integer(kind=iwp), parameter :: Active=2, Inactive=1, Virtual=3
  integer(kind=iwp) :: nSh(8,3)

  integer(kind=iwp) :: iAshI, IB, IBEND, IBGRP, IBSTA, IBUF, IB1, IB2, iSym, iSymA, iSymB, iSymI, iSymJ, &
                       iSymT, iSymU, iSymV, iSymVX, iSymX, IT, ITABS, ITTOT, IU, IUABS, IUTOT, IV, IVABS, IVTOT, IX, IXABS, IXTOT, &
                       jAshI, jSym, kAshI, MXBGRP, MXPIQK, NADDBUF, NBGRP, NCHOBUF, &
                       nKet, NFIMOES, NFNXT, nIshI, nOrbI, NV
  integer(kind=iwp) :: LBRASM, ISYI, NI, iOffi, ISYP, NP, iOffp, NPI, NBRASM, &
                       LKETSM, ISYK, NK, iOffK, ISYQ, NQ, iOffQ, NQK, NKETSM, &
                       NTUVX, NUMERR

  real(kind=wp) :: Scal, tmp

  INT1(:,:) = Zero
  INT2(:,:,:,:) = Zero
  Hbar(:,:) = Zero
  Htilde(:,:) = Zero

  NFNXT=0
  do iSym = 1, nSym
    NFIMOES = NFNXT
    NFNXT = NFNXT + (nOrb(iSym)*(nOrb(iSym)+1))/2
    if (nAsh(iSym) == 0) cycle
    nIshI = nIsh(iSym)
    nOrbI = nOrb(iSym)
    do IU = 1, nAsh(iSym)
      IUTOT = nIsh(iSym) + IU
      IUABS = IU + nAes(iSym)
      do IT = IU, nAsh(iSym)
        ITTOT = nIsh(iSym) + IT
        ITABS = IT + nAes(iSym)
        IBUF = NFIMOES+(ITTOT*(ITTOT-1))/2+IUTOT
        tmp = FIMO(IBUF)
        INT1(ITABS,IUABS) = tmp
        INT1(IUABS,ITABS) = tmp
      end do
    end do
  end do

  Hbar(:,:) = INT1(:,:)
  Htilde(:,:) = INT1(:,:)

  call mma_allocate(WRK1,NBSQT,Label='WRK1')
  call mma_allocate(WRK2,NBSQT,Label='WRK2')
!
! --- Two-Electron Integral
!
  iSymA = 1
  iSymI = 1
  iSymB = 1
  iSymJ = 1
  if (IfChol) then
    NTUVX = NASHT**4
    nSh(1:NSYM,Inactive) = NISH(1:NSYM)
    nSh(1:NSYM,Active  ) = NASH(1:NSYM)
    nSh(1:NSYM,Virtual ) = NSSH(1:NSYM)
    do JSYM=1,NSYM
      IB1=NBTCHES(JSYM)+1
      IB2=NBTCHES(JSYM)+NBTCH(JSYM)

      MXBGRP=IB2-IB1+1
      if (MXBGRP <= 0) cycle
      call mma_allocate(BGRP,2,MXBGRP,Label='BGRP')
      IBGRP=1
      do IB=IB1,IB2
       BGRP(1,IBGRP) = IB
       BGRP(2,IBGRP) = IB
       IBGRP=IBGRP+1
      end do
      NBGRP=MXBGRP

      call MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,NCHOBUF,MXPIQK,NADDBUF)
      call mma_allocate(PIQK,MXPIQK,Label='PIQK')
      call mma_allocate(BRA,NCHOBUF,Label='BRA')
      call mma_allocate(KET,NCHOBUF,Label='KET')
      do IBGRP=1,NBGRP

        IBSTA=BGRP(1,IBGRP)
        IBEND=BGRP(2,IBGRP)

        NV=0
        do IB=IBSTA,IBEND
          NV=NV+NVLOC_CHOBATCH(IB)
        end do

        !! int2(tuvx) = (tu|vx)
        !! This can be computed without frozen orbitals
        call Get_Cholesky_Vectors(Active,Active,JSYM,KET,SIZE(KET),nKet,IBSTA,IBEND)

        LBRASM=1
        do ISYI=1,NSYM
          NI=NASH(ISYI)
          iOffi=NAES(iSYI)
          if (NI == 0) cycle
          ISYP=Mul(ISYI,JSYM)
          NP=NASH(ISYP)
          iOffp=NAES(iSYP)
          if (NP == 0) cycle
          NPI=NP*NI
          NBRASM=NPI*NV
          LKETSM=1
          do ISYK=1,NSYM
            NK=NASH(ISYK)
            iOffK=NAES(iSYK)
            if (NK == 0) cycle
            ISYQ=Mul(ISYK,JSYM)
            NQ=NASH(ISYQ)
            iOffQ=NAES(iSYQ)
            if (NQ == 0) cycle
            NQK=NQ*NK
            NKETSM=NQK*NV

            if (NPI*NQK > mxPIQK) then
              write(u6,*) 'NPIQK larger than mxPIQK in TUVX, bug?'
              call AbEnd()
            end if
            call DGEMM_('N','T',NPI,NQK,NV,One,KET(LBRASM),NPI,KET(LKETSM),NQK,Zero,PIQK,NPI)
            call ADDTUVX(NP,NI,NQ,NK,NASHT,iOffP,iOffI,iOffQ,iOffK,INT2,nTUVX,PIQK,NPI*NQK,NUMERR)
            LKETSM=LKETSM+NKETSM
          end do
          LBRASM=LBRASM+NBRASM
        end do

!       If (IBGRP.EQ.1) SCAL = Zero
!       If (IBGRP.NE.1) SCAL = One
!       Call DGEMM_('N','T',NASH(JSYM)**2,NASH(JSYM)**2,NV, &
!                   One ,KET,NASH(JSYM)**2,KET,NASH(JSYM)**2, &
!                   SCAL,INT2,NASH(JSYM)**2)
      end do
      call mma_deallocate(PIQK)
      call mma_deallocate(BRA)
      call mma_deallocate(KET)
      call mma_deallocate(BGRP)
    end do
#ifdef _MOLCAS_MPP_
    if (is_real_par()) call GADGOP(INT2,nAshT**4,'+')
#endif
  else
    do iSymX = 1, nSym
      do iSymV = 1, nSym
        iSymVX = Mul(iSymV,iSymX)
        do iSymU = 1, nSym
          iSymT = Mul(iSymU,iSymVX)
          do IX = 1, nAsh(iSymX)
            IXTOT = nIsh(iSymX) + IX
            IXABS = IX + nAes(iSymX)
            do IV = 1, nAsh(iSymV)
              IVTOT = nIsh(iSymV) + IV
              IVABS = IV + nAes(iSymV)
              call Coul(iSymT,iSymU,iSymV,iSymX,IVTOT,IXTOT,WRK1,WRK2)
              do IU = 1, nAsh(iSymU)
                IUTOT = nIsh(iSymU) + IU
                IUABS = IU + nAes(iSymU)
                do IT = 1, nAsh(iSymT)
                  ITTOT = nIsh(iSymT) + IT
                  ITABS = IT + nAes(iSymT)
                  IBUF = ITTOT + nOrb(iSymT)*(IUTOT-1)
                  INT2(ITABS,IUABS,IVABS,IXABS) = WRK1(IBUF)
!                 write (*,'(4i3,f20.10)') itabs,iuabs,ivabs,ixabs,wrk1(ibuf)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end if
  call mma_deallocate(WRK1)
  call mma_deallocate(WRK2)

!do isymt = 1, nasht
!do isymu = 1, nasht
!do isymv = 1, nasht
!do isymx = 1, nasht
!write (*,'(4i3,f20.10)') isymt,isymu,isymv,isymx,int2(isymt,isymu,isymv,isymx)
!end do
!end do
!end do
!end do

  ! Construct Hbar^eff (after Eq. (A13)) and Htilde

  do iAshI = 1, NASHT ! nAsh(iSym)
    do jAshI = 1, NASHT ! nAsh(iSym)
      Scal = Zero
      do kAshI = 1, NASHT ! nAsh(iSym)
        Scal = Scal + INT2(kAshI,jAshI,kAshI,iAshI)
      end do
      Hbar(iAshI,jAshI) = Hbar(iAshI,jAshI) - Half*Scal
      Htilde(iAshI,jAshI) = Htilde(iAshI,jAshI) - Scal
    end do
  end do

  return

end subroutine nevint
