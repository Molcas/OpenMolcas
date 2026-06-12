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

subroutine MKBNEV()

use Index_Functions, only: nTri_Elem
use caspt2_global, only: iPrGlb, LUSBT, LUSOLV
use caspt2_module, only: NASHT_ => NASHT, NASUP, NG1, NG2, NG3, NINDEP, NSYM
  use PrintLevel, only: debug, verbose
  use stdalloc, only: mma_allocate, mma_deallocate
  use caspt2_global, only: LUSBT, LUSOLV, SGS
  use EQSOLV, only: IDBMAT, IDSMAT
  use SC_NEVPT2, only: Do_SC, ECORR_SC, IDBMAT_NEVPT2, OVLAPS_SC
  use NEVPT2_mod, only: NASHT
#ifdef _MOLCAS_MPP_
  use NEVPT2_E4, only: MAXBUF
  use caspt2_module, only: MXCI
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6, byte

  implicit none
integer(kind=iwp) :: ICASE, IDISK, iLUID, ISYM, NLEV, NTRI
real(kind=wp) :: cpe, cptf0, cptf10, cput, DUM(1), tioe, tiotf0, tiotf10, wallt
  integer(kind=byte), allocatable :: idxG3(:,:)
real(kind=wp), allocatable :: G1(:,:), G2(:,:,:,:), G3(:), Gact(:,:,:,:), Hact(:,:), Hbar(:,:), Htilde(:,:), WRK(:)

  if (IPRGLB >= VERBOSE) then
    write(u6,*)
    write(u6,*)' Construct B matrices'
  end if

! It seems G2 is not a standard normal-ordered RDM.
! G2(standard)(a,b,c,d) = G2(MOLCAS)(a,c,b,d)
! That is, in MOLCAS, G2(a,b,c,d) = <0|a+ c+ d b|0> = D(ac,bd)
! G3(standard)(tuvxyz) = G3(MOLCAS)(txu,yvz)
! G3(t,u,v,x,y,z) = <0|t+ v+ y+ z x u|0> -> G3(standard)(tvyuxz)

  nAshT = nAshT_ ! define the number of the active orbitals used in NEVPT2_MOD
#ifdef _MOLCAS_MPP_
  MAXBUF = 2000000000/(MXCI*8)
#endif

  if (nAshT /= 0) then
    call mma_allocate(Hact,nAshT,nAshT,Label='Hact')
    call mma_allocate(Gact,nAshT,nAshT,nAshT,nAshT,Label='Gact')
    call mma_allocate(Hbar,nAshT,nAshT,Label='Hbar')
    call mma_allocate(Htilde,nAshT,nAshT,Label='Htilde')
    ! Construct the active part of the integrals
    call nevint(nAshT,Hact,Gact,Hbar,Htilde)

!   check_eigen = .true.
!   if (check_eigen) call check_eigenfunction(nAshT,Hact,Gact)

    call mma_allocate(G1,nAshT,nAshT,Label='G1')
    call mma_allocate(G2,nAshT,nAshT,nAshT,nAshT,Label='G2')
    call mma_allocate(G3,NG3,Label='G3')
    call mma_allocate(idxG3,6,NG3,label='idxG3')
    call PT2_GET(NG1,' GAMMA1',G1)
    call PT2_GET(NG2,' GAMMA2',G2)
    call PT2_GET(NG3,' GAMMA3',G3)
    iLUID=0
    call I1DAFILE(LUSOLV,2,idxG3,6*NG3,iLUID)

    if (IPRGLB >= DEBUG) then
      write(u6,'("DEBUG> ",A)') 'case SYM B-MATRIX NORM'
      write(u6,'("DEBUG> ",A)') '==== === ============='
    end if

    if (IPRGLB >= verbose) call TIMING(CPTF0,CPE,TIOTF0,TIOE)
    call MKBNEVG(nAshT,Hact,Gact,G1,G2)
    call MKBNEVE(nAshT,Hact,Gact,G1,G2)

    call MKBNEVF(nAshT,NG3,Hbar,Gact,G2,G3,idxG3)
    call MKBNEVB(nAshT,NG3,Hbar,Gact,G1,G2,G3,idxG3)

    call MKBNEVD(nAshT,NG3,Hbar,Gact,G1,G2,G3,idxG3)
    call GASync()

    if (IPRGLB >= verbose) then
      call TIMING(CPTF10,CPE,TIOTF10,TIOE)
      CPUT =CPTF10-CPTF0
      WALLT=TIOTF10-TIOTF0
    write(u6,'(a,2f10.3)') ' MKBNEV1 : CPU/WALL TIME=',cput,wallt
      call TIMING(CPTF0,CPE,TIOTF0,TIOE)
    end if

    ! I did not come up with a good parallel strategy for the A and C subspaces

    call MKBNEVAC_E3(nAshT,NG3,Hact,Htilde,Gact,G1,G2,G3,idxG3)
    call GASync()

    if (IPRGLB >= verbose) then
      call TIMING(CPTF10,CPE,TIOTF10,TIOE)
      CPUT =CPTF10-CPTF0
      WALLT=TIOTF10-TIOTF0
    write(u6,'(a,2f10.3)') ' MKBNEV2 : CPU/WALL TIME=',cput,wallt
      call TIMING(CPTF0,CPE,TIOTF0,TIOE)
    end if

    call mma_deallocate(Hact)
    call mma_deallocate(Hbar)
    call mma_deallocate(Htilde)
    call mma_deallocate(G1)
    call mma_deallocate(G2)
    call mma_deallocate(G3)
    call mma_deallocate(idxG3)

    nLev = SGS%nLev
    call MKBNEVAC_E4(nAshT,NLEV,Gact)

    if (IPRGLB >= verbose) then
      call TIMING(CPTF10,CPE,TIOTF10,TIOE)
      CPUT =CPTF10-CPTF0
      WALLT=TIOTF10-TIOTF0
    write(u6,'(a,2f10.3)') ' MKBNEV3 : CPU/WALL TIME=',cput,wallt
    end if

    call mma_deallocate(Gact)
  end if

! For completeness, even case H has formally S and B
! matrices. This costs nothing, and saves conditional
! looping, etc in the rest  of the routines.
  DUM(1)=Zero
  do ISYM=1,NSYM
    do ICASE=12,13
      if (NINDEP(ISYM,ICASE) > 0) then
        IDISK=IDSMAT(ISYM,ICASE)
        call DDAFILE(LUSBT,1,DUM,1,IDISK)
        IDISK=IDBMAT(ISYM,ICASE)
        call DDAFILE(LUSBT,1,DUM,1,IDISK)
      end if
    end do
  end do

  ! Copy the original B matrix
  if (Do_SC) then
  ECORR_SC(:,:) = Zero
  OVLAPS_SC(:,:) = Zero
    do iSym = 1, nSym
      do iCase = 1, 11
      if ((iCase == 1) .or. (iCase == 4)) cycle ! done in MKBNEVA_E3 and E4
        if (NINDEP(iSym,iCase) > 0) then
        NTRI = nTri_Elem(NASUP(iSym,iCase))
        call mma_allocate(WRK,NTRI,Label='WRK')
          idisk = IDBMAT(iSym,iCase)
        call DDAFILE(LUSBT,2,WRK,NTRI,IDISK)
          idisk = IDBMAT_NEVPT2(iSym,iCase,1)
        call DDAFILE(LUSBT,1,WRK,NTRI,IDISK)
          call mma_deallocate(WRK)
        end if
      end do
    end do
  end if

end subroutine MKBNEV
