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

subroutine DERE4(NLEV,iSym0,NASA,NASC,NCONF,BDERA,BDERC,Clag)

use Index_Functions, only: nTri_Elem
use BDerNEV, only: Gact, Gder
use caspt2_global, only: IDTCEX, iPrGlb, LUCIEX
use caspt2_module, only: JSTATE, MXCI, NACTEL, NSYM, STSYM
use sguga, only: CIS, EXS, SGS
use NEVPT2_E4, only: do_xvec, do_yvec, ixyzend, ixyzsta, NEVPT2_E4_derivative1, NEVPT2_E4_derivative2, NEVPT2_E4_derivative3, &
                     NEVPT2_E4_XYder1, NEVPT2_E4_XYder2, NEVPT2_E4_XYVEC, NEVPT2_E4_ZVEC, NXY_work, NXYVEC, nxyzdim, NZVEC
use PrintLevel, only: verbose
use Symmetry_Info, only: Mul
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate, mma_MaxDBLE
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: NLEV, iSym0, NASA, NASC, NCONF
real(kind=wp), intent(in) :: BDERA(NASA,NASA), BDERC(NASC,NASC)
real(kind=wp), intent(inout) :: CLag(nconf)
integer(kind=iwp) :: I, IALEV, IBLEV, IDCI, IDX, IP1, IP2, IP3, IP4, IPXYEND, IPXYSTA, ISAB, ISSG1, ISSG2, ISSG3, ISSG4, ISTU, &
                     ISVX, iSym, ISYZ, ITLEV, IULEV, IVLEV, IXLEV, ixy_local, IYLEV, IZLEV, J, JDX, MEMMAX, MEMMAX_SAFE, &
                     memory_per_xy, memory_per_xy2, memory_per_z, NCI, NLEV2, NSGM1, NSGM2, NSGM3, NSGM4, NTRI1
real(kind=wp) :: cpe, cptf0, cptf10, cput, tioe, tiotf0, tiotf10, wallt
integer(kind=iwp), allocatable :: IJ2IDX(:,:), IDX2IJ(:,:)
real(kind=wp), allocatable :: BUF2(:,:), BUFT(:), CI(:), XYcont(:,:,:), XYcontder(:,:,:), XYDER(:,:,:), XYtmp(:,:,:), &
                              XYVEC(:,:,:), ZDER(:,:), ZVEC(:,:)

! translation tables for levels i,j to and from pair indices idx: IJ2IDX, IDX2IJ
! result buffer, maximum size is the largest possible ip1 range,
! which is set to nbuf1 later, i.e. a maximum of nlev2 <= mxlev**2

! Note
! If RAS, computing E_tu|Psi> for t > u is wrong?
! sigma1 with t > u can be called, but it is used only for <Psi|Eut, then u < t
! only lowering operation is permitted

if (nlev == 0) return
if (NACTEL == 0) return

NCI = CIS%NCSF(STSYM)
! This should not happen, but...
if (NCI == 0) return

! Here, for regular CAS or RAS cases.

call mma_allocate(IJ2IDX,nLev,nLev,Label='IJ2IDX')
call mma_allocate(IDX2IJ,2,nLev**2,Label='IDX2IJ')
IJ2IDX(:,:) = -1
IDX2IJ(:,:) = -1

! Special pair index idx2ij allows true RAS cases to be handled:
nlev2 = nlev**2
ntri1 = nTri_Elem(nlev-1)
!ntri2 = nTri_Elem(nlev)
idx = 0
do i=1,nlev-1
  do j=i+1,nlev
    idx = idx+1
    ij2idx(i,j) = idx
    idx2ij(1,idx) = i
    idx2ij(2,idx) = j
    jdx = nlev2+1-idx
    ij2idx(j,i) = jdx
    idx2ij(1,jdx) = j
    idx2ij(2,jdx) = i
  end do
end do
do i=1,nlev
  idx = ntri1+i
  ij2idx(i,i) = idx
  idx2ij(1,idx) = i
  idx2ij(2,idx) = i
end do

! Dummy values necessary for fooling syntax checkers:
call mma_MaxDBLE(memmax)
! Use *almost* all remaining memory:
memmax_safe = int(real(memmax,kind=wp)*0.95_wp)

call mma_allocate(CI,NCONF,Label='CI')
if (NCONF > 1) then
  IDCI = IDTCEX(JSTATE)
  call DDAFILE(LUCIEX,2,CI,NCONF,IDCI)
else
  CI(1) = One
end if

call mma_allocate(BUFT,MXCI,LABEL='BUFT')
memmax_safe = memmax_safe-MXCI

! memory_min = NLEV2*2 ! allocated in XYder2

memory_per_z = mxci*2+nlev2
memory_per_xy = mxci*2+4*nlev2+2*nTri_Elem(nlev2)
memory_per_xy2 = mxci*4+4*nlev2+2*nTri_Elem(nlev2)
if (.false.) then ! maybe for RAS? idk
  NXY_work = 1
  NXYVEC = 0
  NZVEC = 1
else if (memmax_safe/(memory_per_z+memory_per_xy2) > nlev2) then
  !! sufficient memory strat: construct XVEC and YVEC simultaneously
  NXY_work = nlev
  NXYVEC = nlev2
  NZVEC = nlev2
  do_xvec = .true.
  do_yvec = .true.
else if (memmax_safe/(memory_per_z+memory_per_xy) > nlev2) then
  !! reasonable memory strat: construct XVEC and YVEC one by one
  NXY_work = nlev
  NXYVEC = nlev2
  NZVEC = nlev2
  do_xvec = .true.
  do_yvec = .false.
else if (memmax_safe/(memory_per_z+memory_per_xy) > nlev) then
  !! insufficient memory strat: batching
  NXY_work = max(1,min(nlev,memmax_safe/nlev/(memory_per_z+memory_per_xy)))
  !nxy_work = 3
  NXYVEC = NLEV*NXY_work
  NZVEC = NXYVEC
  do_xvec = .true.
  do_yvec = .false.
else
  write(u6,'(" Too little memory left for MKBNEVAC_E4.")')
  i = (memory_per_z+memory_per_xy)*nlev
  write(u6,'(" Requested memory  : ", i12," words")') i
  write(u6,'(" allocatable memory: ", i12," words")') memmax_safe
  i = i-memmax_safe
  write(u6,'(" Need an additional memory of ",i12," words")') i
  i = int(real(i*RtoB,kind=wp)*1.0e-6_wp,kind=iwp)+1
  write(u6,'(" Please increase MOLCAS_MEM by roughly ",i6," (MB)")') i
  call ABEND()
end if

call mma_allocate(ZVEC,MXCI,NZVEC,LABEL='ZVEC')
call mma_allocate(Zder,MXCI,NZVEC,LABEL='Zder')

if (do_xvec .and. do_yvec) then
  call mma_allocate(XYVEC,MXCI,NXYVEC,2,LABEL='XYVEC')
  call mma_allocate(XYder,MXCI,NXYVEC,2,LABEL='XYder')
else if (do_xvec .and. (.not. do_yvec)) then
  call mma_allocate(XYVEC,MXCI,NXYVEC,1,LABEL='XYVEC')
  call mma_allocate(XYder,MXCI,NXYVEC,1,LABEL='XYder')
end if

call mma_allocate(XYtmp,NXYVEC,NLEV2,2,LABEL='XYtmp')
call mma_allocate(XYcont,nxyvec,2,nlev2,LABEL='XYcont')
call mma_allocate(XYcontder,nxyvec,2,nTri_Elem(nlev2),LABEL='XYcontder')

if (iPrGlb >= VERBOSE) then
  write(u6,*)
  write(u6,'(2X,A)') 'Constructing E4 terms for NEVPT2'
  write(u6,'(2X,A,1X,I1)') 'Symmetry of B matrix:',iSym0
  write(u6,'(2X,A,F16.9,A)') ' memory available:   ',real(memmax*RtoB,kind=wp)*1.0e-9_wp,' GB'
  i = nconf+mxci+memory_per_xy2*nlev2+memory_per_z*nlev2
  write(u6,'(2X,A,F16.9,A)') ' memory sufficient:  ',real(i*RtoB,kind=wp)*1.0e-9_wp,' GB'
  i = nconf+mxci+memory_per_xy*nlev2+memory_per_z*nlev2
  write(u6,'(2X,A,F16.9,A)') ' memory reasonable:  ',real(i*RtoB,kind=wp)*1.0e-9_wp,' GB'
  i = nconf+mxci+memory_per_xy*nlev+memory_per_z*nlev
  write(u6,'(2X,A,F16.9,A)') ' memory insufficient:',real(i*RtoB,kind=wp)*1.0e-9_wp,' GB'
  write(u6,*)
  if (NXYVEC == 0) then
    write(u6,'(" (0) special case")')
  else if (do_xvec .and. do_yvec) then
    write(u6,'(" (1) memory sufficient strat")')
  else if (do_xvec) then
    if ((nxyvec == nlev2) .and. (nzvec == nlev2)) then
      write(u6,'(" (2) memory reasonable strat")')
    else
      write(u6,'(" (3) memory insufficient strat")')
    end if
  end if
  if (do_xvec .and. do_yvec) then
    i = nconf+mxci+memory_per_xy*nxyvec+memory_per_z*nzvec
  else if (do_xvec .and. (.not. do_yvec)) then
    i = nconf+mxci+memory_per_xy2*nxyvec+memory_per_z*nzvec
  end if
  write(u6,'(2X,A,F16.9,A)') ' memory allocated:   ',real(i*RtoB,kind=wp)*1.0e-9_wp,' GB'
  call mma_MaxDBLE(memmax)
  write(u6,'(2X,A,F16.9,A)') ' memory left     :   ',real(memmax*RtoB,kind=wp)*1.0e-9_wp,' GB'
  write(u6,'(2X," allocated/total Zvecs. : ",i3,"/",i3)') nzvec,nlev2
  write(u6,'(2X," allocated/total XYvecs.: ",i3,"/",i3)') nxyvec,nlev2
  call xflush(u6)
end if

call GASync()

if (NXYVEC > 0) then
  SymmetryLoop: do iSym=1,nSym !! Symmetry of CSF (I) for RI
    !! Stride loop of Z_{t,u}^I
    !! Construct Z_{t,u}^I for ixyzsta <= t <= ixyzend and 1 <= u <= NLEV
    XYvecLoop2: do ixy_local=1,NLEV,NXY_work
      ipxysta = 1
      ipxyend = NLEV2
      ixyzsta = ixy_local
      ixyzend = min(NLEV,ixyzsta-1+NXY_work)
      nxyzdim = ixyzend-ixyzsta+1
      if (iPrGlb >= VERBOSE) call TIMING(CPTF0,CPE,TIOTF0,TIOE)
      !! Construct ZVEC
      call NEVPT2_E4_ZVEC(NLEV,idx2ij,Gact,CI,ZVEC,XYVEC)
      !! Construct XVEC and YVEC
      call NEVPT2_E4_XYVEC(iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,ZVEC,XYVEC)
      if (iPrGlb >= VERBOSE) then
        call TIMING(CPTF10,CPE,TIOTF10,TIOE)
        CPUT = CPTF10-CPTF0
        WALLT = TIOTF10-TIOTF0
        write(u6,'(a,2f10.3)') ' AC_E4(1): CPU/WALL TIME=',cput,wallt
        call TIMING(CPTF0,CPE,TIOTF0,TIOE)
      end if
      !! Construct XYcont for X*E and Y*E
      call NEVPT2_E4_derivative1(iSym0,iSym,NLEV,idx2ij,ipxysta,ipxyend,BDERA,BDERC,XYcont)
      if (iPrGlb >= VERBOSE) then
        call TIMING(CPTF10,CPE,TIOTF10,TIOE)
        CPUT = CPTF10-CPTF0
        WALLT = TIOTF10-TIOTF0
        write(u6,'(a,2f10.3)') ' AC_E4(2): CPU/WALL TIME=',cput,wallt
        call TIMING(CPTF0,CPE,TIOTF0,TIOE)
      end if
      !! Construct XYcont and XYder for X*E*E and Y*E*E
      call NEVPT2_E4_derivative2(iSym0,iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,CI,BDERA,BDERC,XYcont,XYcontder,XYder,ZVEC, &
                                 XYtmp)

      if (iPrGlb >= VERBOSE) then
        call TIMING(CPTF10,CPE,TIOTF10,TIOE)
        CPUT = CPTF10-CPTF0
        WALLT = TIOTF10-TIOTF0
        write(u6,'(a,2f10.3)') ' AC_E4(3): CPU/WALL TIME=',cput,wallt
        call TIMING(CPTF0,CPE,TIOTF0,TIOE)
      end if
      !! Construct XYder for X*E and Y*E
      call NEVPT2_E4_derivative3(iSym,NLEV,idx2ij,ipxysta,ipxyend,BUFT,CI,XYcont,XYder)
      if (iPrGlb >= VERBOSE) then
        call TIMING(CPTF10,CPE,TIOTF10,TIOE)
        CPUT = CPTF10-CPTF0
        WALLT = TIOTF10-TIOTF0
        write(u6,'(a,2f10.3)') ' AC_E4(4): CPU/WALL TIME=',cput,wallt
        call TIMING(CPTF0,CPE,TIOTF0,TIOE)
      end if
      !! Construct Zder from XYder
      call NEVPT2_E4_XYder1(iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,CI,XYder,Gder,Zder)
      if (iPrGlb >= VERBOSE) then
        call TIMING(CPTF10,CPE,TIOTF10,TIOE)
        CPUT = CPTF10-CPTF0
        WALLT = TIOTF10-TIOTF0
        write(u6,'(a,2f10.3)') ' AC_E4(5): CPU/WALL TIME=',cput,wallt
        call TIMING(CPTF0,CPE,TIOTF0,TIOE)
      end if
      !! Complete CI and ERI derivatives of Zder
      call NEVPT2_E4_XYder2(iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,CI,Gact,XYvec,XYder,XYcont,XYcontder,ZVEC,Zder,CLag,Gder, &
                            XYtmp)
      if (iPrGlb >= VERBOSE) then
        call TIMING(CPTF10,CPE,TIOTF10,TIOE)
        CPUT = CPTF10-CPTF0
        WALLT = TIOTF10-TIOTF0
        write(u6,'(a,2f10.3)') ' AC_E4(6): CPU/WALL TIME=',cput,wallt
        call TIMING(CPTF0,CPE,TIOTF0,TIOE)
      end if
      if (do_xvec .and. (.not. do_yvec)) then
        call TIMING(CPTF0,CPE,TIOTF0,TIOE)
        do_xvec = .false.
        do_yvec = .true.

        call NEVPT2_E4_ZVEC(NLEV,idx2ij,Gact,CI,ZVEC,XYVEC)
        call NEVPT2_E4_XYVEC(iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,ZVEC,XYVEC)
        call NEVPT2_E4_derivative1(iSym0,iSym,NLEV,idx2ij,ipxysta,ipxyend,BDERA,BDERC,XYcont)
        call NEVPT2_E4_derivative2(iSym0,iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,CI,BDERA,BDERC,XYcont,XYcontder,XYder,ZVEC, &
                                   XYtmp)
        call NEVPT2_E4_derivative3(iSym,NLEV,idx2ij,ipxysta,ipxyend,BUFT,CI,XYcont,XYder)
        call NEVPT2_E4_XYder1(iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,CI,XYder,Gder,Zder)
        call NEVPT2_E4_XYder2(iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,CI,Gact,XYvec,XYder,XYcont,XYcontder,ZVEC,Zder,CLag, &
                              Gder,XYtmp)

        if (iPrGlb >= VERBOSE) then
          call TIMING(CPTF10,CPE,TIOTF10,TIOE)
          CPUT = CPTF10-CPTF0
          WALLT = TIOTF10-TIOTF0
          write(u6,'(a,2f10.3)') ' AC_E4(7): CPU/WALL TIME=',cput,wallt
        end if
        do_xvec = .true.
        do_yvec = .false.
      end if
    end do XYvecLoop2
  end do SymmetryLoop
else
  write(u6,*) 'notyet'
  call abend()
  !! least memory strategy: directly compute E*E*E*E
  call mma_deallocate(BUFT)
  call mma_allocate(BUFT,1,LABEL='BUFT')

  call mma_MaxDBLE(memmax)
  memmax_safe = int(real(memmax,kind=wp)*0.95_wp)
  if (memmax_safe/mxci <= 4) then
    write(u6,*) ' Too little memory left for MKBNEVAC_E4.'
    write(u6,*) ' Need at least 4 vector of length MXCI=',MXCI
    call ABEND()
  end if
  call mma_allocate(BUF2,MXCI,4,LABEL='BUF2')

  do ip4=1,nlev2
    ialev = idx2ij(1,ip4)
    iblev = idx2ij(2,ip4)
    isab = Mul(SGS%ism(ialev),SGS%ism(iblev))
    issg4 = Mul(isab,stsym)
    nsgm4 = CIS%ncsf(issg4)
    !ia = L2ACT(ialev)
    !ib = L2ACT(iblev)
    BUF2(1:nsgm4,4) = Zero
    call SG_Epq_Psi(SGS,CIS,EXS,IALEV,IBLEV,One,STSYM,CI,BUF2(:,4))
    do ip3=ip4,nlev2
      iylev = idx2ij(1,ip3)
      izlev = idx2ij(2,ip3)
      isyz = Mul(SGS%ism(iylev),SGS%ism(izlev))
      issg3 = Mul(isyz,issg4)
      nsgm3 = CIS%ncsf(issg3)
      !iy = L2ACT(iylev)
      !iz = L2ACT(izlev)
      BUF2(1:nsgm3,3) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,IYLEV,IZLEV,One,issg4,BUF2(:,4),BUF2(:,3))
      write(u6,'(a,2f10.3)') ' AC_E4(1): CPU/WALL TIME=',cput,wallt
    end do
    do ip2=ip3,nlev2
      ivlev = idx2ij(1,ip2)
      ixlev = idx2ij(2,ip2)
      isvx = Mul(SGS%ism(ivlev),SGS%ism(ixlev))
      issg2 = Mul(isvx,issg3)
      nsgm2 = CIS%ncsf(issg2)
      !iv = L2ACT(ivlev)
      !ix = L2ACT(ixlev)
      BUF2(1:nsgm2,2) = Zero
      call SG_Epq_Psi(SGS,CIS,EXS,IVLEV,IXLEV,One,issg3,BUF2(:,3),BUF2(:,2))
      do ip1=ip2,nlev2
        itlev = idx2ij(1,ip1)
        iulev = idx2ij(2,ip1)
        istu = Mul(SGS%ism(itlev),SGS%ism(iulev))
        issg1 = Mul(istu,issg2)
        nsgm1 = CIS%ncsf(issg1)
        !it = L2ACT(itlev)
        !iu = L2ACT(iulev)
        BUF2(1:nsgm1,1) = Zero
        call SG_Epq_Psi(SGS,CIS,EXS,ITLEV,IULEV,One,issg2,BUF2(:,2),BUF2(:,1))
      end do
    end do
  end do
end if

call mma_deallocate(BUFT)

call mma_deallocate(XYVEC,safe='*')
call mma_deallocate(XYtmp,safe='*')
call mma_deallocate(XYcont,safe='*')
call mma_deallocate(XYcontder,safe='*')

call mma_deallocate(XYDER,safe='*')
call mma_deallocate(ZVEC,safe='*')
call mma_deallocate(Zder,safe='*')

call mma_deallocate(CI)

call mma_deallocate(ij2idx)
call mma_deallocate(idx2ij)

return

end subroutine DERE4
