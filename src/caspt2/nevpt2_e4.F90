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

module NEVPT2_E4

! Variables and subroutines relevant to the E*E*E*E contributions
!
! 1. NEVPT2_E4_ZVEC: Construct ZVEC
!    Z_{t,u}^I = \sum_{v,x} (tu|vx) <Phi_I|E_{v,x}|Psi>
! 2. NEVPT2_E4_XYVEC: Construct XVEC and YVEC
!    X_{v,u}^I = <Phi_I|E_{x,u}|Phi_J> Z_{v,x}^J
!    Y_{v,t}^I = <Phi_I|E_{t,x}|Phi_J> Z_{v,x}^J
! 3. NEVPT2_E4_contract1: Contract XVEC and YVEC with <Psi|E_{tu}|I>
!    A_{t,u,y,z} = \sum_I <Psi|E_{tu}|I>*X_{y,z}^I
!    B_{t,u,y,z} = \sum_I <Psi|E_{tu}|I>*Y_{y,z}^I
!    (used for some E*E*E contributions and \tilde{E} terms)
! 4. NEVPT2_E4_contract2: Contract XVEC and YVEC with <I|E_{vx}E_{yz}|Psi>
!    C_{t,u,v,x,y,z} = \sum_I <I|E_{vx}E_{yz}|Psi>*X_{t,u}^I
!    D_{t,u,v,x,y,z} = \sum_I <I|E_{vx}E_{yz}|Psi>*Y_{t,u}^I
!    (why bra state? I defined v,x,y,z of A in a different order?)

! Part 1
! First, construct ZVEC (Z_{t,u}^I)
! - cost: nlev2*sigma1 + nlev2*nzvec*daxpy
! Next, construct XVEC and YVEC
! - cost: nlev*(nlev-1)*sigma1 + 2*nlev*(nlev-1)*nxyvec*daxpy
!
! Part 2
! <Psi|E_{t,u}|I>*X_{v,x}^I and <Psi|E_{t,u}|I>*Y_{v,x}^I
! - cost: nlev2*sigma1 + 2*nlev2*dgemv(nxyvec)
!
! Part 3
! <Psi|E_{t,u}E_{v,x}|I>*X_{y,z}^I and <Psi|E_{t,u}E_{v,x}|I>*Y_{y,z}^I
! - cost: nlev2*(nlev2-1)*sigma1 + 2*nlev2*dgemv(nxyvec,nlev2-1) /2

! Part 3 is clearly the most computationally expensive part, so increase NXYVEC as much as possible.
! On the other hand, we cannot ignore the importance of ZVEC construction.
! NBZ*(nlev2*sigma1 + nlev2*nzvec*daxpy) + NBXY*(nlev*(nlev-1)*sigma1 + 2*nlev*(nlev-1)*nxyvec*daxpy)
! \approx NBZ*nlev2*(sigma1 + nzvec*daxpy) NBXY*nlev2*(sigma1 + 2*nxyvec*daxpy)
! =       nlev2 * [ NBZ*(sigma1 + nzvec*daxpy) + NBXY*(sigma1 + 2*Nxyvec*daxpy) ]
! =       nlev2 * [ (NBZ+NBXY)*sigma1 + (NBZ*nzvec + NBXY*NXYVEC)*daxpy ]
! \approx nlev2 * [ (NBZ+NBXY)*sigma1 + 3*nlev2*daxpy ]
! So, sum of the Z and XY batch should be minimized -> NBZ = NBXY

! define three strategies depending on the available memory
! 1. > 3*nlev2*MXCI --> ZVEC -> XVEC + YVEC -> E(*E) contraction
! 2. > 2*nlev2*MXCI --> ZVEC -> XVEC -> E(*E) contraction -> ZVEC -> YVEC -> E(*E) contraction
! 3. < 2*nlev2*MXCI --> batched 2

! Consider distributed memory strategy
! ... at the moment, no

use sguga, only: sg_epq_psi
use Index_Functions, only: iTri, nTri_Elem
use general_data, only: STSYM
use caspt2_module, only: MXCI, NTUVES
use sguga_states, only: CIS, EXS, SGS
use SUPERINDEX, only: KTUV
use Symmetry_Info, only: Mul
use Task_Manager, only: Free_Tsk, Init_Tsk, Rsv_Tsk
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: ixyzend = 0, ixyzsta = 0, NXY_work = 0, NXYVEC = 0, nxyzdim = 0, NZVEC = 0
#ifdef _MOLCAS_MPP_
integer(kind=iwp) :: MAXBUF = 0
#endif
logical(kind=iwp) :: do_xvec = .false., do_yvec = .false.
logical(kind=iwp), parameter :: do_zder = .true.
integer(kind=iwp), parameter :: istate=1

public :: do_xvec, do_yvec, ixyzend, ixyzsta, NEVPT2_E4_contract1, NEVPT2_E4_contract2, NEVPT2_E4_derivative1, &
          NEVPT2_E4_derivative2, NEVPT2_E4_derivative3, NEVPT2_E4_XYder1, NEVPT2_E4_XYder2, NEVPT2_E4_XYVEC, NEVPT2_E4_ZVEC, &
          NXY_work, NXYVEC, nxyzdim, NZVEC
#ifdef _MOLCAS_MPP_
public :: MAXBUF
#endif

contains

!-----------------------------------------------------------------------

subroutine NEVPT2_E4_ZVEC(NLEV,idx2ij,Gact,CI,ZVEC,WRK)

  integer(kind=iwp), intent(in) :: NLEV, idx2ij(2,NLEV**2)
  real(kind=wp), intent(in) :: Gact(:,:,:,:), CI(:)
  real(kind=wp), intent(inout) :: ZVEC(MXCI,NZVEC)
  real(kind=wp), intent(out) :: WRK(MXCI,NZVEC) !! This is XYVEC outside this subroutine
  integer(kind=iwp) :: ibuf, ID, ip1, ip2, it, itlev, iu, iulev, iv, ivlev, ix, ixlev, nlev2, nTask, nxy
  real(kind=wp) :: scal
  real(kind=wp), allocatable :: Gact_sort(:,:)

  ! Construct ZVEC
  ! Z_{t,u}^I = \sum_{v,x} (tu|vx) <Phi_I|E_{v,x}|Psi>

  if (NXY_work == NLEV) then
    nxy = NLEV*NLEV
  else
    nxy = NLEV*nxyzdim
  end if
  NLEV2 = NLEV**2
  call mma_allocate(Gact_sort,nxy,NZVEC,LABEL='Gact_sort')

  nTask = NLEV2
  call Init_Tsk(ID,nTask)

  ibuf = 0
  Scal = Zero
  do while (Rsv_Tsk(ID,ip2))
    ivlev = idx2ij(1,ip2)
    ixlev = idx2ij(2,ip2)
!   isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
    iv=SGS(istate)%L2ACT(ivlev)
    ix=SGS(istate)%L2ACT(ixlev)
    ibuf = ibuf+1
    if (NXY_work == NLEV) then
      do ip1=1,nlev2
        itlev = idx2ij(1,ip1)
        iulev = idx2ij(2,ip1)
        it=SGS(istate)%L2ACT(itlev)
        iu=SGS(istate)%L2ACT(iulev)
        Gact_sort(ip1,ibuf) = Gact(it,iu,iv,ix)
      end do
    else
      do itlev=ixyzsta,ixyzend
        do iulev=1,nlev
          !ip1 = ij2idx(ixy_local,iulev_local)
          !itlev = idx2ij(1,ip1)
          !iulev = idx2ij(2,ip1)
!         it=SGS(istate)%L2ACT(itlev)
!         iu=SGS(istate)%L2ACT(iulev)
          !Gact_sort(iulev_local) = Gact(it,iu,iv,ix)
          Gact_sort(itlev-ixyzsta+1+nxyzdim*(iulev-1),ibuf) = Gact(itlev,iulev,iv,ix)
        end do
      end do
    end if
    WRK(:,ibuf) = Zero
    call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IVLEV,IXLEV,One,STSYM,CI,WRK(:,ibuf))
    if (ibuf == NZVEC) then
      !! Finalize the contribution
      call dgemm_('N','T',MXCI,nxy,ibuf,One,WRK,MXCI,Gact_sort,nxy,Scal,ZVEC,MXCI)
      ibuf = 0
      Scal = One
    end if
  end do
  call Free_Tsk(ID)

  if (ibuf > 0) call dgemm_('N','T',MXCI,nxy,ibuf,One,WRK,MXCI,Gact_sort,nxy,Scal,ZVEC,MXCI)

  ! avoid the 2 GB (?) barrier
# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    do ip1=1,NZVEC,MAXBUF
      call GADGOP(ZVEC(:,ip1:),MXCI*min(NZVEC-ip1+1,MAXBUF),'+')
    end do
  end if
# endif

  call mma_deallocate(Gact_sort)

end subroutine NEVPT2_E4_ZVEC

!-----------------------------------------------------------------------

subroutine NEVPT2_E4_XYVEC(iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,ZVEC,XYVEC)

  integer(kind=iwp), intent(in) :: iSym, NLEV, idx2ij(2,NLEV**2), ij2idx(NLEV,NLEV), ipxysta, ipxyend
  real(kind=wp), intent(inout) :: BUFT(:), ZVEC(MXCI,NZVEC), XYVEC(MXCI,NXYVEC,2)
  integer(kind=iwp) :: ID, ip1, ip2, ipxy, issg1, issg2, istu, isvx, it, itlev, iu, iulev, ivlev, ix, ixlev, locx, locy, nlev2, &
                       nsgm1, nTask

  ! Construct XVEC and YVEC
  ! X_{a,b}^I = \sum_{def} (ae|df) \sum_J <I|E_{eb}|J><J|E_{df}|Psi>
  !           = \sum_e \sum_J <I|E_{eb}|J>*(\sum_{df}(ae|df)<J|E_{df}|Psi>)
  !           = \sum_e \sum_J <I|E_{eb}|J>*Z_{ae}^J
  ! X_{v,u}^I = <Phi_I|E_{x,u}|Phi_J> Z_{v,x}^J
  ! Y_{v,t}^I = <Phi_I|E_{t,x}|Phi_J> Z_{v,x}^J

  locx = 1
  locy = 1
  if (do_xvec .and. do_yvec) then
    XYvec(1:MXCI,1:NXYVEC,1:2) = Zero
    locy = 2
  else if (do_xvec .or. do_yvec) then
    XYvec(1:MXCI,1:NXYVEC,1:1) = Zero
    locy = 1
  end if

  NLEV2 = NLEV**2
  nTask = NLEV2
  call Init_Tsk(ID,nTask)
  do while (Rsv_Tsk(ID,ip2))
    ivlev = idx2ij(1,ip2)
    ixlev = idx2ij(2,ip2)
    isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
    issg2 = Mul(isvx,stsym)
    !iv = SGS(istate)%L2ACT(ivlev)
    ix = SGS(istate)%L2ACT(ixlev)
    if ((NXY_work /= NLEV) .and. ((ivlev < ixyzsta) .or. (ivlev > ixyzend))) cycle
    do ip1=1,nlev2
      itlev = idx2ij(1,ip1)
      iulev = idx2ij(2,ip1)
      istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
      issg1 = Mul(istu,issg2)
      if (issg1 /= iSym) cycle
      nsgm1=CIS(istate)%ncsf(issg1)
      it = SGS(istate)%L2ACT(itlev)
      iu = SGS(istate)%L2ACT(iulev)
      if (((.not. do_xvec) .or. (it /= ix)) .and. ((.not. do_yvec) .or. (iu /= ix))) cycle
      BUFT(1:nsgm1) = Zero
      if (NXY_work == NLEV) then
        call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),ITLEV,IULEV,One,issg2,ZVEC(:,ip2),BUFT)
      else
        call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),ITLEV,IULEV,One,issg2,ZVEC(:,ivlev-ixyzsta+1+nxyzdim*(ixlev-1)),BUFT)
      end if
      ! Etu Evx * (at|vx) -> Xau
      if (do_xvec .and. (it == ix)) then
        if (NXY_work == NLEV) then
          ipxy = ij2idx(ivlev,iulev)
          if ((ipxy >= ipxysta) .and. (ipxy <= ipxyend)) &
            XYvec(1:nsgm1,ipxy-ipxysta+1,locx) = XYvec(1:nsgm1,ipxy-ipxysta+1,locx)+buft(1:nsgm1)
        else
          XYvec(1:nsgm1,ivlev-ixyzsta+1+nxyzdim*(iulev-1),locx) = &
            XYvec(1:nsgm1,ivlev-ixyzsta+1+nxyzdim*(iulev-1),locx)+buft(1:nsgm1)
        end if
      end if
      ! Evx Etu * (ax|tu) -> Yav
      if (do_yvec .and. (iu == ix)) then
        if (NXY_work == NLEV) then
          ipxy = ij2idx(ivlev,itlev)
          if ((ipxy >= ipxysta) .and. (ipxy <= ipxyend)) &
            XYvec(1:nsgm1,ipxy-ipxysta+1,locy) = XYvec(1:nsgm1,ipxy-ipxysta+1,locy)+buft(1:nsgm1)
        else
          XYvec(1:nsgm1,ivlev-ixyzsta+1+nxyzdim*(itlev-1),locy) = &
            XYvec(1:nsgm1,ivlev-ixyzsta+1+nxyzdim*(itlev-1),locy)+buft(1:nsgm1)
        end if
      end if
    end do
  end do
  call Free_Tsk(ID)

# ifdef _MOLCAS_MPP_
  ! avoid the 2 GB (?) barrier
  if (is_real_par()) then
    do ip1=1,NXYVEC,MAXBUF
      call GADGOP(XYvec(:,ip1:,1),MXCI*min(NXYVEC-ip1+1,MAXBUF),'+')
    end do
    if (do_xvec .and. do_yvec) then
      do ip1=1,NXYVEC,MAXBUF
        call GADGOP(XYvec(:,ip1:,2),MXCI*min(NXYVEC-ip1+1,MAXBUF),'+')
      end do
    end if
  end if
# endif

end subroutine NEVPT2_E4_XYVEC

!-----------------------------------------------------------------------

subroutine NEVPT2_E4_contract1(iSym0,iSym,NLEV,idx2ij,ipxysta,ipxyend,BUFT,CI,XYVEC,BA,BC,XYcont)

  integer(kind=iwp), intent(in) :: iSym0, iSym, NLEV, idx2ij(2,NLEV**2), ipxysta, ipxyend
  real(kind=wp), intent(inout) :: BUFT(:), BA(:), BC(:)
  real(kind=wp), intent(in) :: CI(:), XYVEC(MXCI,NXYVEC,2)
  real(kind=wp), intent(out) :: XYcont(:,:,:)
  integer(kind=iwp) :: I, ID, ip1, ipxy, issg1, istu, ISUP, isvx, it, itlev, iu, iulev, iv, ivlev, ix, ixlev, iz, izlev, JSUP, &
                       nlev2, nsgm1, nTask, nxy
  real(kind=wp) :: tmp1, tmp2

  if (NXY_work == NLEV) then
    nxy = ipxyend-ipxysta+1
  else
    nxy = NLEV*nxyzdim
  end if
  NLEV2 = NLEV**2
  XYcont(1:nxy,1:2,1:nlev2) = Zero
  nTask = nlev2
  call Init_Tsk(ID,nTask)
  do while (Rsv_Tsk(ID,ip1))
    itlev = idx2ij(1,ip1)
    iulev = idx2ij(2,ip1)
    istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
    issg1 = Mul(istu,stsym)
    if (issg1 /= iSym) cycle
    nsgm1=CIS(istate)%ncsf(issg1)
    it=SGS(istate)%L2ACT(itlev)
    iu=SGS(istate)%L2ACT(iulev)
    BUFT(1:nsgm1) = Zero
    call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,One,STSYM,CI,BUFT(:))
    if (do_xvec .and. do_yvec) then
      call dgemv_('T',nsgm1,nxy,One,XYvec(:,:,1),mxci,BUFT(:),1,Zero,XYcont(:,1,ip1),1)
      call dgemv_('T',nsgm1,nxy,One,XYvec(:,:,2),mxci,BUFT(:),1,Zero,XYcont(:,2,ip1),1)
    else if (do_xvec .and. (.not. do_yvec)) then
      call dgemv_('T',nsgm1,nxy,One,XYvec(:,:,1),mxci,BUFT(:),1,Zero,XYcont(:,1,ip1),1)
    else if ((.not. do_xvec) .and. do_yvec) then
      call dgemv_('T',nsgm1,nxy,One,XYvec(:,:,1),mxci,BUFT(:),1,Zero,XYcont(:,2,ip1),1)
    end if
    do ipxy=ipxysta,ipxyend
      ivlev = idx2ij(1,ipxy)
      ixlev = idx2ij(2,ipxy)
      iv=SGS(istate)%L2ACT(ivlev)
      ix=SGS(istate)%L2ACT(ixlev)
      isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
      ! ----- E*X and E*Y terms
      ! Etu Xvx = Xvx Eut
      if (do_xvec) then
        if (NXY_work == NLEV) then
          tmp1 = XYcont(ipxy-ipxysta+1,1,ip1)
        else
          !if (ivlev == ixy_local) then
          if ((ivlev >= ixyzsta) .and. (ivlev <= ixyzend)) then
            !tmp1 = XYcont(ixlev,1,ip1)
            tmp1 = XYcont(ivlev-ixyzsta+1+nxyzdim*(ixlev-1),1,ip1)
          else
            tmp1 = Zero
          end if
        end if
        if (abs(tmp1) >= 1.0e-12_wp) then
          do izlev=1,nlev
            iz=SGS(istate)%L2ACT(izlev)
            if ((Mul(SGS(istate)%ism(izlev),istu) == iSym0) .and. (Mul(SGS(istate)%ism(izlev),isvx) == iSym0)) then
              ! A: A_{zut,vzx} <-- -del(zz) Etu Xvx = -Xvx Eut
              ISUP = KTUV(iZ,iU,iT)-nTUVES(iSYM0)
              JSUP = KTUV(iV,iZ,iX)-nTUVES(iSYM0)
              if (isup >= jsup) then
                I = iTri(ISUP,JSUP)
                BA(I) = BA(I)-tmp1
              end if
              ! A: A_{zut,zvx} <-- -del(zz) Etu Xvx (E*X from tilde{E})
              ISUP = KTUV(iZ,iU,iT)-nTUVES(iSYM0)
              JSUP = KTUV(iZ,iV,iX)-nTUVES(iSYM0)
              if (isup >= jsup) then
                I = iTri(ISUP,JSUP)
                BA(I) = BA(I)+Two*tmp1
              end if
            end if
          end do
        end if
      end if

      if (do_yvec) then
        if (NXY_work == NLEV) then
          tmp2 = XYcont(ipxy-ipxysta+1,2,ip1)
        else
          !if (ivlev == ixy_local) then
          if ((ivlev >= ixyzsta) .and. (ivlev <= ixyzend)) then
            !tmp2 = XYcont(ixlev,2,ip1)
            tmp2 = XYcont(ivlev-ixyzsta+1+nxyzdim*(ixlev-1),2,ip1)
          else
            tmp2 = Zero
          end if
        end if
        if (abs(tmp2) >= 1.0e-12_wp) then
          do izlev=1,nlev
            iz=SGS(istate)%L2ACT(izlev)
            if ((Mul(SGS(istate)%ism(izlev),istu) == iSym0) .and. (Mul(SGS(istate)%ism(izlev),isvx) == iSym0)) then
              ! C: A_{zut,vxz} <-- del(zz) Etu Yvx = Yvx Eut
              ISUP = KTUV(iZ,iU,iT)-nTUVES(iSYM0)
              JSUP = KTUV(iV,iX,iZ)-nTUVES(iSYM0)
              if (isup >= jsup) then
                I = iTri(ISUP,JSUP)
                BC(I) = BC(I)+tmp2
              end if
              ! A: A_{zut,zxv} <-- -del(zz) Etu Xvx
              ISUP = KTUV(iZ,iU,iT)-nTUVES(iSYM0)
              JSUP = KTUV(iZ,iX,iV)-nTUVES(iSYM0)
              if (isup >= jsup) then
                I = iTri(ISUP,JSUP)
                BA(I) = BA(I)-Two*tmp2
              end if
            end if
          end do
        end if
      end if
    end do
  end do

  call Free_Tsk(ID)

# ifdef _MOLCAS_MPP_
  if (is_real_par()) call GADGOP(XYcont,nxyvec*2*nlev2,'+')
# endif

end subroutine NEVPT2_E4_contract1

!-----------------------------------------------------------------------

subroutine NEVPT2_E4_contract2(iSym0,iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,CI,XYVEC,BA,BC,XYcont,XYtmp,ZVEC)

  integer(kind=iwp), intent(in) :: iSym0, iSym, NLEV, idx2ij(2,NLEV**2), ij2idx(NLEV,NLEV), ipxysta, ipxyend
  real(kind=wp), intent(inout) :: BUFT(:), BA(:), BC(:)
  real(kind=wp), intent(in) :: CI(:), XYVEC(MXCI,NXYVEC,2), XYcont(:,:,:)
  real(kind=wp), intent(out) :: XYtmp(:,:,:), ZVEC(:,:)
  integer(kind=iwp) :: I, ibufxy, ID, ip1, ip2, ipxy, issg1, issg2, istu, ISUP, isvx, isyz, it, itlev, iu, iulev, iv, ivlev, ix, &
                       ixlev, iy, iylev, iz, izlev, JSUP, nbufxy, nlev2, nsgm1, nsgm2, nTask, nxy
  real(kind=wp) :: tmp1, tmp2
  logical(kind=iwp) :: do_dgemm

  if (NXY_work == NLEV) then
    nxy = ipxyend-ipxysta+1
  else
    !nxy = NLEV
    nxy = NLEV*nxyzdim
  end if
  NLEV2 = NLEV*NLEV
  nsgm2 = 0

  nTask = nlev2
  call Init_Tsk(ID,nTask)
  do while (Rsv_Tsk(ID,ip1))
    iylev = idx2ij(1,ip1)
    izlev = idx2ij(2,ip1)
    isyz=Mul(SGS(istate)%ism(iylev),SGS(istate)%ism(izlev))
    issg1 = Mul(isyz,stsym)
    nsgm1=CIS(istate)%ncsf(issg1)
    iy=SGS(istate)%L2ACT(iylev)
    iz=SGS(istate)%L2ACT(izlev)
    BUFT(1:nsgm1) = Zero
    call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IYLEV,IZLEV,One,STSYM,CI,BUFT(:))

    if (NXY_work == NLEV) then
      Zvec(1:MXCI,1:nlev2-ip1+1) = Zero
    else
      Zvec(1:MXCI,1:nxy) = Zero
    end if
    !if (NXY_work == 1) XYtmp(:,:,:) = Zero
    !Zvec(1:MXCI,1:nlev2-ip1+1) = Zero
    do_dgemm = .false.
    nbufxy = 0
    ibufxy = ip1
    do ip2=ip1,nlev2
      ivlev = idx2ij(1,ip2)
      ixlev = idx2ij(2,ip2)
      isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
      issg2 = Mul(isvx,issg1) ! symmetry of <I|EvxEyz|Psi>
      nbufxy = nbufxy+1
      if (issg2 /= iSym) then
        if ((NXY_work /= NLEV) .and. (nbufxy == nlev)) then
          nsgm2=CIS(istate)%ncsf(iSym)
          if (do_xvec .and. do_yvec) then
            call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,1),NXYVEC)
            call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,2),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,2),NXYVEC)
          else if (do_xvec .and. (.not. do_yvec)) then
            call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,1),NXYVEC)
          else if ((.not. do_xvec) .and. do_yvec) then
            call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,2),NXYVEC)
          end if
          nbufxy = 0
          ibufxy = ip2+1
          Zvec(1:MXCI,1:nxy) = Zero
        end if
        cycle
      end if
      nsgm2=CIS(istate)%ncsf(issg2)
      iv=SGS(istate)%L2ACT(ivlev)
      ix=SGS(istate)%L2ACT(ixlev)
      if (NXY_work == NLEV) then
        call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IVLEV,IXLEV,One,issg1,BUFT(:),ZVEC(1:nsgm2,ip2-ip1+1))
      else
        call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IVLEV,IXLEV,One,issg1,BUFT(:),ZVEC(1:nsgm2,nbufxy))
        !if (do_xvec .and. do_yvec) then
        !  call dgemv_('T',nsgm2,nxy,One,XYvec(:,:,1),mxci,ZVEC(:,1),1,Zero,XYtmp(:,ip2-ip1+1,1),1)
        !  call dgemv_('T',nsgm2,nxy,One,XYvec(:,:,2),mxci,ZVEC(:,1),1,Zero,XYtmp(:,ip2-ip1+1,2),1)
        !else if (do_xvec .and. (.not. do_yvec)) then
        !  call dgemv_('T',nsgm2,nxy,One,XYvec(:,:,1),mxci,ZVEC(:,1),1,Zero,XYtmp(:,ip2-ip1+1,1),1)
        !else if ((.not. do_xvec) .and. do_yvec) then
        !  call dgemv_('T',nsgm2,nxy,One,XYvec(:,:,1),mxci,ZVEC(:,1),1,Zero,XYtmp(:,ip2-ip1+1,2),1)
        !end if
        if (nbufxy == nlev) then
          if (do_xvec .and. do_yvec) then
            call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,1),NXYVEC)
            call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,2),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,2),NXYVEC)
          else if (do_xvec .and. (.not. do_yvec)) then
            call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,1),NXYVEC)
          else if ((.not. do_xvec) .and. do_yvec) then
            call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,2),NXYVEC)
          end if
          nbufxy = 0
          ibufxy = ip2+1
          Zvec(1:MXCI,1:nxy) = Zero
        end if
      end if
      do_dgemm = .true.
    end do

    if (NXY_work == NLEV) then
      if ((nlev2-ip1+1 == 0) .or. (nsgm2 == 0) .or. (.not. do_dgemm)) then
        XYtmp(:,:,:) = Zero
      else
        if (do_xvec .and. do_yvec) then
          call dgemm_('T','N',nxy,nlev2-ip1+1,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,:,1),NXYVEC)
          call dgemm_('T','N',nxy,nlev2-ip1+1,nsgm2,One,XYvec(:,:,2),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,:,2),NXYVEC)
        else if (do_xvec .and. (.not. do_yvec)) then
          call dgemm_('T','N',nxy,nlev2-ip1+1,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,:,1),NXYVEC)
        else if ((.not. do_xvec) .and. do_yvec) then
          call dgemm_('T','N',nxy,nlev2-ip1+1,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,:,2),NXYVEC)
        end if
      end if
    else if ((NXY_work /= NLEV) .and. (nbufxy > 0)) then
      if (do_xvec .and. do_yvec) then
        call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,1),NXYVEC)
        call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,2),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,2),NXYVEC)
      else if (do_xvec .and. (.not. do_yvec)) then
        call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,1),NXYVEC)
      else if ((.not. do_xvec) .and. do_yvec) then
        call dgemm_('T','N',nxy,nbufxy,nsgm2,One,XYvec(:,:,1),mxci,ZVEC(:,:),mxci,Zero,XYtmp(:,ibufxy-ip1+1,2),NXYVEC)
      end if
    end if
    do ip2=ip1,nlev2
      ivlev = idx2ij(1,ip2)
      ixlev = idx2ij(2,ip2)
      isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
      issg2 = Mul(isvx,issg1) ! symmetry of <I|EvxEyz|Psi>
      if (issg2 /= iSym) cycle
      nsgm2=CIS(istate)%ncsf(issg2)
      iv=SGS(istate)%L2ACT(ivlev)
      ix=SGS(istate)%L2ACT(ixlev)
      do ipxy=ipxysta,ipxyend
        itlev = idx2ij(1,ipxy)
        iulev = idx2ij(2,ipxy)
        istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
        it=SGS(istate)%L2ACT(itlev)
        iu=SGS(istate)%L2ACT(iulev)
        ! ----- E*E*X terms
        if (do_xvec) then
          if (NXY_work == NLEV) then
            tmp1 = XYtmp(ipxy-ipxysta+1,ip2-ip1+1,1)
          else
            !if (itlev == ixy_local) then
            if ((itlev >= ixyzsta) .and. (itlev <= ixyzend)) then
              !tmp1 = XYtmp(iulev,ip2-ip1+1,1)
              tmp1 = XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip1+1,1)
            else
              tmp1 = Zero
            end if
          end if
          if (abs(tmp1) >= 1.0e-12_wp) call Add_MKBNEVAC_Xvec(iv,ix,iy,iz,ivlev,ixlev,isvx,isyz)
        end if
        ! ----- E*E*Y terms
        if (do_yvec) then
          if (NXY_work == NLEV) then
            tmp2 = XYtmp(ipxy-ipxysta+1,ip2-ip1+1,2)
          else! if (NXY_work == 1) then
            !if (itlev == ixy_local) then
            if ((itlev >= ixyzsta) .and. (itlev <= ixyzend)) then
              !tmp2 = XYtmp(iulev,ip2-ip1+1,2)
              tmp2 = XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip1+1,2)
            else
              tmp2 = Zero
            end if
          end if
          if (abs(tmp2) >= 1.0e-12_wp) call Add_MKBNEVAC_Yvec(iv,ix,iy,iz,ivlev,ixlev,isvx,isyz)
        end if
      end do

      if (ip2 /= ip1) then
        ! EyzEvx = EvxEyz + Eyx del(vz) - Evz del(yx)
        ! Note that the index of ipxy for XYcont is reversed
        if ((iv == iz) .and. (iSym == Mul(Mul(SGS(istate)%ism(iylev),SGS(istate)%ism(ixlev)),STSYM))) then
          ipxy = ij2idx(ixlev,iylev)
          if (do_xvec) XYtmp(1:nxy,ip2-ip1+1,1) = XYtmp(1:nxy,ip2-ip1+1,1)+XYcont(1:nxy,1,ipxy)
          if (do_yvec) XYtmp(1:nxy,ip2-ip1+1,2) = XYtmp(1:nxy,ip2-ip1+1,2)+XYcont(1:nxy,2,ipxy)
        end if
        if ((iy == ix) .and. (iSym == Mul(Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(izlev)),STSYM))) then
          ipxy = ij2idx(izlev,ivlev)
          if (do_xvec) XYtmp(1:nxy,ip2-ip1+1,1) = XYtmp(1:nxy,ip2-ip1+1,1)-XYcont(1:nxy,1,ipxy)
          if (do_yvec) XYtmp(1:nxy,ip2-ip1+1,2) = XYtmp(1:nxy,ip2-ip1+1,2)-XYcont(1:nxy,2,ipxy)
        end if

        do ipxy=ipxysta,ipxyend
          itlev = idx2ij(1,ipxy)
          iulev = idx2ij(2,ipxy)
          istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
          it=SGS(istate)%L2ACT(itlev)
          iu=SGS(istate)%L2ACT(iulev)
          ! ----- E*E*X terms
          if (do_xvec) then
            if (NXY_work == NLEV) then
              tmp1 = XYtmp(ipxy-ipxysta+1,ip2-ip1+1,1)
            else! if (NXY_work == 1) then
              !if (itlev == ixy_local) then
              if ((itlev >= ixyzsta) .and. (itlev <= ixyzend)) then
                !tmp1 = XYtmp(iulev,ip2-ip1+1,1)
                tmp1 = XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip1+1,1)
              else
                tmp1 = Zero
              end if
            end if
            if (abs(tmp1) >= 1.0e-12_wp) call Add_MKBNEVAC_Xvec(iy,iz,iv,ix,iylev,izlev,isyz,isvx)
          end if
          ! ----- E*E*Y terms
          if (do_yvec) then
            if (NXY_work == NLEV) then
              tmp2 = XYtmp(ipxy-ipxysta+1,ip2-ip1+1,2)
            else! if (NXY_work == 1) then
              !if (itlev == ixy_local) then
              if ((itlev >= ixyzsta) .and. (itlev <= ixyzend)) then
                !tmp2 = XYtmp(iulev,ip2-ip1+1,2)
                tmp2 = XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip1+1,2)
              else
                tmp2 = Zero
              end if
            end if
            if (abs(tmp2) >= 1.0e-12_wp) call Add_MKBNEVAC_Yvec(iy,iz,iv,ix,iylev,izlev,isyz,isvx)
          end if
        end do
      end if
    end do
  end do

  call Free_Tsk(ID)

  contains

  subroutine Add_MKBNEVAC_Xvec(iv_,ix_,iy_,iz_,ivlev_,ixlev_,isvx_,isyz_)

    integer(kind=iwp), intent(in) :: iv_, ix_, iy_, iz_, ivlev_, ixlev_, isvx_, isyz_

    if ((Mul(SGS(istate)%ism(ixlev_),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(ivlev_),istu) == iSym0)) then
      ! C: A_{xyz,vtu} <-- Xtu Evx Eyz
      ISUP = KTUV(iX_,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iV_,iT,iU)-nTUVES(iSYM0)
      if (isup >= jsup) then
        I = iTri(ISUP,JSUP)
        BC(I) = BC(I)+tmp1
      end if
    end if
    if ((Mul(SGS(istate)%ism(ivlev_),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(ixlev_),istu) == iSym0)) then
      ! A: A_{vyz,xtu} <-- -Xtu Evx Eyz
      ISUP = KTUV(iV_,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iX_,iT,iU)-nTUVES(iSYM0)
      if (isup >= jsup) then
        I = iTri(ISUP,JSUP)
        BA(I) = BA(I)-tmp1
      end if
    end if
    if ((Mul(SGS(istate)%ism(iulev),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(itlev),isvx_) == iSym0)) then
      ! A: A_{uyz,txv} <-- -Xtu Evx Eyz
      ISUP = KTUV(iU,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iT,iX_,iV_)-nTUVES(iSYM0)
      if (isup >= jsup) then
        I = iTri(ISUP,JSUP)
        BA(I) = BA(I)-tmp1
      end if
    end if

    return

  end subroutine Add_MKBNEVAC_Xvec

  subroutine Add_MKBNEVAC_Yvec(iv_,ix_,iy_,iz_,ivlev_,ixlev_,isvx_,isyz_)

    integer(kind=iwp), intent(in) :: iv_, ix_, iy_, iz_, ivlev_, ixlev_, isvx_, isyz_

    if ((Mul(SGS(istate)%ism(ixlev_),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(ivlev_),istu) == iSym0)) then
      ! C: A_{xyz,vut} <-- -Ytu Evx Eyz
      ISUP = KTUV(iX_,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iV_,iU,iT)-nTUVES(iSYM0)
      if (isup >= jsup) then
        I = iTri(ISUP,JSUP)
        BC(I) = BC(I)-tmp2
      end if
    end if
    if ((Mul(SGS(istate)%ism(iulev),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(itlev),isvx_) == iSym0)) then
      ! C: A_{uyz,txv} <-- -Ytu Evx Eyz
      ISUP = KTUV(iU,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iT,iX_,iV_)-nTUVES(iSYM0)
      if (isup >= jsup) then
        I = iTri(ISUP,JSUP)
        BC(I) = BC(I)-tmp2
      end if
    end if
    if ((Mul(SGS(istate)%ism(ivlev_),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(ixlev_),istu) == iSym0)) then
      ! A: A_{vyz,xut} <-- +Xtu Evx Eyz
      ISUP = KTUV(iV_,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iX_,iU,iT)-nTUVES(iSYM0)
      if (isup >= jsup) then
        I = iTri(ISUP,JSUP)
        BA(I) = BA(I)+tmp2
      end if
    end if

    return

  end subroutine Add_MKBNEVAC_Yvec

end subroutine NEVPT2_E4_contract2

!-----------------------------------------------------------------------

subroutine NEVPT2_E4_derivative1(iSym0,iSym,NLEV,idx2ij,ipxysta,ipxyend,BDERA,BDERC,XYcont)

  integer(kind=iwp), intent(in) :: iSym0, iSym, NLEV, idx2ij(2,NLEV**2), ipxysta, ipxyend
  real(kind=wp), intent(in) :: BDERA(:,:), BDERC(:,:)
  real(kind=wp), intent(out) :: XYcont(:,:,:)
  integer(kind=iwp) :: ID, ip1, ipxy, issg1, istu, ISUP, isvx, it, itlev, iu, iulev, iv, ivlev, ix, ixlev, iz, izlev, JSUP, nlev2, &
                       nTask, nxy
  real(kind=wp) :: tmp1, tmp2

  if (NXY_work == NLEV) then
    nxy = ipxyend-ipxysta+1
  else
    nxy = NLEV*nxyzdim
  end if
  NLEV2 = NLEV**2
  XYcont(1:nxy,1:2,1:nlev2) = Zero

  nTask = nlev2
  call Init_Tsk(ID,nTask)
  do while (Rsv_Tsk(ID,ip1))
    itlev = idx2ij(1,ip1)
    iulev = idx2ij(2,ip1)
    istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
    issg1 = Mul(istu,stsym)
    if (issg1 /= iSym) cycle
    it=SGS(istate)%L2ACT(itlev)
    iu=SGS(istate)%L2ACT(iulev)
    do ipxy=ipxysta,ipxyend
      ivlev = idx2ij(1,ipxy)
      ixlev = idx2ij(2,ipxy)
      iv=SGS(istate)%L2ACT(ivlev)
      ix=SGS(istate)%L2ACT(ixlev)
      isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
      if ((NXY_work /= NLEV) .and. ((ivlev < ixyzsta) .or. (ivlev > ixyzend))) cycle
      ! ----- E*X and E*Y terms
      if (do_xvec) then
        tmp1 = Zero
        do izlev=1,nlev
          iz = SGS(istate)%L2ACT(izlev)
          if ((Mul(SGS(istate)%ism(izlev),istu) == iSym0) .and. (Mul(SGS(istate)%ism(izlev),isvx) == iSym0)) then
            ! A: A_{zut,vzx} <-- -del(zz) Etu Xvx = -Xvx Eut
            ISUP = KTUV(iZ,iU,iT)-nTUVES(iSYM0)
            JSUP = KTUV(iV,iZ,iX)-nTUVES(iSYM0)
            tmp1 = tmp1-BDERA(ISUP,JSUP)
            ! A: A_{zut,zvx} <-- -del(zz) Etu Xvx (E*X from tilde{E})
            ISUP = KTUV(iZ,iU,iT)-nTUVES(iSYM0)
            JSUP = KTUV(iZ,iV,iX)-nTUVES(iSYM0)
            tmp1 = tmp1+Two*BDERA(ISUP,JSUP)
          end if
        end do
        if (NXY_work == NLEV) then
          XYcont(ipxy-ipxysta+1,1,ip1) = tmp1
        else if ((ivlev >= ixyzsta) .and. (ivlev <= ixyzend)) then
          XYcont(ivlev-ixyzsta+1+nxyzdim*(ixlev-1),1,ip1) = tmp1
        end if
      end if

      if (do_yvec) then
        tmp2 = Zero
        do izlev=1,nlev
          iz = SGS(istate)%L2ACT(izlev)
          if ((Mul(SGS(istate)%ism(izlev),istu) == iSym0) .and. (Mul(SGS(istate)%ism(izlev),isvx) == iSym0)) then
            ! C: A_{zut,vxz} <-- del(zz) Etu Yvx = Yvx Eut
            ISUP = KTUV(iZ,iU,iT)-nTUVES(iSYM0)
            JSUP = KTUV(iV,iX,iZ)-nTUVES(iSYM0)
            tmp2 = tmp2+BDERC(ISUP,JSUP)
            ! A: A_{zut,zxv} <-- -del(zz) Etu Xvx
            ISUP = KTUV(iZ,iU,iT)-nTUVES(iSYM0)
            JSUP = KTUV(iZ,iX,iV)-nTUVES(iSYM0)
            tmp2 = tmp2-Two*BDERA(ISUP,JSUP)
          end if
        end do
        if (NXY_work == NLEV) then
          XYcont(ipxy-ipxysta+1,2,ip1) = tmp2
        else if ((ivlev >= ixyzsta) .and. (ivlev <= ixyzend)) then
          XYcont(ivlev-ixyzsta+1+nxyzdim*(ixlev-1),2,ip1) = tmp2
        end if
      end if
    end do
  end do

  call Free_Tsk(ID)

end subroutine NEVPT2_E4_derivative1

!-----------------------------------------------------------------------

subroutine NEVPT2_E4_derivative2(iSym0,iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,CI,BDERA,BDERC,XYcont,XYcontder,XYder,ZVEC, &
                                 XYtmp)

  integer(kind=iwp), intent(in) :: iSym0, iSym, NLEV, idx2ij(2,NLEV**2), ij2idx(NLEV,NLEV), ipxysta, ipxyend
  real(kind=wp), intent(inout) :: BUFT(:), XYcont(:,:,:), XYcontder(:,:,:), XYder(:,:,:), ZVEC(:,:), XYtmp(:,:,:)
  real(kind=wp), intent(in) :: CI(:), BDERA(:,:), BDERC(:,:)
  integer(kind=iwp) :: ibufxy, ID, ip1, ip2, ip23, ip3, ipxy, issg2, issg3, istu, ISUP, isvx, isyz, it, itlev, iu, iulev, iv, &
                       ivlev, ix, ixlev, iy, iylev, iz, izlev, JSUP, nbufxy, nlev2, nsgm2, nsgm3, nTask, nxy
  real(kind=wp) :: tmp1, tmp2

  if (NXY_work == NLEV) then
    nxy = ipxyend-ipxysta+1
  else
    nxy = NLEV*nxyzdim
  end if
  NLEV2 = NLEV**2

  if (do_xvec .and. do_yvec) then
    XYder(1:mxci,1:nxy,1:2) = Zero
  else if (do_xvec .or. do_yvec) then
    XYder(1:mxci,1:nxy,1:1) = Zero
  end if

  XYcontder(1:nxyvec,1:2,1:nTri_Elem(nlev2)) = Zero
  nTask = nlev2
  call Init_Tsk(ID,nTask)
  do while (Rsv_Tsk(ID,ip3))
    iylev = idx2ij(1,ip3)
    izlev = idx2ij(2,ip3)
    isyz=Mul(SGS(istate)%ism(iylev),SGS(istate)%ism(izlev))
    issg3 = Mul(isyz,stsym)
    nsgm3=CIS(istate)%ncsf(issg3)
    iy=SGS(istate)%L2ACT(iylev)
    iz=SGS(istate)%L2ACT(izlev)
    XYtmp(1:nxy,1:nlev2-ip3+1,1:2) = Zero
    do ip2=ip3,nlev2
      ivlev = idx2ij(1,ip2)
      ixlev = idx2ij(2,ip2)
      isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
      issg2 = Mul(isvx,issg3) ! symmetry of <I|EvxEyz|Psi>
      if (issg2 /= iSym) cycle
      nsgm2=CIS(istate)%ncsf(issg2)
      iv=SGS(istate)%L2ACT(ivlev)
      ix=SGS(istate)%L2ACT(ixlev)

      do ipxy=ipxysta,ipxyend
        itlev = idx2ij(1,ipxy)
        iulev = idx2ij(2,ipxy)
        istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
        it = SGS(istate)%L2ACT(itlev)
        iu = SGS(istate)%L2ACT(iulev)
        if ((NXY_work /= NLEV) .and. ((itlev < ixyzsta) .or. (itlev > ixyzend))) cycle
        ! ----- E*E*X terms
        if (do_xvec) then
          call Add_MKBNEVAC_Xvec(iv,ix,iy,iz,ivlev,ixlev,isvx,isyz,tmp1)
          if (NXY_work == NLEV) then
            XYtmp(ipxy-ipxysta+1,ip2-ip3+1,1) = XYtmp(ipxy-ipxysta+1,ip2-ip3+1,1)+tmp1
          else
            XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip3+1,1) = XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip3+1,1)+tmp1
          end if
        end if
        ! ----- E*E*Y terms
        if (do_yvec) then
          call Add_MKBNEVAC_Yvec(iv,ix,iy,iz,ivlev,ixlev,isvx,isyz,tmp2)
          if (NXY_work == NLEV) then
            XYtmp(ipxy-ipxysta+1,ip2-ip3+1,2) = XYtmp(ipxy-ipxysta+1,ip2-ip3+1,2)+tmp2
          else
            XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip3+1,2) = XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip3+1,2)+tmp2
          end if
        end if
      end do
      if (ip2 /= ip3) then
        ! EyzEvx = EvxEyz + Eyx del(vz) - Evz del(yx)
        ! Note that the index of ipxy for XYcont is reversed
        do ipxy=ipxysta,ipxyend
          itlev = idx2ij(1,ipxy)
          iulev = idx2ij(2,ipxy)
          istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
          it=SGS(istate)%L2ACT(itlev)
          iu=SGS(istate)%L2ACT(iulev)
          if ((NXY_work /= NLEV) .and. ((itlev < ixyzsta) .or. (itlev > ixyzend))) cycle
          ! ----- E*E*X terms
          if (do_xvec) then
            call Add_MKBNEVAC_Xvec(iy,iz,iv,ix,iylev,izlev,isyz,isvx,tmp1)
            if (NXY_work == NLEV) then
              XYtmp(ipxy-ipxysta+1,ip2-ip3+1,1) = XYtmp(ipxy-ipxysta+1,ip2-ip3+1,1)+tmp1
            else
              XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip3+1,1) = XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip3+1,1)+tmp1
            end if
          end if
          ! ----- E*E*Y terms
          if (do_yvec) then
            call Add_MKBNEVAC_Yvec(iy,iz,iv,ix,iylev,izlev,isyz,isvx,tmp2)
            if (NXY_work == NLEV) then
              XYtmp(ipxy-ipxysta+1,ip2-ip3+1,2) = XYtmp(ipxy-ipxysta+1,ip2-ip3+1,2)+tmp2
            else
              XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip3+1,2) = XYtmp(itlev-ixyzsta+1+nxyzdim*(iulev-1),ip2-ip3+1,2)+tmp2
            end if
          end if
          if ((iv == iz) .and. (iSym == Mul(Mul(SGS(istate)%ism(iylev),SGS(istate)%ism(ixlev)),STSYM))) then
            ip1 = ij2idx(ixlev,iylev)
            if (NXY_work == NLEV) then
              if (do_xvec) XYcont(ipxy-ipxysta+1,1,ip1) = XYcont(ipxy-ipxysta+1,1,ip1)+tmp1
              if (do_yvec) XYcont(ipxy-ipxysta+1,2,ip1) = XYcont(ipxy-ipxysta+1,2,ip1)+tmp2
            else if ((itlev >= ixyzsta) .and. (itlev <= ixyzend)) then
              if (do_xvec) XYcont(itlev-ixyzsta+1+nxyzdim*(iulev-1),1,ip1) = XYcont(itlev-ixyzsta+1+nxyzdim*(iulev-1),1,ip1)+tmp1
              if (do_yvec) XYcont(itlev-ixyzsta+1+nxyzdim*(iulev-1),2,ip1) = XYcont(itlev-ixyzsta+1+nxyzdim*(iulev-1),2,ip1)+tmp2
            end if
          end if
          if ((iy == ix) .and. (iSym == Mul(Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(izlev)),STSYM))) then
            ip1 = ij2idx(izlev,ivlev)
            if (NXY_work == NLEV) then
              if (do_xvec) XYcont(ipxy-ipxysta+1,1,ip1) = XYcont(ipxy-ipxysta+1,1,ip1)-tmp1
              if (do_yvec) XYcont(ipxy-ipxysta+1,2,ip1) = XYcont(ipxy-ipxysta+1,2,ip1)-tmp2
            else if ((itlev >= ixyzsta) .and. (itlev <= ixyzend)) then
              if (do_xvec) XYcont(itlev-ixyzsta+1+nxyzdim*(iulev-1),1,ip1) = XYcont(itlev-ixyzsta+1+nxyzdim*(iulev-1),1,ip1)-tmp1
              if (do_yvec) XYcont(itlev-ixyzsta+1+nxyzdim*(iulev-1),2,ip1) = XYcont(itlev-ixyzsta+1+nxyzdim*(iulev-1),2,ip1)-tmp2
            end if
          end if
        end do
      end if
    end do

    BUFT(1:nsgm3) = Zero
    call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IYLEV,IZLEV,One,STSYM,CI,BUFT(:))

    if (NXY_work == NLEV) then
      Zvec(1:MXCI,1:nlev2-ip3+1) = Zero
    else
      Zvec(1:MXCI,1:nxy) = Zero
    end if
    nbufxy = 0
    ibufxy = ip3
    do ip2=ip3,nlev2
      ivlev = idx2ij(1,ip2)
      ixlev = idx2ij(2,ip2)
      isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
      issg2 = Mul(isvx,issg3) ! symmetry of <I|EvxEyz|Psi>
      nsgm2=CIS(istate)%ncsf(issg2)
      nbufxy = nbufxy+1
      !! see contract2, if symmetries are to be actived
      if (issg2 /= iSym) cycle
      iv=SGS(istate)%L2ACT(ivlev)
      ix=SGS(istate)%L2ACT(ixlev)
      if (NXY_work == NLEV) then
        call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IVLEV,IXLEV,One,issg3,BUFT(:),ZVEC(1:nsgm2,ip2-ip3+1))
      else
        call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IVLEV,IXLEV,One,issg3,BUFT(:),ZVEC(1:nsgm2,nbufxy))
        if (nbufxy == nlev) then
          if (do_xvec .and. do_yvec) then
            call dgemm_('N','T',MXCI,nxy,nbufxy,One,ZVEC(:,:),MXCI,XYtmp(:,ibufxy-ip3+1,1),NXYVEC,One,XYder(:,:,1),MXCI)
            call dgemm_('N','T',MXCI,nxy,nbufxy,One,ZVEC(:,:),MXCI,XYtmp(:,ibufxy-ip3+1,2),NXYVEC,One,XYder(:,:,2),MXCI)
          else if (do_xvec .and. (.not. do_yvec)) then
            call dgemm_('N','T',MXCI,nxy,nbufxy,One,ZVEC(:,:),MXCI,XYtmp(:,ibufxy-ip3+1,1),NXYVEC,One,XYder(:,:,1),MXCI)
          else if ((.not. do_xvec) .and. do_yvec) then
            call dgemm_('N','T',MXCI,nxy,nbufxy,One,ZVEC(:,:),MXCI,XYtmp(:,ibufxy-ip3+1,2),NXYVEC,One,XYder(:,:,1),MXCI)
          end if
          nbufxy = 0
          ibufxy = ip2+1
          Zvec(1:MXCI,1:nxy) = Zero
        end if
      end if
      !! Compute Xder and Yder and left pre-derivative (something like <I|Evx|J>*Xyz^J)
      ip23 = iTri(ip2,ip3)
      if (do_xvec) XYcontder(1:nxy,1,ip23) = XYcontder(1:nxy,1,ip23)+XYtmp(1:nxy,ip2-ip3+1,1)
      if (do_yvec) XYcontder(1:nxy,2,ip23) = XYcontder(1:nxy,2,ip23)+XYtmp(1:nxy,ip2-ip3+1,2)
    end do

    if (nbufxy > 0) then
      if (do_xvec .and. do_yvec) then
        call dgemm_('N','T',MXCI,nxy,nbufxy,One,ZVEC(:,:),MXCI,XYtmp(:,ibufxy-ip3+1,1),NXYVEC,One,XYder(:,:,1),MXCI)
        call dgemm_('N','T',MXCI,nxy,nbufxy,One,ZVEC(:,:),MXCI,XYtmp(:,ibufxy-ip3+1,2),NXYVEC,One,XYder(:,:,2),MXCI)
      else if (do_xvec .and. (.not. do_yvec)) then
        call dgemm_('N','T',MXCI,nxy,nbufxy,One,ZVEC(:,:),MXCI,XYtmp(:,ibufxy-ip3+1,1),NXYVEC,One,XYder(:,:,1),MXCI)
      else if ((.not. do_xvec) .and. do_yvec) then
        call dgemm_('N','T',MXCI,nxy,nbufxy,One,ZVEC(:,:),MXCI,XYtmp(:,ibufxy-ip3+1,2),NXYVEC,One,XYder(:,:,1),MXCI)
      end if
    end if
  end do

  call Free_Tsk(ID)

# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    call GADGOP(XYcont,nxyvec*2*nlev2,'+')
    call GADGOP(XYcontder,nxyvec*2*nTri_Elem(nlev2),'+')
  end if
# endif

 contains

  subroutine Add_MKBNEVAC_Xvec(iv_,ix_,iy_,iz_,ivlev_,ixlev_,isvx_,isyz_,output)

    integer(kind=iwp), intent(in) :: iv_, ix_, iy_, iz_, ivlev_, ixlev_, isvx_, isyz_
    real(kind=wp), intent(out) :: output
    real(kind=wp) :: val

    val = Zero

    if ((Mul(SGS(istate)%ism(ixlev_),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(ivlev_),istu) == iSym0)) then
      ! C: A_{xyz,vtu} <-- Xtu Evx Eyz
      ISUP = KTUV(iX_,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iV_,iT,iU)-nTUVES(iSYM0)
      val = val+BDERC(ISUP,JSUP)
    end if
    if ((Mul(SGS(istate)%ism(ivlev_),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(ixlev_),istu) == iSym0)) then
      ! A: A_{vyz,xtu} <-- -Xtu Evx Eyz
      ISUP = KTUV(iV_,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iX_,iT,iU)-nTUVES(iSYM0)
      val = val-BDERA(ISUP,JSUP)
    end if
    if ((Mul(SGS(istate)%ism(iulev),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(itlev),isvx_) == iSym0)) then
      ! A: A_{uyz,txv} <-- -Xtu Evx Eyz
      ISUP = KTUV(iU,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iT,iX_,iV_)-nTUVES(iSYM0)
      val = val-BDERA(ISUP,JSUP)
    end if

    output = val

    return

  end subroutine Add_MKBNEVAC_Xvec

  subroutine Add_MKBNEVAC_Yvec(iv_,ix_,iy_,iz_,ivlev_,ixlev_,isvx_,isyz_,output)

    integer(kind=iwp), intent(in) :: iv_, ix_, iy_, iz_, ivlev_, ixlev_, isvx_, isyz_
    real(kind=wp), intent(out) :: output
    real(kind=wp) :: val

    val = Zero

    if ((Mul(SGS(istate)%ism(ixlev_),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(ivlev_),istu) == iSym0)) then
      ! C: A_{xyz,vut} <-- -Ytu Evx Eyz
      ISUP = KTUV(iX_,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iV_,iU,iT)-nTUVES(iSYM0)
      val = val-BDERC(ISUP,JSUP)
    end if
    if ((Mul(SGS(istate)%ism(iulev),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(itlev),isvx_) == iSym0)) then
      ! C: A_{uyz,txv} <-- -Ytu Evx Eyz
      ISUP = KTUV(iU,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iT,iX_,iV_)-nTUVES(iSYM0)
      val = val-BDERC(ISUP,JSUP)
    end if
    if ((Mul(SGS(istate)%ism(ivlev_),isyz_) == iSym0) .and. (Mul(SGS(istate)%ism(ixlev_),istu) == iSym0)) then
      ! A: A_{vyz,xut} <-- +Xtu Evx Eyz
      ISUP = KTUV(iV_,iY_,iZ_)-nTUVES(iSYM0)
      JSUP = KTUV(iX_,iU,iT)-nTUVES(iSYM0)
      val = val+BDERA(ISUP,JSUP)
    end if

    output = val

    return

  end subroutine Add_MKBNEVAC_Yvec

end subroutine NEVPT2_E4_derivative2

!-----------------------------------------------------------------------

subroutine NEVPT2_E4_derivative3(iSym,NLEV,idx2ij,ipxysta,ipxyend,BUFT,CI,XYcont,XYder)

  integer(kind=iwp), intent(in) :: iSym, NLEV, idx2ij(2,NLEV**2), ipxysta, ipxyend
  real(kind=wp), intent(inout) :: BUFT(:), XYcont(:,:,:), XYder(:,:,:)
  real(kind=wp), intent(in) :: CI(:)
  integer(kind=iwp) :: ID, ip1, ipxy, issg1, istu, itlev, iulev, nbufxy, nlev2, nsgm1, nTask, nxy

  if (NXY_work == NLEV) then
    nxy = ipxyend-ipxysta+1
  else
    nxy = NLEV*nxyzdim
  end if
  NLEV2 = NLEV**2

  nTask = nlev2
  call Init_Tsk(ID,nTask)
  nbufxy = 0
  do while (Rsv_Tsk(ID,ip1))
    itlev = idx2ij(1,ip1)
    iulev = idx2ij(2,ip1)
    istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
    issg1 = Mul(istu,stsym)
    nbufxy = nbufxy+1
    if (issg1 /= iSym) cycle
    nsgm1=CIS(istate)%ncsf(issg1)
!   it=SGS(istate)%L2ACT(itlev)
!   iu=SGS(istate)%L2ACT(iulev)
    BUFT(1:nsgm1) = Zero
    call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,One,STSYM,CI,BUFT(:))
    do ipxy=1,nxy
      ! ----- E*X and E*Y terms
      if (do_xvec .and. do_yvec) then
        XYder(:,ipxy,1) = XYder(:,ipxy,1)+XYcont(ipxy,1,ip1)*BUFT(:)
        XYder(:,ipxy,2) = XYder(:,ipxy,2)+XYcont(ipxy,2,ip1)*BUFT(:)
      else if (do_xvec .and. (.not. do_yvec)) then
        XYder(:,ipxy,1) = XYder(:,ipxy,1)+XYcont(ipxy,1,ip1)*BUFT(:)
      else if ((.not. do_xvec) .and. do_yvec) then
        XYder(:,ipxy,1) = XYder(:,ipxy,1)+XYcont(ipxy,2,ip1)*BUFT(:)
      end if
    end do
  end do

  call Free_Tsk(ID)

# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    do ip1=1,NXYVEC,MAXBUF
      call GADGOP(XYder(:,ip1,1),MXCI*min(NXYVEC-ip1+1,MAXBUF),'+')
    end do
    if (do_xvec .and. do_yvec) then
      do ip1=1,NXYVEC,MAXBUF
        call GADGOP(XYder(:,ip1,2),MXCI*min(NXYVEC-ip1+1,MAXBUF),'+')
      end do
    end if
  end if
# endif

end subroutine NEVPT2_E4_derivative3

!-----------------------------------------------------------------------

subroutine NEVPT2_E4_XYder1(iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,CI,XYder,Gder,Zder)

  integer(kind=iwp), intent(in) :: iSym, NLEV, idx2ij(2,NLEV**2), ij2idx(NLEV,NLEV), ipxysta, ipxyend
  real(kind=wp), intent(inout) :: BUFT(:), Gder(:,:,:,:), Zder(:,:)
  real(kind=wp), intent(in) :: CI(:), XYder(:,:,:)
  integer(kind=iwp) :: ID, ip1, ip2, ipxy, issg1, issg2, istu, isvx, it, itlev, iu, iulev, iv, ivlev, ix, ixlev, iy, iylev, locx, &
                       locy, nlev2, nsgm1, nsgm2, nTask, nxy
  real(kind=wp), allocatable :: BUF1(:), BUF2(:)
  real(kind=wp), external :: ddot_

  locx = 1
  locy = 1
  if (do_xvec .and. do_yvec) then
    locy = 2
  else if (do_xvec .or. do_yvec) then
    locy = 1
  end if

  if (NXY_work == NLEV) then
    nxy = ipxyend-ipxysta+1
  else
    nxy = NLEV*nxyzdim
  end if
  NLEV2 = NLEV**2

  if (do_zder) then
    Zder(1:MXCI,1:nxy) = Zero
    nTask = NLEV2
    call Init_Tsk(ID,nTask)
    do while (Rsv_Tsk(ID,ip2))
      ivlev = idx2ij(1,ip2)
      ixlev = idx2ij(2,ip2)
      isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
      issg2 = Mul(isvx,stsym)
      nsgm2=CIS(istate)%ncsf(issg2)
      iv=SGS(istate)%L2ACT(ivlev)
      ix=SGS(istate)%L2ACT(ixlev)
      if ((NXY_work /= NLEV) .and. ((ivlev < ixyzsta) .or. (ivlev > ixyzend))) cycle
      do ip1=1,nlev2
        itlev = idx2ij(1,ip1)
        iulev = idx2ij(2,ip1)
        istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
        issg1 = Mul(istu,issg2)
        if (issg1 /= iSym) cycle
        nsgm1=CIS(istate)%ncsf(issg1)
        it = SGS(istate)%L2ACT(itlev)
        iu = SGS(istate)%L2ACT(iulev)
        if (((.not. do_xvec) .or. (it /= ix)) .and. ((.not. do_yvec) .or. (iu /= ix))) cycle
        BUFT(1:nsgm2) = Zero
        ! X_{v,u}^I = <Phi_I|E_{x,u}|Phi_J> Z_{v,x}^J
        if (do_xvec .and. (it == ix)) then
          if (NXY_work == NLEV) then
            ipxy = ij2idx(ivlev,iulev)
            if ((ipxy >= ipxysta) .and. (ipxy <= ipxyend)) BUFT(1:nsgm2) = BUFT(1:nsgm2)+XYder(1:nsgm2,ipxy-ipxysta+1,locx)
          else
            if ((ivlev >= ixyzsta) .and. (ivlev <= ixyzend)) &
              BUFT(1:nsgm2) = BUFT(1:nsgm2)+XYder(1:nsgm2,ivlev-ixyzsta+1+nxyzdim*(iulev-1),locx)
          end if
        end if
        ! Y_{v,t}^I = <Phi_I|E_{t,x}|Phi_J> Z_{v,x}^J
        if (do_yvec .and. (iu == ix)) then
          if (NXY_work == NLEV) then
            ipxy = ij2idx(ivlev,itlev)
            if ((ipxy >= ipxysta) .and. (ipxy <= ipxyend)) BUFT(1:nsgm2) = BUFT(1:nsgm2)+XYder(1:nsgm2,ipxy-ipxysta+1,locy)
          else
            if ((ivlev >= ixyzsta) .and. (ivlev <= ixyzend)) &
              BUFT(1:nsgm2) = BUFT(1:nsgm2)+XYder(1:nsgm2,ivlev-ixyzsta+1+nxyzdim*(itlev-1),locy)
          end if
        end if
        if (NXY_work == NLEV) then
          call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,One,issg2,BUFT(:),Zder(:,ip2))
        else
          call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,One,issg2,BUFT(:), &
                          Zder(:,ivlev-ixyzsta+1+nxyzdim*(ixlev-1)))
        end if
      end do
    end do

  else

    call mma_allocate(BUF1,MXCI,Label='BUF1')
    call mma_allocate(BUF2,MXCI,Label='BUF2')

    !! derivatives of Xvec and Yvec
    nTask = nlev2
    call Init_Tsk(ID,nTask)
    do while (Rsv_Tsk(ID,ip2))
      ivlev = idx2ij(1,ip2)
      ixlev = idx2ij(2,ip2)
      isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
      issg2 = Mul(isvx,stsym)
      nsgm2=CIS(istate)%ncsf(issg2)
      iv=SGS(istate)%L2ACT(ivlev)
      ix=SGS(istate)%L2ACT(ixlev)
      buf1(1:nsgm2) = Zero
      call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IVLEV,IXLEV,One,STSYM,CI,buf1(:))
      !! electron repulsion terms
      do ip1=ip2,nlev2
        itlev = idx2ij(1,ip1)
        iulev = idx2ij(2,ip1)
        istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
        issg1 = Mul(istu,issg2)
        if (issg1 /= iSym) cycle
        nsgm1=CIS(istate)%ncsf(issg1)
        it=SGS(istate)%L2ACT(itlev)
        iu=SGS(istate)%L2ACT(iulev)
        buf2(1:nsgm1) = Zero
        call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),ITLEV,IULEV,One,issg2,buf1(:),BUF2(:))
        do iylev=1,nlev
          iy=SGS(istate)%L2ACT(iylev)
          ! Etu Evx * (at|vx) -> Xau
          ipxy = ij2idx(iylev,iulev)
          if (do_xvec) then
            if (NXY_work == NLEV) then
              Gder(iy,it,iv,ix) = Gder(iy,it,iv,ix)+ddot_(nsgm1,XYder(:,ipxy-ipxysta+1,locx),1,buf2(1:nsgm1),1)
            else if ((iylev >= ixyzsta) .and. (iylev <= ixyzsta)) then
              Gder(iy,it,iv,ix) = Gder(iy,it,iv,ix)+ddot_(nsgm1,XYder(:,iylev-ixyzsta+1+nxyzdim*(iulev-1),locx),1,buf2(1:nsgm1),1)
            end if
          end if
          ! Etu Evx * (au|vx) -> Yat
          ipxy = ij2idx(iylev,itlev)
          if (do_yvec) then
            if (NXY_work == NLEV) then
              Gder(iy,iu,iv,ix) = Gder(iy,iu,iv,ix)+ddot_(nsgm1,XYder(:,ipxy-ipxysta+1,locy),1,buf2(1:nsgm1),1)
            else if ((iylev >= ixyzsta) .and. (iylev <= ixyzsta)) then
              Gder(iy,iu,iv,ix) = Gder(iy,iu,iv,ix)+ddot_(nsgm1,XYder(:,iylev-ixyzsta+1+nxyzdim*(itlev-1),locy),1,buf2(1:nsgm1),1)
            end if
          end if
        end do
        if (ip1 /= ip2) then
          ! EvxEtu = EtuEvx + Evu del(tx) - Etx del(vu)
          buft(1:nsgm1) = buf2(1:nsgm1)
          if ((itlev == ixlev) .and. (iSym == Mul(Mul(SGS(istate)%ism(iulev),SGS(istate)%ism(ivlev)),STSYM))) &
            call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IVLEV,IULEV,+One,STSYM,CI,BUFT(1:nsgm1))
          if ((iulev == ivlev) .and. (iSym == Mul(Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(ixlev)),STSYM))) &
            call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),ITLEV,IXLEV,-One,STSYM,CI,BUFT(1:nsgm1))
          do iylev=1,nlev
            iy=SGS(istate)%L2ACT(iylev)
            ! Evx Etu * (av|tu) -> Xat
            ipxy = ij2idx(iylev,ixlev)
            if (do_xvec) then
              if (NXY_work == NLEV) then
                Gder(iy,iv,it,iu) = Gder(iy,iv,it,iu)+ddot_(nsgm1,XYder(:,ipxy-ipxysta+1,locx),1,buft(1:nsgm1),1)
              else if (NXY_work /= NLEV .and. (iylev >= ixyzsta .and. iylev <= ixyzsta)) then
                Gder(iy,iv,it,iu) = Gder(iy,iv,it,iu)+ddot_(nsgm1,XYder(:,iylev-ixyzsta+1+nxyzdim*(ixlev-1),locx),1,buft(1:nsgm1),1)
              end if
            end if
            ! Evx Etu * (ax|tu) -> Yav
            ipxy = ij2idx(iylev,ivlev)
            if (do_yvec) then
              if (NXY_work == NLEV) then
                Gder(iy,ix,it,iu) = Gder(iy,ix,it,iu)+ddot_(nsgm1,XYder(:,ipxy-ipxysta+1,locy),1,buft(1:nsgm1),1)
              else if (NXY_work /= NLEV .and. (iylev >= ixyzsta .and. iylev <= ixyzsta)) then
                Gder(iy,ix,it,iu) = Gder(iy,ix,it,iu)+ddot_(nsgm1,XYder(:,iylev-ixyzsta+1+nxyzdim*(ivlev-1),locy),1,buft(1:nsgm1),1)
              end if
            end if
          end do
        end if
      end do
    end do
    call mma_deallocate(BUF1)
    call mma_deallocate(BUF2)
  end if

  call Free_Tsk(ID)

# ifdef _MOLCAS_MPP_
  if (is_real_par() .and. do_zder) then
    do ip1=1,NXYVEC,MAXBUF
      call GADGOP(Zder(:,ip1),MXCI*min(NXYVEC-ip1+1,MAXBUF),'+')
    end do
  end if
# endif

end subroutine NEVPT2_E4_XYder1

!-----------------------------------------------------------------------

subroutine NEVPT2_E4_XYder2(iSym,NLEV,idx2ij,ij2idx,ipxysta,ipxyend,BUFT,CI,Gact,XYvec,XYder,XYcont,XYcontder,ZVEC,Zder,CLag,Gder, &
                            XYtmp)

  integer(kind=iwp), intent(in) :: iSym, NLEV, idx2ij(2,NLEV**2), ij2idx(NLEV,NLEV), ipxysta, ipxyend
  real(kind=wp), intent(inout) :: BUFT(:), ZVEC(:,:), CLag(:), Gder(:,:,:,:), XYtmp(:,:,:)
  real(kind=wp), intent(in) :: CI(:), Gact(:,:,:,:), XYvec(:,:,:), XYder(:,:,:), XYcont(:,:,:), XYcontder(:,:,:), Zder(:,:)
  integer(kind=iwp) :: ID, ip1, ip1end, ip1sta, ip2, ip23, ip2_rev, ip3, ip3_rev, ipxy, issg1, issg2, istu, isvx, it, itlev, iu, &
                       iulev, iv, ivlev, ix, ixlev, izlev, locx, locy, nbufxy, nlev2, nsgm1, nsgm2, nTask, nxy
  real(kind=wp), allocatable :: Gact_sort(:), Gder_sort(:)

  locx = 1
  locy = 1
  if (do_xvec .and. do_yvec) then
    locy = 2
  else if (do_xvec .or. do_yvec) then
    locy = 1
  end if

  if (NXY_work == NLEV) then
    nxy = ipxyend-ipxysta+1
  else
    nxy = NLEV*nxyzdim
  end if
  NLEV2 = NLEV**2
  call mma_allocate(Gact_sort,NLEV2,LABEL='Gact_sort')
  call mma_allocate(Gder_sort,NLEV2,LABEL='Gder_sort')

  !! complete the CI derivative
  issg2 = 1
  if (do_zder) then

    nTask = NLEV2
    call Init_Tsk(ID,nTask)

    !ibuf = 0
    !Scal = Zero
    ! Z_{t,u}^I = \sum_{v,x} (tu|vx) <Phi_I|E_{v,x}|Psi>
    do while (Rsv_Tsk(ID,ip2))
      ivlev = idx2ij(1,ip2)
      ixlev = idx2ij(2,ip2)
      isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
      issg2 = Mul(isvx,stsym)
      nsgm2=CIS(istate)%ncsf(issg2)
      iv=SGS(istate)%L2ACT(ivlev)
      ix=SGS(istate)%L2ACT(ixlev)
      !ibuf = ibuf + 1
      BUFT(1:MXCI) = Zero
      call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IVLEV,IXLEV,One,STSYM,CI,BUFT(:))
      if (NXY_work == NLEV) then
        !call dgemv_('T',nsgm1,nxy,One,XYvec(:,:,1),mxci,BUFT(:),1,Zero,XYcont(:,1,ip1),1)
        call dgemv_('T',nsgm2,nxy,One,Zder,MXCI,BUFT,1,Zero,Gder_sort,1)
        BUFT(1:nsgm2) = Zero
        do ip1=1,nlev2
          itlev = idx2ij(1,ip1)
          iulev = idx2ij(2,ip1)
          it=SGS(istate)%L2ACT(itlev)
          iu=SGS(istate)%L2ACT(iulev)
          !Gact_sort(ip1,ibuf) = Gact(it,iu,iv,ix)
          Gder(it,iu,iv,ix) = Gder(it,iu,iv,ix)+Gder_sort(ip1)
          !call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IvLEV,IxLEV,Gact(itlev,iulev,ivlev,ixlev),STSYM,Zder(:,ip1),CLag(:))
          Gact_sort(ip1) = Gact(itlev,iulev,ivlev,ixlev)
          !BUFT(1:nsgm2) = BUFT(1:nsgm2) + Gact(itlev,iulev,ivlev,ixlev)*Zder(1:nsgm2,ip1)
        end do
        call dgemv_('N',nsgm2,nlev2,One,Zder(:,:),MXCI,Gact_sort(:),1,Zero,BUFT,1)
        !call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IvLEV,IxLEV,One,STSYM,BUFT(:),CLag(:))
      else
        call dgemv_('T',nsgm2,nxy,One,Zder,MXCI,BUFT,1,Zero,Gder_sort,1)
        do itlev=ixyzsta,ixyzend
          do iulev=1,nlev
            !ip1 = ij2idx(ixy_local,iulev_local)
            !itlev = idx2ij(1,ip1)
            !iulev = idx2ij(2,ip1)
            it=SGS(istate)%L2ACT(itlev)
            iu=SGS(istate)%L2ACT(iulev)
            !Gact_sort(iulev_local) = Gact(it,iu,iv,ix)
            !Gact_sort(itlev-ixyzsta+1+nxyzdim*(iulev-1),ibuf) = Gact(itlev,iulev,iv,ix)
            Gder(it,iu,iv,ix) = Gder(it,iu,iv,ix)+Gder_sort(itlev-ixyzsta+1+nxyzdim*(iulev-1))
            Gact_sort(itlev-ixyzsta+1+nxyzdim*(iulev-1)) = Gact(itlev,iulev,ivlev,ixlev)
          end do
        end do
        BUFT(1:nsgm2) = Zero
        call dgemv_('N',nsgm2,nxy,One,Zder(:,:),MXCI,Gact_sort(:),1,Zero,BUFT,1)
      end if

      ip2_rev = ij2idx(ixlev,ivlev)
      Zbufloop: do ip1sta=ip2_rev,nlev2,nxy
        nbufxy = 0
        ip1end = min(nlev2,ip1sta-1+nxy)
        !if (NXY_work /= NLEV) XYtmp(1:nxy,1:ip1end-ip1sta+1,1:2) = Zero
        !! XYtmp = XYtmp2 outside: XYtmp2(NXYVEC,NLEV2,2)
        if (NXY_work /= NLEV) XYtmp(:,:,:) = Zero
        do ip1=ip1sta,ip1end
          itlev = idx2ij(1,ip1)
          iulev = idx2ij(2,ip1)
          istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
          issg1 = Mul(istu,issg2)
          if (issg1 /= iSym) cycle
          nsgm1=CIS(istate)%ncsf(issg1)
          it=SGS(istate)%L2ACT(itlev)
          iu=SGS(istate)%L2ACT(iulev)
          !if (NXY_work == NLEV) then
          nbufxy = nbufxy+1
          ip23 = iTri(ip1,ip2_rev)
          XYtmp(1:nxy,nbufxy,1) = XYcontder(1:nxy,1,ip23)
          XYtmp(1:nxy,nbufxy,2) = XYcontder(1:nxy,2,ip23)
          !else if ((itlev >= ixyzsta) .and. (iulev <= ixyzend)) then
          !  nbufxy = nbufxy+1
          !  ip23 = iTri(ip1,ip2_rev)
          !  !XYtmp(1:nxy,nbufxy,1) = XYcontder(1:nxy,1,ip23)
          !  !XYtmp(1:nxy,nbufxy,2) = XYcontder(1:nxy,2,ip23)
          !  XYtmp(1:nxy,itlev-ixyzsta+1+nxyzdim*(iulev-1),1) = XYcontder(1:nxy,1,ip23)
          !  XYtmp(1:nxy,itlev-ixyzsta+1+nxyzdim*(iulev-1),2) = XYcontder(1:nxy,2,ip23)
          !  nbufxy = nxy
          !end if
        end do

        ZVEC(1:nsgm1,1:nbufxy) = Zero
        if (do_xvec .and. do_yvec) then
          call dgemm_('N','N',nsgm1,nbufxy,nxy,One,XYvec(:,:,1),MXCI,XYtmp(:,:,1),NXYVEC,One,ZVEC(:,:),MXCI)
          call dgemm_('N','N',nsgm1,nbufxy,nxy,One,XYvec(:,:,2),MXCI,XYtmp(:,:,2),NXYVEC,One,ZVEC(:,:),MXCI)
        else if (do_xvec .and. (.not. do_yvec)) then
          call dgemm_('N','N',nsgm1,nbufxy,nxy,One,XYvec(:,:,1),MXCI,XYtmp(:,:,1),NXYVEC,One,ZVEC(:,:),MXCI)
        else if ((.not. do_xvec) .and. do_yvec) then
          call dgemm_('N','N',nsgm1,nbufxy,nxy,One,XYvec(:,:,1),MXCI,XYtmp(:,:,2),NXYVEC,One,ZVEC(:,:),MXCI)
        end if

        nbufxy = 0
        do ip1=ip1sta,ip1end
          itlev = idx2ij(1,ip1)
          iulev = idx2ij(2,ip1)
          istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
          issg1 = Mul(istu,issg2)
          if (issg1 /= iSym) cycle
          nsgm1=CIS(istate)%ncsf(issg1)
          it=SGS(istate)%L2ACT(itlev)
          iu=SGS(istate)%L2ACT(iulev)
          !if (NXY_work == NLEV) then
          nbufxy = nbufxy+1
          call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,One,issg1,ZVEC(:,nbufxy),BUFT(:))
          !else if ((itlev >= ixyzsta) .and. (iulev <= ixyzend)) then
          !  nbufxy = nbufxy + 1
          !  !nbufxy = itlev-ixyzsta+1+nxyzdim*(iulev-1)
          !  call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,One,issg1,ZVEC(:,nbufxy),BUFT(:))
          !end if
        end do
      end do Zbufloop

      if (do_xvec .and. do_yvec) then
        call dgemv_('N',nsgm1,nxy,One,XYvec(:,:,1),mxci,XYcont(:,1,ip2),1,One,BUFT(:),1)
        call dgemv_('N',nsgm1,nxy,One,XYvec(:,:,2),mxci,XYcont(:,2,ip2),1,One,BUFT(:),1)
      else if (do_xvec .and. (.not. do_yvec)) then
        call dgemv_('N',nsgm1,nxy,One,XYvec(:,:,1),mxci,XYcont(:,1,ip2),1,One,BUFT(:),1)
      else if ((.not. do_xvec) .and. do_yvec) then
        call dgemv_('N',nsgm1,nxy,One,XYvec(:,:,1),mxci,XYcont(:,2,ip2),1,One,BUFT(:),1)
      end if

      call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IVLEV,IXLEV,One,issg2,BUFT(:),CLag)
    end do

  else

    nTask = nlev2
    call Init_Tsk(ID,nTask)
    do while (Rsv_Tsk(ID,ip3))
      ivlev = idx2ij(1,ip3)
      ixlev = idx2ij(2,ip3)
      isvx=Mul(SGS(istate)%ism(ivlev),SGS(istate)%ism(ixlev))
      issg1 = Mul(isvx,issg2)
      if (issg1 /= iSym) cycle
      nsgm1=CIS(istate)%ncsf(issg1)
      iv=SGS(istate)%L2ACT(ivlev)
      ix=SGS(istate)%L2ACT(ixlev)
      BUFT(1:mxci) = Zero
      do ip2=1,nlev2
        itlev = idx2ij(1,ip2)
        iulev = idx2ij(2,ip2)
        istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
        issg1 = Mul(istu,issg2)
        if (issg1 /= iSym) cycle
        nsgm1=CIS(istate)%ncsf(issg1)
        it=SGS(istate)%L2ACT(itlev)
        iu=SGS(istate)%L2ACT(iulev)
        do izlev=1,NLEV
          !! X(z,u) = (zt|vx)*<I|EtuEvx|Psi>
          ipxy = ij2idx(izlev,iulev)
          if (do_xvec .and. (ipxy >= ipxysta) .and. (ipxy <= ipxyend)) &
            call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,Gact(izlev,itlev,ivlev,ixlev),issg2, &
                          XYder(:,ipxy-ipxysta+1,locx),BUFT(:))
          !! Y(z,t) = (zu|vx)*<I|EtuEvx|Psi>
          ipxy = ij2idx(izlev,itlev)
          if (do_yvec .and. (ipxy >= ipxysta) .and. (ipxy <= ipxyend)) &
            call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,Gact(izlev,iulev,ivlev,ixlev),issg2, &
                          XYder(:,ipxy-ipxysta+1,locy),BUFT(:))
        end do
      end do

      ip3_rev = ij2idx(ixlev,ivlev)
      do ip2=ip3_rev,nlev2
        itlev = idx2ij(1,ip2)
        iulev = idx2ij(2,ip2)
        istu=Mul(SGS(istate)%ism(itlev),SGS(istate)%ism(iulev))
        issg1 = Mul(istu,issg2)
        if (issg1 /= iSym) cycle
        nsgm1=CIS(istate)%ncsf(issg1)
        it = SGS(istate)%L2ACT(itlev)
        iu = SGS(istate)%L2ACT(iulev)
        ip23 = iTri(ip2,ip3_rev)
        do ipxy=1,nxy
          if (do_xvec .and. do_yvec) then
            call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,XYcontder(ipxy,1,ip23),issg1,XYvec(:,ipxy,1),BUFT(:))
            call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,XYcontder(ipxy,2,ip23),issg1,XYvec(:,ipxy,2),BUFT(:))
          else if (do_xvec .and. (.not. do_yvec)) then
            call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,XYcontder(ipxy,1,ip23),issg1,XYvec(:,ipxy,1),BUFT(:))
          else if ((.not. do_xvec) .and. do_yvec) then
            call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IULEV,ITLEV,XYcontder(ipxy,2,ip23),issg1,XYvec(:,ipxy,1),BUFT(:))
          end if
        end do
      end do

      if (do_xvec .and. do_yvec) then
        call dgemv_('N',nsgm1,nxy,One,XYvec(:,:,1),mxci,XYcont(:,1,ip3),1,One,BUFT(:),1)
        call dgemv_('N',nsgm1,nxy,One,XYvec(:,:,2),mxci,XYcont(:,2,ip3),1,One,BUFT(:),1)
      else if (do_xvec .and. (.not. do_yvec)) then
        call dgemv_('N',nsgm1,nxy,One,XYvec(:,:,1),mxci,XYcont(:,1,ip3),1,One,BUFT(:),1)
      else if ((.not. do_xvec) .and. do_yvec) then
        call dgemv_('N',nsgm1,nxy,One,XYvec(:,:,1),mxci,XYcont(:,2,ip3),1,One,BUFT(:),1)
      end if

      call SG_Epq_Psi(SGS(istate),CIS(istate),EXS(istate),IVLEV,IXLEV,One,issg2,BUFT(:),CLag)
    end do
  end if

  call Free_Tsk(ID)

  call mma_deallocate(Gact_sort)
  call mma_deallocate(Gder_sort)

end subroutine NEVPT2_E4_XYder2

end module NEVPT2_E4
