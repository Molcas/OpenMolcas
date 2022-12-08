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

subroutine David5(nDet,mxItr,nItr,CI_Conv,ThrEne,iSel,ExplE,ExplV,HTUTRI,GTUVXTRI)

use citrans, only: citrans_csf2sd, citrans_sd2csf, citrans_sort

use csfbas, only: CONF, CTS, KCFTP, KDTOC
use faroald, only: my_norb, ndeta, ndetb, sigma_update
use davctl_mod, only: istart, n_Roots, nkeep, nvec
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
#include "rasdim.fh"
#include "rasrc.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "timers.fh"
#include "rasscf_lucia.fh"
#include "output_ras.fh"
! lroots, maxjt in rasscf.fh
! nsel in general.fh
integer(kind=iwp), intent(in) :: nDet, iSel(nSel)
integer(kind=iwp), intent(inout) :: mxItr
integer(kind=iwp), intent(out) :: nItr
real(kind=wp), intent(out) :: CI_Conv(2,lRoots,MAXJT)
real(kind=wp), intent(in) :: ThrEne, ExplE(nSel), ExplV(nSel,nSel), HTUTRI(*), GTUVXTRI(*)
integer(kind=iwp) :: i, iConf, iConv, idelta, iDummy, ij, IPRLEV, iskipconv, it, it_ci, itu, ituvx, iu, iv, ix, ixmax, jRoot, &
                     kRoot, l1, l2, l3, lPrint, mRoot, nBasVec, nconverged, nleft, nnew, ntrial
real(kind=wp) :: Alpha(mxRoot), Beta(mxRoot), Cik, Dummy(1), E0, E1, ECORE_HEX, FP, Hji, ovl, R, RR, scl, Sji, ThrRes, updsiz, Z
logical(kind=iwp) :: Skip
integer(kind=iwp), allocatable :: vkcnf(:)
real(kind=wp), allocatable :: Cs(:), ctemp(:), Es(:), gtuvx(:,:,:,:), Hs(:), htu(:,:), psi(:,:), Scr1(:,:), Scr2(:,:), Scr3(:,:), &
                              sigtemp(:), sgm(:,:), Ss(:), Vec1(:), Vec3(:), VECSVC(:)
real(kind=wp), allocatable, target :: Tmp(:)
real(kind=wp), pointer, contiguous :: Vec2(:)
integer(kind=iwp), external :: ip_of_Work
real(kind=wp), external :: dDot_, dnrm2_, GET_ECORE

!-----------------------------------------------------------------------
! MGD dec 2017 : When optimizing many states, the lowest ones tend to
! converge much faster than the rest. Changed the code so that the converged states
! are not optimize further, saving potentially a lot of time.

if (DoFaro) then
  ! fill in the integrals from their triangular storage
  call mma_allocate(htu,my_norb,my_norb,label='htu')
  call mma_allocate(gtuvx,my_norb,my_norb,my_norb,my_norb,label='gtuvx')
  htu(:,:) = Zero
  gtuvx(:,:,:,:) = Zero
  itu = 0
  ituvx = 0
  do it=1,my_norb
    do iu=1,it
      itu = itu+1
      !write(u6,'(1x,3I4,F21.14)') it,iu,itu,htutri(itu)
      htu(iu,it) = htutri(itu)
      htu(it,iu) = htutri(itu)
      do iv=1,it
        ixmax = iv
        if (it == iv) ixmax = iu
        do ix=1,ixmax
          ituvx = ituvx+1
          !write(u6,'(1x,5I4,F21.14)') it,iu,iv,ix,ituvx,gtuvxtri(ituvx)
          GTUVX(IT,IU,IV,IX) = GTUVXTRI(ITUVX)
          GTUVX(IU,IT,IV,IX) = GTUVXTRI(ITUVX)
          GTUVX(IT,IU,IX,IV) = GTUVXTRI(ITUVX)
          GTUVX(IU,IT,IX,IV) = GTUVXTRI(ITUVX)
          GTUVX(IV,IX,IT,IU) = GTUVXTRI(ITUVX)
          GTUVX(IX,IV,IT,IU) = GTUVXTRI(ITUVX)
          GTUVX(IV,IX,IU,IT) = GTUVXTRI(ITUVX)
          GTUVX(IX,IV,IU,IT) = GTUVXTRI(ITUVX)
        end do
      end do
    end do
  end do
  ! Euhm, stuff needed for awkward conversions from a
  ! non-specified SYG to GUGA format befor converting to
  ! determinants. This is because for Lucia, CSFs have been
  ! converted to SYG format somewhere up in cistart.
  call mma_allocate(VECSVC,nconf,label='CIVEC')
  call mma_allocate(vkcnf,nactel,label='kcnf')
end if

call Timing(Alfex_1,Swatch,Swatch,Swatch)
Rc_CI = 0
IPRLEV = IPRLOC(3)

! allocate space for CI-vectors
l1 = nKeep
l2 = l1*l1
l3 = (l2+l1)/2
! Trying to avoid writing out of bounds in CSDTVC :::: JESPER :::: CHEAT
call mma_allocate(Vec1,nConf,label='Vector1')
call mma_allocate(Tmp,ndet,label='Vector2')
Vec2(1:nConf) => Tmp(1:nConf)
call mma_allocate(Vec3,nConf,label='Vector3')
call mma_allocate(Es,l1,label='Esmall')
call mma_allocate(Hs,l3,label='Hsmall')
call mma_allocate(Ss,l3,label='Ssmall')
call mma_allocate(Cs,l2,label='Csmall')
call mma_allocate(Scr1,nSel,lRoots,label='Scr1')
call mma_allocate(Scr2,nSel,lRoots,label='Scr2')
call mma_allocate(Scr3,nSel,lRoots,label='Scr3')
if (.not. DoFaro) then
  call mma_allocate(ctemp,ndet,label='CTEMP')
  call mma_allocate(sigtemp,ndet,label='SIGTEM')
end if
!-----------------------------------------------------------------------

! Print convergence thresholds in ITERFILE
ThrRes = max(0.2e-6_wp,sqrt(ThrEne))
write(IterFile,'(19X,A,F18.10)') '- Threshold for energy   ...:',ThrEne
write(IterFile,'(19X,A,F18.10)') '- Threshold for Residual ...:',ThrRes
write(IterFile,'(20A4)') ('****',i=1,20)
write(IterFile,*)
write(IterFile,'(1X,A4,4X,A4,4X,A18,4X,A14,4X,A14)') 'Iter','Root','Energy','dE','Residual'
write(IterFile,'(72A1)') ('=',i=1,72)
call xFlush(IterFile)
!=======================================================================
! start long loop over iterations
Skip = .false.
nconverged = 0
iskipconv = 1
nnew = 0
nvec = lRoots
do it_ci=1,mxItr
  !---------------------------------------------------------------------
  ! MGD for stability purposes recompute sigma vec from time to time
  idelta = 1
  if ((mod(it_ci-1,24) == 0)) idelta = 0
  do mRoot=lRoots*idelta+1,lRoots+nnew
    ! New CI vectors (it_ci,mroot) are available.
    ! compute new sigma vectors
    call Load_CI_vec(mRoot,nConf,Vec1,LuDavid)
    if (iprlev >= DEBUG) then
      lPrint = min(nConf,200)
      write(u6,'(1X,A,I2,A,I2)') 'CI vector, iter =',it_ci,' mRoot =',mRoot
      write(u6,'(1X,A)') '(max. 200 elements)'
      write(u6,'(1X,A)') '-----------------------------'
      call dVcPrt(' ',' ',Vec1,lPrint)
    end if

    call Timing(Rolex_1,Swatch,Swatch,Swatch)
    if (DOFARO) then
      ! determinant wavefunctions
      call mma_allocate(sgm,ndeta,ndetb,label='sgm')
      call mma_allocate(psi,ndeta,ndetb,label='psi')

      VECSVC(:) = Zero
      call REORD2(MY_NORB,NACTEL,1,0,CONF,IWORK(KCFTP),VEC1,VECSVC,VKCNF)
      call CITRANS_SORT('C',VECSVC,VEC2)
      PSI = Zero
      call CITRANS_CSF2SD(VEC2,PSI)
      SGM = Zero
      call SIGMA_UPDATE(HTU,GTUVX,SGM,PSI)
      call CITRANS_SD2CSF(SGM,VEC2)
      call CITRANS_SORT('O',VEC2,VECSVC)
      call Reord2(my_norb,NACTEL,1,1,CONF,iWork(KCFTP),VECSVC,VEC2,VKCNF)

      if (iprlev >= DEBUG) then
        FP = DNRM2_(NCONF,VEC2,1)
        write(u6,'(1X,A,F21.14)') 'sigma dnrm2_(faroald): ',FP
      end if

      ! free the arrays
      call mma_deallocate(sgm)
      call mma_deallocate(psi)
    else
      ! Convert the CI-vector from CSF to Det. basis
      ctemp(1:nConf) = Vec1(:)
      sigtemp(:) = Zero
      call csdtvc(ctemp,sigtemp,1,work(kdtoc),cts,stSym,1)
      c_pointer = ip_of_Work(ctemp(1))
      ! Calling Lucia to determine the sigma vector
      call Lucia_Util('Sigma',iDummy,iDummy,Dummy)
      ! Set mark so densi_master knows that the Sigma-vector exists on disk.
      iSigma_on_disk = 1
      call CSDTVC(Tmp,ctemp,2,work(kdtoc),cts,stSym,1)

      if (iprlev >= DEBUG) then
        FP = DNRM2_(NCONF,VEC2,1)
        write(u6,'(1X,A,F21.14)') 'sigma dnrm2_(lucia):   ',FP
      end if
    end if

    ! Add ECORE_HEX (different from zero when particle-hole formalism used)
    ECORE_HEX = GET_ECORE()
    Vec1(:) = Vec1(:)+ecore_hex*Vec2(:)
    ! Timings on generation of the sigma vector
    call Timing(Rolex_2,Swatch,Swatch,Swatch)
    Rolex_2 = Rolex_2-Rolex_1
    Rolex_3 = Rolex_3+Rolex_2

    if (iprlev >= DEBUG) then
      lPrint = min(nConf,200)
      write(u6,*)
      write(u6,'(1X,A,I2,A,I2)') 'sigma vector, iter =',it_ci,' mRoot =',mRoot
      write(u6,'(1X,A)') '(max. 200 elements)'
      write(u6,'(1X,A)') '--------------------------------'
      call dVcPrt(' ',' ',Vec2,lPrint)
    end if
    call Save_Sig_vec(mRoot,nConf,Vec2,LuDavid)
  end do
  ! Sigma vectors (it_ci,mroot) have been computed, for mroot=1..lroots
  !---------------------------------------------------------------------
  ! compute Hsmall and Ssmall
  ! These are Hsmall(jtrial,ktrial), where jtrial is (jter,jroot), and
  ! ktrial is (kter,kroot), and similar Ssmall.
  ! jtrial=1..mxKeep*lroots correspond to jter=it_ci-mxKeep+1..it_ci
  ! (Fewer, at the beginning)

  do jRoot=1,nvec
    call Load_CI_vec(jRoot,nConf,Vec1,LuDavid)
    call Load_Sig_vec(jRoot,nConf,Vec2,LuDavid)
    do kRoot=1,jRoot
      call Load_CI_vec(kRoot,nConf,Vec3,LuDavid)
      ij = kRoot+(jRoot*jRoot-jRoot)/2
      Sji = dDot_(nConf,Vec1,1,Vec3,1)
      Hji = dDot_(nConf,Vec2,1,Vec3,1)
      Ss(ij) = Sji
      Hs(ij) = Hji
    end do
  end do
  ntrial = nvec
  if (iprlev >= DEBUG) then
    call TriPrt('Hsmall',' ',Hs,ntrial)
    call TriPrt('Ssmall',' ',Ss,ntrial)
  end if
  ! Hsmall and Ssmall have been computed (ntrial x ntrial, in triangular
  ! storage.)
  !---------------------------------------------------------------------

  ! solve secular equation HC=SCE.

  ! PAM2009 nBasVec on input = min(ntrial,nconf)
  ! nBasVec returned as nr of orthonormal solutions to HC=SCE
  nBasVec = nConf
  call HCSCE(ntrial,Hs,Ss,Cs,Es,nBasVec)
  if (nBasVec < lRoots) then
    write(u6,*) 'David: nBasVec less than lRoots'
    write(u6,*) 'nBasvec, lRoots = ',nBasVec,lRoots
    if (ICIRST == 1) write(u6,*) 'CIREstart was used. Check the number of roots in the previous calculation'
    call Abend()
  end if
  if (iprlev >= DEBUG) then
    call dVcPrt('Eigenvalues of Hsmall',' ',Es,ntrial)
    call RecPrt('Eigenvectors of Hsmall',' ',Cs,ntrial,ntrial)
  end if
  !---------------------------------------------------------------------
  ! compute the current 'best' CI, sigma and residual vector

  ! CI vector is Vec1
  ! sigma vector is saved in Vec2
  ! residual vector is saved in Vec3
  do mRoot=1,lRoots
    ! initialize 'best' CI and sigma vector
    Vec1(:) = Zero
    Vec2(:) = Zero
    ! accumulate contributions
    do jRoot=1,nvec
      Cik = Cs(jRoot+(mRoot-1)*ntrial)
      call Load_CI_vec(jRoot,nConf,Vec3,LuDavid)
      Vec1(:) = Vec1(:)+Cik*Vec3(:)
      call Load_Sig_vec(jRoot,nConf,Vec3,LuDavid)
      Vec2(:) = Vec2(:)+Cik*Vec3(:)
    end do
    RR = dDot_(nConf,Vec1,1,Vec1,1)
    scl = One/sqrt(RR)
    Vec1(:) = scl*Vec1(:)
    Vec2(:) = scl*Vec2(:)
    call Save_tmp_CI_vec(mRoot,nConf,Vec1,LuDavid)
    call Save_tmp_Sig_vec(mRoot,nConf,Vec2,LuDavid)
    ! compute residual vector
    E0 = Es(mRoot)
    Vec3(:) = Vec2(:)-E0*Vec1(:)
    ! save current best energy and residual
    RR = dDot_(nConf,Vec3,1,Vec3,1)
    CI_conv(1,mroot,it_ci) = E0
    CI_conv(2,mroot,it_ci) = sqrt(RR)
    ! print vectors
    if (iprlev >= DEBUG) then
      lPrint = min(nConf,200)
      write(u6,'(1X,A,I2,A,I2)') 'new best CI vector, iter =',it_ci,' mRoot =',mRoot
      write(u6,'(1X,A)') '(max. 200 elements)'
      write(u6,'(1X,A)') '--------------------------------------'
      call dVcPrt(' ',' ',Vec1,lPrint)
      write(u6,'(1X,A,I2,A,I2)') 'new best sigma vector, iter =',it_ci,' mRoot =',mRoot
      write(u6,'(1X,A)') '(max. 200 elements)'
      write(u6,'(1X,A)') '-----------------------------------------'
      call dVcPrt(' ',' ',Vec2,lPrint)
      write(u6,'(1X,A,I2,A,I2)') 'new residual vector, iter =',it_ci,' mRoot =',mRoot
      write(u6,'(1X,A)') '(max. 200 elements)'
      write(u6,'(1X,A)') '-----------------------------------------'
      call dVcPrt(' ',' ',Vec3,lPrint)
    end if
    ! to improve the preconditioner select all elements in the
    ! subspace of the explicit Hamiltonian
    if (nSel > 1) then
      do i=1,nSel
        iConf = iSel(i)
        Scr1(i,mRoot) = Vec3(iConf)
        Scr2(i,mRoot) = Vec1(iConf)
      end do
    end if
  end do
  ! Current best CI & Sigma vectors have been stored in a temporary place
  ! for mroot=1..lroots.
  ! Also, the selected elements of the CI and Sigma vectors have been
  ! saved at Scr1 (Sigma) and Scr2 (CI)
  !---------------------------------------------------------------------

  ! check for convergence
  nItr = it_ci
  if (it_ci > 1) then
    dE = CI_conv(1,1,it_ci-1)-CI_conv(1,1,it_ci)
  else
    dE = Zero
  end if
  write(IterFile,'(1X,I4,4X,I4,4X,F18.10,4X,F14.10,4X,F14.10)') it_ci,1,CI_conv(1,1,it_ci),dE,CI_conv(2,1,it_ci)
  do jRoot=2,lRoots
    if (it_ci > 1) then
      dE = CI_conv(1,jroot,it_ci-1)-CI_conv(1,jroot,it_ci)
    end if
    write(IterFile,'(9X,I4,4X,F18.10,4X,F14.10,4X,F14.10)') jRoot,CI_conv(1,jroot,it_ci),dE,CI_conv(2,jroot,it_ci)
  end do
  if (lRoots > 1) write(IterFile,*)
  call xFlush(IterFile)
  if (iprlev > DEBUG) then
    write(u6,*)
    write(u6,'(1X,120A1)') ('*',i=1,120)
    write(u6,'(1X,A,I2)') 'CI iteration ',it_ci
    ThrRes = max(0.2e-6_wp,sqrt(ThrEne))
    write(u6,'(1X,A,2F18.10)') 'ThrEne,ThrRes=',ThrEne,ThrRes
    do jRoot=1,lRoots
      if (it_ci > 1) then
        dE = CI_conv(1,jroot,it_ci-1)-CI_conv(1,jroot,it_ci)
      else
        dE = Zero
      end if
      write(u6,'(1X,A,I2,A,F18.10,2(A,F14.10))') ' root ',jRoot,' energy =',CI_conv(1,jroot,it_ci),' dE =',dE,' residual =', &
                                                 CI_conv(2,jroot,it_ci)
    end do
    write(u6,'(1X,120A1)') ('*',i=1,120)
    write(u6,*)
  end if
  ThrRes = max(0.2e-6_wp,sqrt(ThrEne))
  iConv = 0
  nconverged = 0
  ! Do not check for convergence of hidden roots
  do jRoot=1,lRoots-hroots
    if (it_ci > 1) then
      dE = CI_conv(1,jroot,it_ci-1)-CI_conv(1,jroot,it_ci)
    else
      dE = Zero
    end if
    dE = abs(dE)
    R = CI_conv(2,jroot,it_ci)
    if ((dE < ThrEne) .and. (R < ThrRes)) then
      iConv = iConv+1
      if (jRoot == nconverged+1) nconverged = nconverged+1
    end if
  end do
  if (iskipconv == 0) nconverged = 0
  if (iConv >= lRoots-hroots) then
    Skip = .true.
    exit
  end if
  !---------------------------------------------------------------------
  ! compute correction vectors q1 = r/(E0-H) and q2 = c/(E0-H)

  nleft = lRoots-nconverged
  if (nSel > 1) then
    call DGEMM_('T','N',nSel,nleft,nSel,One,ExplV,nSel,Scr1(:,nconverged+1),nSel,Zero,Scr3,nSel)
    do mRoot=nconverged+1,lRoots
      E0 = Es(mRoot)
      do i=1,nSel
        Z = E0-ExplE(i)
        if (abs(Z) < 0.001_wp) Z = 0.001_wp
        Scr3(i,mRoot-nconverged) = Scr3(i,mRoot-nconverged)/Z
      end do
    end do
    call DGEMM_('N','N',nSel,nleft,nSel,One,ExplV,nSel,Scr3,nSel,Zero,Scr1,nSel)
    call DGEMM_('T','N',nSel,nleft,nSel,One,ExplV,nSel,Scr2(:,nconverged+1),nSel,Zero,Scr3,nSel)
    do mRoot=nconverged+1,lRoots
      E0 = Es(mRoot)
      do i=1,nSel
        Z = E0-ExplE(i)
        if (abs(Z) < 0.001_wp) Z = 0.001_wp
        Scr3(i,mRoot-nconverged) = Scr3(i,mRoot-nconverged)/Z
      end do
    end do
    call DGEMM_('N','N',nSel,nleft,nSel,One,ExplV,nSel,Scr3,nSel,Zero,Scr2,nSel)
  end if
  !---------------------------------------------------------------------
  do mRoot=nconverged+1,lRoots
    E0 = -Es(mRoot)
    call Load_tmp_Sig_vec(mRoot,nConf,Vec1,LuDavid)
    call Load_tmp_CI_vec(mRoot,nConf,Vec2,LuDavid)
    Vec1(:) = Vec1(:)+E0*Vec2(:)
    call Load_H_diag(nConf,Vec3,LuDavid)
    E0 = Es(mRoot)
    do i=1,nConf
      Z = E0-Vec3(i)
      if (abs(Z) < 1.0e-4_wp) Z = 1.0e-4_wp
      Vec3(i) = Vec1(i)/Z
    end do
    if (nSel > 1) then
      do i=1,nSel
        iConf = iSel(i)
        Vec3(iConf) = Scr1(i,mRoot-nconverged)
      end do
    end if
    Alpha(mRoot) = dDot_(nConf,Vec3,1,Vec2,1)
    call Load_H_diag(nConf,Vec3,LuDavid)
    E0 = Es(mRoot)
    do i=1,nConf
      Z = E0-Vec3(i)
      if (abs(Z) < 1.0e-4_wp) Z = 1.0e-4_wp
      Vec3(i) = Vec2(i)/Z
    end do
    if (nSel > 1) then
      do i=1,nSel
        iConf = iSel(i)
        Vec3(iConf) = Scr2(i,mRoot-nconverged)
      end do
    end if
    Beta(mRoot) = dDot_(nConf,Vec3,1,Vec2,1)
  end do
  !---------------------------------------------------------------------

  ! compute correction vectors q3 = (r-E1*q2)/(E0-H)
  if (nSel > 1) then
    do mRoot=nconverged+1,lRoots
      call Load_tmp_Sig_vec(mRoot,nConf,Vec1,LuDavid)
      call Load_tmp_CI_vec(mRoot,nConf,Vec2,LuDavid)
      E0 = -Es(mRoot)
      E1 = -Alpha(mRoot)/Beta(mRoot)
      Vec1(:) = Vec1(:)+(E0+E1)*Vec2(:)
      do i=1,nSel
        iConf = iSel(i)
        Scr1(i,mRoot-nconverged) = Vec1(iConf)
      end do
    end do
    call DGEMM_('T','N',nSel,nleft,nSel,One,ExplV,nSel,Scr1,nSel,Zero,Scr3,nSel)
    do mRoot=nconverged+1,lRoots
      E0 = Es(mRoot)
      do i=1,nSel
        Z = E0-ExplE(i)
        if (abs(Z) < 0.001_wp) Z = 0.001_wp
        Scr3(i,mRoot-nconverged) = Scr3(i,mRoot-nconverged)/Z
      end do
    end do
    call DGEMM_('N','N',nSel,nleft,nSel,One,ExplV,nSel,Scr3,nSel,Zero,Scr1,nSel)
  end if
  ! move the index of CI_vec
  istart = istart+nnew
  istart = mod(istart,nkeep-n_Roots)

  nnew = 0
  do mRoot=nconverged+1,lRoots
    call Load_tmp_Sig_vec(mRoot,nConf,Vec1,LuDavid)
    call Load_tmp_CI_vec(mRoot,nConf,Vec2,LuDavid)
    E0 = -Es(mRoot)
    E1 = -Alpha(mRoot)/Beta(mRoot)
    Vec1(:) = Vec1(:)+(E0+E1)*Vec2(:)
    call Load_H_diag(nConf,Vec3,LuDavid)
    E0 = Es(mRoot)
    do i=1,nConf
      Z = E0-Vec3(i)
      if (abs(Z) < 1.0e-4_wp) Z = 1.0e-4_wp
      Vec3(i) = Vec1(i)/Z
    end do
    if (nSel > 1) then
      do i=1,nSel
        iConf = iSel(i)
        Vec3(iConf) = Scr1(i,mRoot-nconverged)
      end do
    end if
    ! Orthonormalize wrt previous vectors
    updsiz = dnrm2_(nconf,Vec3,1)
    Vec3(:) = Vec3(:)/updsiz
    do jRoot=lRoots+1,min(nvec,nkeep-nconverged)
      call Load_CI_vec(jRoot,nConf,Vec2,LuDavid)
      ovl = dDot_(nConf,Vec3,1,Vec2,1)
      Vec3(:) = Vec3(:)-ovl*Vec2(:)
    end do
    updsiz = dnrm2_(nconf,Vec3,1)
    if (updsiz > 1.0e-6_wp) then
      Vec3(:) = Vec3(:)/updsiz
      nnew = nnew+1
      nvec = nvec+1
      nvec = min(nvec,nkeep)
      call Save_CI_vec(lRoots+mRoot-nconverged,nConf,Vec3,LuDavid)
    end if
  end do
  !---------------------------------------------------------------------
  ! move the current best CI and sigma vectors to the first place
  ! in the list of retained CI vectors
  do mRoot=1,lRoots
    call Load_tmp_CI_vec(mRoot,nConf,Vec1,LuDavid)
    call Save_CI_vec(mRoot,nConf,Vec1,LuDavid)
    call Load_tmp_Sig_vec(mRoot,nConf,Vec1,LuDavid)
    call Save_Sig_vec(mRoot,nConf,Vec1,LuDavid)
  end do

! end of the long loop over iterations
end do
!=======================================================================

if (.not. Skip) then
  mxItr = min(mxCiIt,mxItr+12)
  if (IPRLEV >= USUAL) then
    write(u6,*) '       No convergence in the CI section: MAXJT will be increased to ',mxItr
  end if
  Rc_CI = 16
  nItr = nItr-1
end if

! deallocate local temporary vectors
call mma_deallocate(Vec1)
nullify(Vec2)
call mma_deallocate(Tmp)
call mma_deallocate(Vec3)
call mma_deallocate(Es)
call mma_deallocate(Hs)
call mma_deallocate(Ss)
call mma_deallocate(Cs)
call mma_deallocate(Scr1)
call mma_deallocate(Scr2)
call mma_deallocate(Scr3)
if (DoFaro) then
  call mma_deallocate(htu)
  call mma_deallocate(gtuvx)
  call mma_deallocate(VECSVC)
  call mma_deallocate(vkcnf)
else
  call mma_deallocate(ctemp)
  call mma_deallocate(sigtemp)
end if

call Timing(Alfex_2,Swatch,Swatch,Swatch)
Alfex_2 = Alfex_2-Alfex_1
Alfex_3 = Alfex_3+Alfex_2

return

end subroutine David5
