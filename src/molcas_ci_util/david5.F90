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

use citrans
use faroald
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

implicit none
#include "rasdim.fh"
#include "rasrc.fh"
#include "rasscf.fh"
#include "general.fh"
#include "csfbas.fh"
#include "davctl.fh"
#include "WrkSpc.fh"
#include "timers.fh"
#include "rasscf_lucia.fh"
#include "output_ras.fh"
! lroots, maxjt in rasscf.fh
! nsel in general.fh
integer(kind=iwp) :: nDet, mxItr, nItr, iSel(nSel)
real(kind=wp) :: CI_Conv(2,lRoots,MAXJT), ThrEne, ExplE(nSel), ExplV(nSel,nSel), HTUTRI(*), GTUVXTRI(*)
integer(kind=iwp) :: i, iConf, iConv, iCs, idelta, iDummy, iEs, iHs, ij, iOff, IPRLEV, iScr1, iScr2, iScr3, iScr4, iScr5, &
                     iskipconv, iSs, it, it_ci, itu, ituvx, iu, iv, iVec1, iVec2, iVec3, IVECSVC, ivkcnf, ix, ixmax, jRoot, &
                     kctemp, kRoot, ksigtemp, l1, l2, l3, lPrint, mRoot, nBasVec, nconverged, nleft, nnew, ntrial
real(kind=wp) :: Alpha(mxRoot), Beta(mxRoot), Cik, Dummy(1), E0, E1, ECORE_HEX, FP, Hji, ovl, R, RR, scl, Sji, ThrRes, updsiz, Z
real(kind=wp), allocatable :: gtuvx(:,:,:,:), htu(:,:), sgm(:,:), psi(:,:)
real(kind=r8), external :: dDot_, dnrm2_, GET_ECORE

!-----------------------------------------------------------------------
! MGD dec 2017 : When optimizing many states, the lowest ones tend to
! converge much faster than the rest. Changed the code so that the converged states
! are not optimize further, saving potentially a lot of time.

if (DoFaro) then
  ! fill in the integrals from their triangular storage
  allocate(htu(my_norb,my_norb))
  allocate(gtuvx(my_norb,my_norb,my_norb,my_norb))
  htu = Zero
  gtuvx = Zero
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
  call GetMem('CIVEC','Allo','Real',IVECSVC,nconf)
  call getmem('kcnf','allo','inte',ivkcnf,nactel)
end if

call Timing(Alfex_1,Swatch,Swatch,Swatch)
Rc_CI = 0
IPRLEV = IPRLOC(3)

! allocate space for CI-vectors
l1 = nKeep
l2 = l1*l1
l3 = (l2+l1)/2
! Trying to avoid writing out of bounds in CSDTVC :::: JESPER :::: CHEAT
call GetMem('Vector1','Allo','Real',iVec1,ndet)
call GetMem('Vector2','Allo','Real',iVec2,ndet)
call GetMem('Vector3','Allo','Real',iVec3,ndet)
call GetMem('Esmall','Allo','Real',iEs,l1)
call GetMem('Hsmall','Allo','Real',iHs,l3)
call GetMem('Ssmall','Allo','Real',iSs,l3)
call GetMem('Csmall','Allo','Real',iCs,l2)
call GetMem('Scr1','Allo','Real',iScr1,l2)
call GetMem('Scr2','Allo','Real',iScr2,l2)
call GetMem('Scr3','Allo','Real',iScr3,lRoots*nSel)
call GetMem('Scr4','Allo','Real',iScr4,lRoots*nSel)
call GetMem('Scr5','Allo','Real',iScr5,lRoots*nSel)
call GetMem('CTEMP','ALLO','REAL',kctemp,ndet)
call GetMem('SIGTEM','ALLO','REAL',ksigtemp,ndet)
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
    call Load_CI_vec(mRoot,nConf,Work(iVec1),LuDavid)
    if (iprlev >= DEBUG) then
      lPrint = min(nConf,200)
      write(u6,'(1X,A,I2,A,I2)') 'CI vector, iter =',it_ci,' mRoot =',mRoot
      write(u6,'(1X,A)') '(max. 200 elements)'
      write(u6,'(1X,A)') '-----------------------------'
      call dVcPrt(' ',' ',Work(iVec1),lPrint)
    end if

    call Timing(Rolex_1,Swatch,Swatch,Swatch)
    if (DOFARO) then
      ! determinant wavefunctions
      allocate(sgm(ndeta,ndetb))
      allocate(psi(ndeta,ndetb))

      call DCOPY_(NCONF,[Zero],0,WORK(IVECSVC),1)
      call REORD2(MY_NORB,NACTEL,1,0,IWORK(KICONF(1)),IWORK(KCFTP),WORK(IVEC1),WORK(IVECSVC),IWORK(IVKCNF))
      call CITRANS_SORT('C',WORK(IVECSVC),WORK(IVEC2))
      PSI = Zero
      call CITRANS_CSF2SD(WORK(IVEC2),PSI)
      SGM = Zero
      call SIGMA_UPDATE(HTU,GTUVX,SGM,PSI)
      call CITRANS_SD2CSF(SGM,WORK(IVEC2))
      call CITRANS_SORT('O',WORK(IVEC2),WORK(IVECSVC))
      call Reord2(my_norb,NACTEL,1,1,iWork(KICONF(1)),iWork(KCFTP),Work(IVECSVC),Work(IVEC2),iWork(ivkcnf))

      if (iprlev >= DEBUG) then
        FP = DNRM2_(NCONF,WORK(IVEC2),1)
        write(u6,'(1X,A,F21.14)') 'sigma dnrm2_(faroald): ',FP
      end if

      ! free the arrays
      deallocate(sgm,psi)
    else
      ! Convert the CI-vector from CSF to Det. basis
      call dcopy_(nconf,work(ivec1),1,work(kctemp),1)
      call dcopy_(ndet,[Zero],0,work(ksigtemp),1)
      call csdtvc(work(kctemp),work(ksigtemp),1,work(kdtoc),iwork(kicts(1)),stSym,1)
      call dcopy_(ndet,[Zero],0,work(ksigtemp),1)
      c_pointer = kctemp
      ! Calling Lucia to determine the sigma vector
      call Lucia_Util('Sigma',iDummy,iDummy,Dummy)
      ! Set mark so densi_master knows that the Sigma-vector exists on disk.
      iSigma_on_disk = 1
      call CSDTVC(work(iVec2),work(kctemp),2,work(kdtoc),iWork(kicts(1)),stSym,1)

      if (iprlev >= DEBUG) then
        FP = DNRM2_(NCONF,WORK(IVEC2),1)
        write(u6,'(1X,A,F21.14)') 'sigma dnrm2_(lucia):   ',FP
      end if
    end if

    ! Add ECORE_HEX (different from zero when particle-hole formalism used)
    ECORE_HEX = GET_ECORE()
    call daxpy_(nconf,ecore_hex,work(iVec1),1,work(iVec2),1)
    ! imings on generation of the sigma vector
    call Timing(Rolex_2,Swatch,Swatch,Swatch)
    Rolex_2 = Rolex_2-Rolex_1
    Rolex_3 = Rolex_3+Rolex_2

    if (iprlev >= DEBUG) then
      lPrint = min(nConf,200)
      write(u6,*) ' '
      write(u6,'(1X,A,I2,A,I2)') 'sigma vector, iter =',it_ci,' mRoot =',mRoot
      write(u6,'(1X,A)') '(max. 200 elements)'
      write(u6,'(1X,A)') '--------------------------------'
      call dVcPrt(' ',' ',Work(iVec2),lPrint)
    end if
    call Save_Sig_vec(mRoot,nConf,Work(iVec2),LuDavid)
  end do
  ! Sigma vectors (it_ci,mroot) have been computed, for mroot=1..lroots
  !---------------------------------------------------------------------
  ! compute Hsmall and Ssmall
  ! These are Hsmall(jtrial,ktrial), where jtrial is (jter,jroot), and
  ! ktrial is (kter,kroot), and similar Ssmall.
  ! jtrial=1..mxKeep*lroots correspond to jter=it_ci-mxKeep+1..it_ci
  ! (Fewer, at the beginning)

  do jRoot=1,nvec
    call Load_CI_vec(jRoot,nConf,Work(iVec1),LuDavid)
    call Load_Sig_vec(jRoot,nConf,Work(iVec2),LuDavid)
    do kRoot=1,jRoot
      call Load_CI_vec(kRoot,nConf,Work(iVec3),LuDavid)
      ij = kRoot+(jRoot*jRoot-jRoot)/2
      Sji = dDot_(nConf,Work(iVec1),1,Work(iVec3),1)
      Hji = dDot_(nConf,Work(iVec2),1,Work(iVec3),1)
      Work(iSs+ij-1) = Sji
      Work(iHs+ij-1) = Hji
    end do
  end do
  ntrial = nvec
  if (iprlev >= DEBUG) then
    call TriPrt('Hsmall',' ',Work(iHs),ntrial)
    call TriPrt('Ssmall',' ',Work(iSs),ntrial)
  end if
  ! Hsmall and Ssmall have been computed (ntrial x ntrial, in triangular
  ! storage.)
  !---------------------------------------------------------------------

  ! solve secular equation HC=SCE.

  ! PAM2009 nBasVec on input = min(ntrial,nconf)
  ! nBasVec returned as nr of orthonormal solutions to HC=SCE
  nBasVec = nConf
  call HCSCE(ntrial,Work(iHs),Work(iSs),Work(iCs),Work(iEs),nBasVec)
  if (nBasVec < lRoots) then
    write(u6,*) 'David: nBasVec less than lRoots'
    write(u6,*) 'nBasvec, lRoots = ',nBasVec,lRoots
    if (ICIRST == 1) write(u6,*) 'CIREstart was used. Check the number of roots in the previous calculation'
    call Abend()
  end if
  if (iprlev >= DEBUG) then
    call dVcPrt('Eigenvalues of Hsmall',' ',Work(iEs),ntrial)
    call RecPrt('Eigenvectors of Hsmall',' ',Work(iCs),ntrial,ntrial)
  end if
  !---------------------------------------------------------------------
  ! compute the current 'best' CI, sigma and residual vector

  ! CI vector is Work(iVec1)
  ! sigma vector is saved in Work(iVec2)
  ! residual vector is saved in Work(iVec3)
  do mRoot=1,lRoots
    ! initialize 'best' CI and sigma vector
    call dCopy_(nConf,[Zero],0,Work(iVec1),1)
    call dCopy_(nConf,[Zero],0,Work(iVec2),1)
    ! accumulate contributions
    do jRoot=1,nvec
      Cik = Work(iCs-1+jRoot+(mRoot-1)*ntrial)
      call Load_CI_vec(jRoot,nConf,Work(iVec3),LuDavid)
      call Daxpy_(nConf,Cik,Work(iVec3),1,Work(iVec1),1)
      call Load_Sig_vec(jRoot,nConf,Work(iVec3),LuDavid)
      call Daxpy_(nConf,Cik,Work(iVec3),1,Work(iVec2),1)
    end do
    RR = dDot_(nConf,Work(iVec1),1,Work(iVec1),1)
    scl = One/sqrt(RR)
    call DScal_(nConf,scl,Work(iVec1),1)
    call DScal_(nConf,scl,Work(iVec2),1)
    call Save_tmp_CI_vec(mRoot,nConf,Work(iVec1),LuDavid)
    call Save_tmp_Sig_vec(mRoot,nConf,Work(iVec2),LuDavid)
    ! compute residual vector
    E0 = Work(iEs+mRoot-1)
    call dCopy_(nConf,Work(iVec2),1,Work(iVec3),1)
    call daxpy_(nConf,-E0,Work(iVec1),1,Work(iVec3),1)
    ! save current best energy and residual
    RR = dDot_(nConf,Work(iVec3),1,Work(iVec3),1)
    CI_conv(1,mroot,it_ci) = E0
    CI_conv(2,mroot,it_ci) = sqrt(RR)
    ! print vectors
    if (iprlev >= DEBUG) then
      lPrint = min(nConf,200)
      write(u6,'(1X,A,I2,A,I2)') 'new best CI vector, iter =',it_ci,' mRoot =',mRoot
      write(u6,'(1X,A)') '(max. 200 elements)'
      write(u6,'(1X,A)') '--------------------------------------'
      call dVcPrt(' ',' ',Work(iVec1),lPrint)
      write(u6,'(1X,A,I2,A,I2)') 'new best sigma vector, iter =',it_ci,' mRoot =',mRoot
      write(u6,'(1X,A)') '(max. 200 elements)'
      write(u6,'(1X,A)') '-----------------------------------------'
      call dVcPrt(' ',' ',Work(iVec2),lPrint)
      write(u6,'(1X,A,I2,A,I2)') 'new residual vector, iter =',it_ci,' mRoot =',mRoot
      write(u6,'(1X,A)') '(max. 200 elements)'
      write(u6,'(1X,A)') '-----------------------------------------'
      call dVcPrt(' ',' ',Work(iVec3),lPrint)
    end if
    ! to improve the preconditioner select all elements in the
    ! subspace of the explicit Hamiltonian
    if (nSel > 1) then
      iOff = (mRoot-1)*nSel
      do i=1,nSel
        iConf = iSel(i)
        Work(iScr3+iOff+i-1) = Work(iVec3+iConf-1)
        Work(iScr4+iOff+i-1) = Work(iVec1+iConf-1)
      end do
    end if
  end do
  ! Current best CI & Sigma vectors have been stored in a temporary place
  ! for mroot=1..lroots.
  ! Also, the selected elements of the CI and Sigma vectors have been
  ! saved at Work(iScr3)(Sigma)  and Work(iScr4)(CI)
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
  if (iConv >= lRoots-hroots) goto 100
  !---------------------------------------------------------------------
  ! compute correction vectors q1 = r/(E0-H) and q2 = c/(E0-H)

  nleft = lRoots-nconverged
  if (nSel > 1) then
    ioff = nconverged*nSel
    call DGEMM_('T','N',nSel,nleft,nSel,One,ExplV,nSel,Work(iScr3+ioff),nSel,Zero,Work(iScr5),nSel)
    do mRoot=nconverged+1,lRoots
      E0 = Work(iEs+mRoot-1)
      iOff = (mRoot-nconverged-1)*nSel
      do i=1,nSel
        Z = E0-ExplE(i)
        if (abs(Z) < 0.001_wp) Z = 0.001_wp
        Work(iScr5+iOff+i-1) = Work(iScr5+iOff+i-1)/Z
      end do
    end do
    call DGEMM_('N','N',nSel,nleft,nSel,One,ExplV,nSel,Work(iScr5),nSel,Zero,Work(iScr3),nSel)
    ioff = nconverged*nSel
    call DGEMM_('T','N',nSel,nleft,nSel,One,ExplV,nSel,Work(iScr4+ioff),nSel,Zero,Work(iScr5),nSel)
    do mRoot=nconverged+1,lRoots
      E0 = Work(iEs+mRoot-1)
      iOff = (mRoot-nconverged-1)*nSel
      do i=1,nSel
        Z = E0-ExplE(i)
        if (abs(Z) < 0.001_wp) Z = 0.001_wp
        Work(iScr5+iOff+i-1) = Work(iScr5+iOff+i-1)/Z
      end do
    end do
    call DGEMM_('N','N',nSel,nleft,nSel,One,ExplV,nSel,Work(iScr5),nSel,Zero,Work(iScr4),nSel)
  end if
  !---------------------------------------------------------------------
  do mRoot=nconverged+1,lRoots
    E0 = -Work(iEs+mRoot-1)
    call Load_tmp_Sig_vec(mRoot,nConf,Work(iVec1),LuDavid)
    call Load_tmp_CI_vec(mRoot,nConf,Work(iVec2),LuDavid)
    call daxpy_(nConf,E0,Work(iVec2),1,Work(iVec1),1)
    call Load_H_diag(nConf,Work(iVec3),LuDavid)
    E0 = Work(iEs+mRoot-1)
    do i=0,nConf-1
      Z = E0-Work(iVec3+i)
      if (abs(Z) < 1.0e-4_wp) Z = 1.0e-4_wp
      Work(iVec3+i) = Work(iVec1+i)/Z
    end do
    if (nSel > 1) then
      iOff = (mRoot-nconverged-1)*nSel
      do i=1,nSel
        iConf = iSel(i)
        Work(iVec3+iConf-1) = Work(iScr3+iOff+i-1)
      end do
    end if
    Alpha(mRoot) = dDot_(nConf,Work(iVec3),1,Work(iVec2),1)
    call Load_H_diag(nConf,Work(iVec3),LuDavid)
    E0 = Work(iEs+mRoot-1)
    do i=0,nConf-1
      Z = E0-Work(iVec3+i)
      if (abs(Z) < 1.0e-4_wp) Z = 1.0e-4_wp
      Work(iVec3+i) = Work(iVec2+i)/Z
    end do
    if (nSel > 1) then
      iOff = (mRoot-nconverged-1)*nSel
      do i=1,nSel
        iConf = iSel(i)
        Work(iVec3+iConf-1) = Work(iScr4+iOff+i-1)
      end do
    end if
    Beta(mRoot) = dDot_(nConf,Work(iVec3),1,Work(iVec2),1)
  end do
  !---------------------------------------------------------------------

  ! compute correction vectors q3 = (r-E1*q2)/(E0-H)
  if (nSel > 1) then
    do mRoot=nconverged+1,lRoots
      call Load_tmp_Sig_vec(mRoot,nConf,Work(iVec1),LuDavid)
      call Load_tmp_CI_vec(mRoot,nConf,Work(iVec2),LuDavid)
      E0 = -Work(iEs+mRoot-1)
      call daxpy_(nConf,E0,Work(iVec2),1,Work(iVec1),1)
      E1 = -Alpha(mRoot)/Beta(mRoot)
      call daxpy_(nConf,E1,Work(iVec2),1,Work(iVec1),1)
      iOff = (mRoot-nconverged-1)*nSel
      do i=1,nSel
        iConf = iSel(i)
        Work(iScr3+iOff+i-1) = Work(iVec1+iConf-1)
      end do
    end do
    call DGEMM_('T','N',nSel,nleft,nSel,One,ExplV,nSel,Work(iScr3),nSel,Zero,Work(iScr5),nSel)
    do mRoot=nconverged+1,lRoots
      E0 = Work(iEs+mRoot-1)
      iOff = (mRoot-nconverged-1)*nSel
      do i=1,nSel
        Z = E0-ExplE(i)
        if (abs(Z) < 0.001_wp) Z = 0.001_wp
        Work(iScr5+iOff+i-1) = Work(iScr5+iOff+i-1)/Z
      end do
    end do
    call DGEMM_('N','N',nSel,nleft,nSel,One,ExplV,nSel,Work(iScr5),nSel,Zero,Work(iScr3),nSel)
  end if
  ! move the index of CI_vec
  istart = istart+nnew
  istart = mod(istart,nkeep-n_Roots)

  nnew = 0
  do mRoot=nconverged+1,lRoots
    call Load_tmp_Sig_vec(mRoot,nConf,Work(iVec1),LuDavid)
    call Load_tmp_CI_vec(mRoot,nConf,Work(iVec2),LuDavid)
    E0 = -Work(iEs+mRoot-1)
    call daxpy_(nConf,E0,Work(iVec2),1,Work(iVec1),1)
    E1 = -Alpha(mRoot)/Beta(mRoot)
    call daxpy_(nConf,E1,Work(iVec2),1,Work(iVec1),1)
    call Load_H_diag(nConf,Work(iVec3),LuDavid)
    E0 = Work(iEs+mRoot-1)
    do i=0,nConf-1
      Z = E0-Work(iVec3+i)
      if (abs(Z) < 1.0e-4_wp) Z = 1.0e-4_wp
      Work(iVec3+i) = Work(iVec1+i)/Z
    end do
    if (nSel > 1) then
      iOff = (mRoot-nconverged-1)*nSel
      do i=1,nSel
        iConf = iSel(i)
        Work(iVec3+iConf-1) = Work(iScr3+iOff+i-1)
      end do
    end if
    ! Orthonormalize wrt previous vectors
    updsiz = dnrm2_(nconf,Work(iVec3),1)
    scl = One/updsiz
    call DScal_(nConf,scl,Work(iVec3),1)
    do jRoot=lRoots+1,min(nvec,nkeep-nconverged)
      call Load_CI_vec(jRoot,nConf,Work(iVec2),LuDavid)
      ovl = dDot_(nConf,Work(iVec3),1,Work(iVec2),1)
      call daxpy_(nConf,-ovl,Work(iVec2),1,Work(iVec3),1)
    end do
    updsiz = dnrm2_(nconf,Work(iVec3),1)
    if (updsiz > 1.0e-6_wp) then
      scl = One/updsiz
      call DScal_(nConf,scl,Work(iVec3),1)
      nnew = nnew+1
      nvec = nvec+1
      nvec = min(nvec,nkeep)
      call Save_CI_vec(lRoots+mRoot-nconverged,nConf,Work(iVec3),LuDavid)
    end if
  end do
  !---------------------------------------------------------------------
  ! move the current best CI and sigma vectors to the first place
  ! in the list of retained CI vectors
  do mRoot=1,lRoots
    call Load_tmp_CI_vec(mRoot,nConf,Work(iVec1),LuDavid)
    call Save_CI_vec(mRoot,nConf,Work(iVec1),LuDavid)
    call Load_tmp_Sig_vec(mRoot,nConf,Work(iVec1),LuDavid)
    call Save_Sig_vec(mRoot,nConf,Work(iVec1),LuDavid)
  end do

! end of the long loop over iterations
end do
!=======================================================================

mxItr = min(mxCiIt,mxItr+12)
if (IPRLEV >= USUAL) then
  write(u6,*) '       No convergence in the CI section: MAXJT will be increased to ',mxItr
end if
Rc_CI = 16
nItr = nItr-1
! deallocate local temporary vectors
100 continue
call GetMem('CTEMP','Free','REAL',kctemp,ndet)
call GetMem('SIGTEM','Free','REAL',ksigtemp,ndet)
call GetMem('Vector1','Free','Real',iVec1,ndet)
call GetMem('Vector2','Free','Real',iVec2,ndet)
call GetMem('Vector3','Free','Real',iVec3,ndet)
call GetMem('Esmall','Free','Real',iEs,l1)
call GetMem('Hsmall','Free','Real',iHs,l3)
call GetMem('Ssmall','Free','Real',iSs,l3)
call GetMem('Csmall','Free','Real',iCs,l2)
call GetMem('Scr1','Free','Real',iScr1,l2)
call GetMem('Scr2','Free','Real',iScr2,l2)
call GetMem('Scr3','Free','Real',iScr3,lRoots*nSel)
call GetMem('Scr4','Free','Real',iScr4,lRoots*nSel)
call GetMem('Scr5','Free','Real',iScr5,lRoots*nSel)
if (DoFaro) then
  call GetMem('CIVEC','Free','Real',IVECSVC,nconf)
  call getmem('kcnf','free','inte',ivkcnf,nactel)
end if

call Timing(Alfex_2,Swatch,Swatch,Swatch)
Alfex_2 = Alfex_2-Alfex_1
Alfex_3 = Alfex_3+Alfex_2

return

end subroutine David5
