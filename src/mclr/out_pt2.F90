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

subroutine Out_Pt2(iKapDisp,iCIDisp)

use Index_Functions, only: iTri, nTri_Elem
use ipPage, only: ipclose, ipget, ipin, W
use MCLR_Data, only: CMO, ESTERR, ipCI, ipCM, ipMat, IRLXROOT, ISNAC, ISTATE, LuJob, LuPT2, LuTEMP, n1Dens, n2Dens, nA, NACSTATES, &
                     nConf1, nDens, nDensC, OVERRIDE
use input_mclr, only: iRoot, iTOC, nAsh, nBas, nCSF, nDisp, nIsh, nRoots, nSym, ntAsh, ntBtri, PT2, State_Sym
use PCM_grad, only: DSCFAO, DSSAO, PCM_grad_D2v, PCMSCFAO, PCMSCFMO
use rctfld_module, only: lRF
use ISRotation, only: ISR_final2
use dmrginfo, only: DoDMRG, LRRAS2, RGRAS2
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iKapDisp(nDisp), iCiDisp(nDisp)
integer(kind=iwp) :: i, iA, iAA, iBas, iDisk, idum(7,8), ii, iikl, ij, ij1, ij2, ijkl, ijlk, iLen, iMax, iOff, ipCIP, iR, iS, &
                     iSym, j, jA, jAA, jDisk, ji2, jikl, jilk, k, kl, kl1, kl2, klRow, l, lk2, LuDens, LuTmp, nBas_Tot, nBasI, &
                     nBuf, nDim, nDLMO, nG1, nG2, nLCMO, nNac, nPLMO, nTot1
real(kind=wp) :: DM(3), rdum(1), Val
logical(kind=iwp) :: CI, Is_Roots_Set
character(len=80) :: Note
character(len=16) :: mstate
character(len=8) :: Method
real(kind=wp), allocatable :: CMON(:), Conn(:), D1(:), D_CI(:), D_K(:), DAO(:), DMs(:,:), DTmp(:), G1m(:), G1q(:), K1(:), K2(:), &
                              OCCU(:), P1(:), P_CI(:), Temp(:), Tmp(:), tmpDe(:,:), tmpDeM(:,:), tmpP(:), tmpPM(:,:,:,:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
isym = 1
CI = .true.
call Setup_MCLR(iSym)
nbas_tot = sum(nbas(1:nsym))
nDLMO = sum(nash(1:nsym))
nLCMO = sum(nbas(1:nsym)**2)
ntot1 = 0
do is=1,nsym
  ntot1 = ntot1+nTri_Elem(nbas(is))
end do
nNAC = nTri_Elem(nDLMO)
nDLMO = nTri_Elem(nDLMO)
nPLMO = nTri_Elem(nDLMO)

call mma_allocate(K1,nDens,Label='K1')
call mma_allocate(K2,nDens,Label='K2')
call mma_allocate(DAO,nDens,Label='DAO')
call mma_allocate(D_CI,n1Dens,Label='D_CI')
call mma_allocate(D1,n1Dens,Label='D1')
call mma_allocate(P_CI,n2Dens,Label='P_CI')
call mma_allocate(P1,n2Dens,Label='P1')
call mma_allocate(Conn,nDens,Label='Conn')
Conn(:) = Zero
call mma_allocate(OCCU,nbas_tot,Label='OCCU')
call mma_allocate(CMON,nDens,Label='CMON')
! OBS nBuf might not be def.
call mma_MaxDBLE(nBuf)
call mma_allocate(Dtmp,nDens,Label='DTmp')

! 1)   CI Part

! All multipliers are introduced as densities

if (CI) then
  nconf1 = ncsf(State_sym)
  ilen = nconf1*nroots ! nroot = # of roots in SA
  ipcip = ipget(nconf1*nroots)
  iDisk = iCIDisp(1)
  call ipin(ipCIp)
  call dDaFile(LuTemp,2,W(ipCIp)%A,iLen,iDisk)

  ! Calculate the densities that correct the nonvariational CI stuff

  call CIDens_sa(.true.,ipCIp,ipCI,State_sym,State_sym,P_CI,D_CI) ! \bar{d} and \bar{D}

  ! ====================================================================
  if (doDMRG) then  ! yma
    call dmrg_dim_change_mclr(LRras2,ntash,0)
    call dmrg_dim_change_mclr(RGras2,ndim,0)

    call mma_allocate(tmpDe,ndim,ndim,Label='TmpDe')
    call mma_allocate(tmpP,nTri_Elem(ndim**2),Label='tmpP')
    call mma_allocate(tmpDeM,ntash,ntash,Label='tmpDeM')
    call mma_allocate(tmpPM,ntash,ntash,ntash,ntash,Label='tmpPM')
    tmpDe = Zero
    tmpP = Zero
    tmpDeM = Zero
    tmpPM = Zero

    ij = 0
    do i=1,ntash
      do j=1,ntash
        ij = ij+1
        if (abs(D_CI(ij)) < 1.0e-12_wp) D_CI(ij) = Zero
        tmpDeM(i,j) = D_CI(ij)
      end do
    end do

    ij = 0
    do i=1,ndim
      do j=1,ndim
        ij = ij+1
        if ((i > ntash) .or. (j > ntash)) then
          tmpDe(i,j) = Zero
        else
          tmpDe(i,j) = tmpDeM(i,j)
        end if
      end do
    end do

    do i=1,ntash
      do j=1,ntash
        do k=1,ntash
          do l=1,ntash
            ij1 = ntash*(i-1)+j
            ij2 = ntash*(j-1)+i
            kl1 = ntash*(k-1)+l
            kl2 = ntash*(l-1)+k
            if (ij1 >= kl1) then
              if (abs(P_CI(itri(ij1,kl1))) < 1.0e-12_wp) P_CI(iTri(ij1,kl1)) = Zero
              tmpPM(i,j,k,l) = P_CI(iTri(ij1,kl1))
            end if
          end do
        end do
      end do
    end do

    do i=1,ndim
      do j=1,ndim
        do k=1,ndim
          do l=1,ndim
            ij1 = ndim*(i-1)+j
            ij2 = ndim*(j-1)+i
            kl1 = ndim*(k-1)+l
            kl2 = ndim*(l-1)+k
            if (ij1 >= kl1) then
              if ((i > ntash) .or. (j > ntash) .or. (k > ntash) .or. (l > ntash)) then
                tmpP(iTri(ij1,kl1)) = Zero
              else
                tmpP(iTri(ij1,kl1)) = tmpPM(i,j,k,l)
              end if
            end if
          end do
        end do
      end do
    end do

    D_CI(1:ndim**1) = pack(tmpDe(:,:),.true.)
    P_CI(1:n2Dens) = tmpP(1:n2Dens)
    call mma_deallocate(tmpDe)
    call mma_deallocate(tmpDeM)
    call mma_deallocate(tmpP)
    call mma_deallocate(tmpPM)
    call dmrg_dim_change_mclr(RGras2,ntash,0)
  end if
  ! ====================================================================

  ! Some administrative shit

  ! Store densities in triangular form

  do i=1,ntAsh
    do j=1,i
      D1(iTri(i,j)) = D_CI((i-1)*ntash+j)
    end do
  end do

  do i=1,ntAsh
    do j=1,i
      ij = iTri(i,j)
      ij2 = i+(j-1)*ntash
      ji2 = j+(i-1)*ntash
      do k=1,i
        do l=1,k
          kl = iTri(k,l)
          kl2 = k+(l-1)*ntash
          lk2 = l+(k-1)*ntash
          ijkl = iTri(ij2,kl2)
          jikl = iTri(ji2,kl2)
          ijlk = iTri(ij2,lk2)
          jilk = iTri(ji2,lk2)
          P1(iTri(ij,kl)) = Quart*(P_CI(ijkl)+P_CI(jikl)+P_CI(ijlk)+P_CI(jilk))
        end do
      end do
    end do
  end do

  do K=1,NTASH
    do L=1,K
      KL = iTri(K,L)
      KLROW = nTri_Elem(KL-1)
      if (L == K) then
        IMAX = K
      else
        IMAX = K-1
      end if
      do I=1,IMAX
        II = nTri_Elem(I)
        IIKL = KLROW+II
        P1(IIKL) = P1(IIKL)*Half
      end do
    end do
  end do

  !do i=1,ntAsh
  !  do j=1,i
  !    ij = iTri(i,j)
  !    ij2 = i+(j-1)*ntash
  !    ji2 = j+(i-1)*ntash
  !    do k=1,ntAsh
  !      do l=1,k
  !        kl = iTri(k,l)
  !        kl2 = k+(l-1)*ntash
  !        ijkl = iTri(ij2,kl2)
  !        jikl = iTri(ji2,kl2)
  !        fact = Half
  !        if ((ij >= kl) .and. (k == l)) fact = Quart
  !        if ((ij < kl) .and. (i == j)) fact = Quart
  !        P1(iTri(ij,kl)) = fact*(P_CI(ijkl)+P_CI(jikl))
  !      end do
  !    end do
  !  end do
  !end do
  !if (debug) Call triprt('P1',' ',P1,nTri_Elem(ntash))

  ! Write the 'bar' densities to disk,  not symmetry blocked.

  !call Put_dArray('DLMO',D1,ndim1) ! \bar{D} triangular  ! yma
  !call Put_dArray('PLMO',P1,ndim2) ! \bar{d} triangular  ! yma

  call Put_dArray('DLMO',D1,nDLMO) ! \bar{D} triangular
  call Put_dArray('PLMO',P1,nPLMO) ! \bar{d} triangular

end if

! 2) Orbital response
!    ================

! Read in from disk

iDisk = iKapDisp(1)
call dDaFile(LuTemp,2,K1,nDensC,iDisk) ! Read \bar{kappa}
call Uncompress(K1,K2,1)

! If we want to estimate the error

if (esterr) then
  !do iestate=1,lroots
  !  call calcerr(K2,iestate)
  !end Do
  call calcerr(K2,istate)
end if

! First we fix the renormalization contribution

call mma_allocate(D_K,nLCMO,Label='D_K')
call Get_dArray_chk('FockOcc',D_K,nLCMO)
! Calculates the effective Fock matrix
call Make_Conn(Conn,K2,P_CI,D_CI)   !D_CI not changed
Conn(1:nTot1) = Conn(1:nTot1)+D_K(1:nTot1)
!Conn(:) = D_K(1:nDens)
if (PT2) then
  ! Add the WLag term (will be contracted with overlap
  ! derivative) computed in CASPT2
  do i=1,nTot1
    read(LuPT2,*) Val
    Conn(i) = Conn(i)+Val
  end do
end if
call Put_dArray('FockOcc',Conn,nTot1)

! Transposed one index transformation of the density
! (only the inactive density to store it separately)

call OITD(K2,1,DAO,Dtmp,.false.)

if (PT2) then
  ! For gradient calculation. D^var couples with inactive
  ! orbitals only, so D0(1,4) (in integral_util/prepp.f), which
  ! couples with active orbitals has to be modified,
  do iSym=1,nSym
    nBasI = nBas(iSym)
    do iI=1,nBasI
      do iJ=1,nBasI
        read(LuPT2,*) Val
        DAO(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI) = DAO(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)+Val
      end do
    end do
  end do
  ! The PT2 density will be used later again.
  do iI=1,sum(nBas(1:nSym)**2)
    backspace LuPT2
  end do
end if

! Transformation to AO basis (covariant)

! Transforms to AO differently dep on last arg.

call TCMO(DAO,1,-2)

! Fold AO density and write to disk
! Mult all terms that are not diag by 2

call FOLD2(nsym,nbas,DAO,K1)

call Put_dArray('DLAO',K1,ntot1)

! Now with active density too, to form the variational density

! gives \tilde{D}
call OITD(K2,1,D_K,Dtmp,.true.)

do iS=1,nsym

  ! C*\tilde{\kappa} --> ipDAO

  if (nBas(is) >= 1) call DGEMM_('N','N',NBAS(is),NBAS(is),NBAS(is),One,CMO(ipCM(is)),NBAS(is),K2(ipmat(is,is)),NBAS(is),Zero, &
                                 DAO(ipCM(is)),NBAS(is))
end do

call Put_dArray('LCMO',DAO,nLCMO)

if (doDMRG) then  ! yma
  call dmrg_dim_change_mclr(RGras2,ntash,0)
  nash(:) = RGras2(:)
end if

if (isNAC) then
  ng1 = nNAC
  call mma_allocate(G1q,ng1,Label='G1q')
  call Get_cArray('Relax Method',Method,8)
  if (Method == 'MSPDFT') then
    call Get_DArray('D1MOt',G1q,ng1)
  else
    call Get_dArray_chk('D1mo',G1q,ng1)
  end if
  iR = 0 ! set to dummy value.
else
  iR = iroot(istate)
  jdisk = itoc(3)
  ng1 = nTri_Elem(ntash)
  ng2 = nTri_Elem(ng1)
  call mma_allocate(G1q,n1Dens,Label='G1q')

  ! Read active one el dens for state j from JOBIPH and store in G1q

  call Get_cArray('Relax Method',Method,8)
  if (Method == 'MSPDFT') then
    call Get_DArray('D1MOt',G1q,ng1)
  else
    do i=1,iR-1  ! Dummy read until state j
      call dDaFile(LUJOB,0,rdum,ng1,jDisk)
      call dDaFile(LUJOB,0,rdum,ng1,jDisk)
      call dDaFile(LUJOB,0,rdum,ng2,jDisk)
      call dDaFile(LUJOB,0,rdum,ng2,jDisk)
    end do
    call dDaFile(LUJOB,2,G1q,ng1,jDisk)

    if (PT2 .and. (nRoots > 1)) call Get_dArray('D1mo',G1q,ng1)

  end if
end if

! Construct a variationally stable density matrix. In MO

! D_eff = D^j + \tilde{D} +\bar{D}
! D_K = (G1q + inact) + D_K + D_CI

if (isNAC) then

  ! For NAC, first build DAO and then DAO_var

  do is=1,nSym
    ! Note: no inactive part for transition densities
    do iA=1,nash(is)
      do jA=1,nash(is)
        i = iA+nish(is)
        j = jA+nish(is)
        iAA = iA+na(is)
        jAA = jA+na(is)
        D_K(ipmat(is,is)+i-1+(j-1)*nbas(is)) = D_K(ipmat(is,is)+i-1+(j-1)*nbas(is))+D_CI(iAA+(jAA-1)*ntash)+G1q(iTri(iAA,jAA))
      end do
    end do
  end do

  if (PT2) then
    ! PT2 density (in MO)
    do iSym=1,nSym
      nBasI = nBas(iSym)
      do iI=1,nBasI
        do iJ=1,nBasI
          read(LuPT2,*) Val
          D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI) = D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)+Val
        end do
      end do
    end do
    ! PT2C density (in MO)
    do iSym=1,nSym
      nBasI = nBas(iSym)
      do iI=1,nBasI
        do iJ=1,nBasI
          read(LuPT2,*) Val
          D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI) = D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)+Val*Quart
          D_K(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI) = D_K(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI)+Val*Quart
        end do
      end do
    end do
  end if
  call mma_allocate(Temp,nBuf/2,Label='Temp')
  call NatOrb_MCLR(D_K,CMO,CMON,OCCU)
  call dmat_MCLR(CMON,OCCU,Temp)
  call Put_dArray('D1aoVar',Temp,nTot1)
  call mma_deallocate(Temp)

  ! Transform the antisymmetric transition density matrix to AO
  !  (there is no guarantee the symmetry will work here)

  iDisk = 0
  LuDens = 20
  call DaName(LuDens,'MCLRDENS')
  call dDaFile(LuDens,2,G1q,ng1,iDisk)
  call DaClos(LuDens)
  call mma_allocate(G1m,nDens,Label='G1m')
  G1m(:) = Zero
  ! Reconstruct the square matrix
  do is=1,nSym
    do iA=1,nash(is)
      i = iA+nish(is)
      iAA = iA+na(is)
      do jA=1,iA-1
        j = jA+nish(is)
        jAA = jA+na(is)
        G1m(ipmat(is,is)+i-1+(j-1)*nbas(is)) = G1q(iTri(iAA,jAA))
        G1m(ipmat(is,is)+j-1+(i-1)*nbas(is)) = -G1q(iTri(iAA,jAA))
      end do
      G1m(ipmat(is,is)+i-1+(i-1)*nbas(is)) = Zero
    end do
  end do

  if (PT2) then
    ! PT2C density (in MO) for CSF derivative
    do iSym=1,nSym
      nBasI = nBas(iSym)
      do iI=1,nBasI
        do iJ=1,nBasI
          read(LuPT2,*) Val
          G1m(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI) = G1m(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI)+Val
        end do
      end do
    end do
  end if
  ! Transform
  call TCMO(G1m,1,-2)
  ! Save the triangular form
  iOff = 0
  do is=1,nSym
    ibas = nbas(is)
    do i=1,ibas
      do j=1,i
        G1m(iOff+iTri(i,j)) = G1m(ipmat(is,is)+j-1+(i-1)*nbas(is))
      end do
    end do
    iOff = iOff+nTri_Elem(ibas)
  end do
  call Put_dArray('D1ao-',G1m,nTot1)
  call mma_deallocate(G1m)

else

  ! Normal SA gradient (no NAC)

  do is=1,nSym
    do i=1,nish(is)

      ! The inactive density

      D_K(ipmat(is,is)+i-1+(i-1)*nbas(is)) = D_K(ipmat(is,is)+i-1+(i-1)*nbas(is))+Two
    end do
    do iA=1,nash(is)
      do jA=1,nash(is)
        i = iA+nish(is)
        j = jA+nish(is)
        iAA = iA+na(is)
        jAA = jA+na(is)

        ! The active density G1q and \bar{D}

        D_K(ipmat(is,is)+i-1+(j-1)*nbas(is)) = D_K(ipmat(is,is)+i-1+(j-1)*nbas(is))+D_CI(iAA+(jAA-1)*ntash)+G1q(iTri(iAA,jAA))
      end do
    end do
  end do

  if (PT2) then
    ! Add PT2 density (in MO)
    do iSym=1,nSym
      nBasI = nBas(iSym)
      do iI=1,nBasI
        do iJ=1,nBasI
          read(LuPT2,*) Val
          D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI) = D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)+Val
        end do
      end do
    end do
    ! Also, PT2C density (in MO)
    ! This density couples with inactive orbitals only,
    ! while the above PT2 density couples with inactive+active
    ! orbitals.
    do iSym=1,nSym
      nBasI = nBas(iSym)
      do iI=1,nBasI
        do iJ=1,nBasI
          read(LuPT2,*) Val
          D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI) = D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)+Val*Quart
          D_K(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI) = D_K(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI)+Val*Quart
        end do
      end do
    end do
  end if

  ! Diagonalize the effective density to be able to use Prpt
  ! OCCU eigenvalues of eff dens
  ! CMON eigenvectors (new orb coef)

  call NatOrb_MCLR(D_K,CMO,CMON,OCCU)
  call mma_Allocate(Tmp,nBuf/2,Label='Tmp')
  call dmat_MCLR(CMON,OCCU,Tmp)
  call Put_dArray('D1aoVar',Tmp,nTot1)
  call mma_deallocate(Tmp)

  !call mma_allocate(TEMP,nNac,Label='TEMP')
  !call mma_allocate(tTmp,nNac,Label='tTmp')
  !call get_dArray_chk('D1mo',TEMP,nNac)
  !call get_dArray_chk('DLMO',tTmp,nNac)
  !TEMP(:) = TEMP(:)+tTmp(:)
  !call mma_deallocate(TEMP)
  !call mma_deallocate(tTmp)

  Note = 'var'
  LuTmp = 50
  LuTmp = IsFreeUnit(LuTmp)
  call WrVec('TMPORB',LuTmp,'O',nSym,nBas,nBas,rDum,OCCU,rDum,iDum,Note)
  call Prpt()

  !*********************************************************************
  ! There should now be dipole moments on the runfile which
  ! corresponds to the gradient of the energy w.r.t. the
  ! electric field. Let's update the list of values stored
  ! on the runfile.

  Is_Roots_Set = .false.
  call Qpg_iScalar('Number of roots',Is_Roots_Set)
  nRoots = 1
  if (Is_Roots_Set) call Get_iScalar('Number of roots',nRoots)

  if (nRoots /= 1) then
    !write(u6,*) 'iR=',iR
    call mma_allocate(DMs,3,nROOTS,Label='DMs')
    call Get_dArray('Last Dipole Moments',DMs,3*nRoots)
    !call RecPrt('Last Dipole Moments',' ',DMS,3,nRoots)
    call Get_dArray('Dipole Moment',DM,3)
    !call RecPrt('Dipole Moment',' ',DM,1,3)
    DMs(:,iR) = DM(:)
    !call RecPrt('Last Dipole Moments',' ',DMs,3,nRoots)
    call Put_dArray('Last Dipole Moments',DMs,3*nRoots)
    call mma_deallocate(DMs)
  end if
  !*********************************************************************
end if
call mma_deallocate(G1q)

!----- debug -----

if (doDMRG) then ! yma
  call dmrg_dim_change_mclr(LRras2,ntash,0)
  nash(:) = LRras2(:)
end if

! Write the effective active one el density to disk in the same format as g1q

!call mma_allocate(Deff_act,nDens,Label='Deff_act')
!Deff_act(:) = D_K(1:nDens)
!do is=1,nSym
!  do i=1,nish(is)
!
!    ! Subtract the inactive density
!
!    Deff_act(ipmat(is,is)+i-1+(i-1)*nbas(is)) = D_K(ipmat(is,is)+i-1+(i-1)*nbas(is))-Two
!  end do
!end do

!call Put_DEff(Deff_act,nDens)

!call mma_deallocate(Deff_act)

!--------------------------------------------------

! Diagonalize the effective density to be able to use Prpt
! OCCU eigenvalues of eff dens
! CMON eigenvectors (new orb coef)

!call NatOrb_MCLR(D_K,CMO,CMON,OCCU)
!call mma_allocate(Temp,nBuf/2,Label='Temp')
!call dmat_MCLR(CMON,OCCU,Temp)
!call Put_dArray('D1aoVar',Temp,nTot1)
!Note = 'var'
!LuTmp = 50
!LuTmp = IsFreeUnit(LuTmp)
!call WrVec('TMPORB',LuTmp,'O',nSym,nBas,nBas,Dum,OCCU,Dum,iDum,Note)
!call Prpt()

! Standard routine, Temp effective dens in AO

!call dmat_MCLR(CMON,OCCU,Temp)

!call Put_dArray('D1aoVar',Temp,nTot1)
!call mma_deallocate(Temp)

call Put_iScalar('SA ready',1)
if (isNAC) then
  write(mstate,'(1X,I7,",",I7)') NACStates(1),NACStates(2)
else
  write(mstate,'(I16)') irlxroot
end if
if (override) mstate(1:1) = '+'
call Put_cArray('MCLR Root',mstate,16)

!! In the end, it may be useful to generate ASCs using the relaxed
!! density matrix for non-equilibrium solvation model
if (lRF) then
  call Get_dArray('D1aoVar',DSSAO,ntBtri)
  call unfold(DSSAO,ntot1,DSCFAO,nLCMO,nSym,nBas)
  call PCM_grad_D2V(DSCFAO,PCMSCFMO,PCMSCFAO,.false.,.false.,.false.,.false.)
end if

call mma_deallocate(K1)
call mma_deallocate(K2)
call mma_deallocate(DAO)
call mma_deallocate(D_CI)
call mma_deallocate(D1)
call mma_deallocate(P_CI)
call mma_deallocate(P1)
call mma_deallocate(Conn)
call mma_deallocate(OCCU)
call mma_deallocate(CMON)
call mma_deallocate(Dtmp)
call mma_deallocate(D_K)
call ISR_final2()

call ipclose(-1)

end subroutine Out_Pt2
