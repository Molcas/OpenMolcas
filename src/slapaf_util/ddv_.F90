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

subroutine ddV_(Cart,mTtAtm,Hess,iANr,iTabBonds,iTabAtoms,nBonds,nMax,nHidden)

use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: nIrrep, iOper, VarR, VarT
use Slapaf_Info, only: ddV_Schlegel, iOptC, Magic_Bond
use ddvdt, only: A_Bend, A_Str, A_StrH, A_Trsn, aAV, alpha_vdW, B_Str, f_Const_Min, r_ref_vdW, rAV, rko, rkf, rkr, rkr_vdW, rkt, &
                 Rot_Const, Trans_Const
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Ten, Half, Angstrom, deg2rad
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mTtAtm, iANr(mTtAtm), nBonds, iTabBonds(3,nBonds), nMax, iTabAtoms(2,0:nMax,mTtAtm), nHidden
real(kind=wp), intent(in) :: Cart(3,mTtAtm)
real(kind=wp), intent(out) :: Hess(nTri_Elem(3*mTtAtm))
integer(kind=iwp) :: i, iAtom, iBond, iBondType, icoor, ij, iNb0, iNb1, iNb2, iNeighbor, ir, iSym, iTest, ixyz, jAtom, jBond, &
                     jBondType, jCoor, jNeighbor, jr, kAtom, kBond, kBondType, kNeighbor, kr, lAtom, lBond, lBondType, lr, mAtom, &
                     mr, nCoBond_j, nNeighbor, nNeighbor_i, nNeighbor_j, nNeighbor_k, nOrder
real(kind=wp) :: A35, aij, aik, ail, ajk, akl, alpha, ami, amj, beta, C(3,4), CosFi2, CosFi3, CosFi4, CosFi_Max, CosPhi, &
                 cosThetax, cosThetay, cosThetaz, Diff, dO1_dx1, dO1_dx2, dO1_dy1, dO1_dy2, dO1_dz1, dO1_dz2, dO2_dx1, dO2_dx2, &
                 dO2_dy1, dO2_dy2, dO2_dz1, dO2_dz2, dO3_dx1, dO3_dx2, dO3_dy1, dO3_dy2, dO3_dz1, dO3_dz2, Dum(1), f1, &
                 f2, f_const, f_const_min_, Fact, g_ij, g_jk, g_kl, g_vdW, g_vdW_ij, g_vdW_im, g_vdW_jk, g_vdW_jm, g_vdW_kl, gij, &
                 gim, gjm, gmm, Hxx, Hxy, Hxz, Hyy, Hyz, Hzz, r0, r0_vdW, r0_vdW_ij, r0_vdW_im, r0_vdW_jk, r0_vdW_jm, r0_vdW_kl, &
                 r0mi, r0mj, r1, Rab, RabCov, Rbc, RbcCov, rij(3), rij0, rij2, rik(3), rik0, rik2, ril(3), ril0, ril2, rjk(3), &
                 rjk0, rjk2, rkl(3), rkl0, rkl2, rL, rL2, rmi, rmi2, rmidotrmj, rmj, rmj2, RotAng, RotMat(3,3), RotVec(3), rrij, &
                 rZero, si(3), SinPhi, sj(3), sk(3), sl(3), sm(3), Tau, Test_zero, Thr_Line, ThrFi1, ThrFi2, tij, TMass, Trans(3), &
                 x(2), xij, xkl, xmi, xmj, xyz(3,4), y(2), yij, ykl, ymi, ymj, z(2), zij, zkl, zmi, zmj
logical(kind=iwp) :: Help, Invariant(3), MinBas
real(kind=wp), allocatable :: CurrXYZ(:,:), Grad(:,:,:), xMass(:)
integer(kind=iwp), external :: iTabRow, LHR, nCoBond
real(kind=wp), external :: CovRad, CovRadT, rMass
logical(kind=iwp), external :: Torsion_Check

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) 'ddV: nBonds=',nBonds
nqR = 0
nqB = 0
nqT = 0
nqO = 0
do iAtom=1,mTtAtm
  nNeighbor_i = iTabAtoms(1,0,iAtom)
  write(u6,*) 'iAtom,nNeighbor=',iAtom,nNeighbor_i
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
f_const_min_ = f_const_min*0.1_wp
f_const = Zero

MinBas = .false.
if (MinBas) then
  Fact = 1.3_wp
else
  Fact = One
end if
rZero = 1.0e-10_wp

Hess(:) = Zero
#ifdef _DEBUGPRINT_
n3 = 3*mTtAtm
call TriPrt(' In LNM: Hessian at start','(12f8.3)',Hess,n3)
call DiagMtrx_T(Hess,n3,iNeg)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(xMass,mTtAtm,Label='xMass')
call Get_Mass_All(xMass,mTtAtm-nHidden)
do iAtom=mTtAtm-nHidden+1,mTtAtm
  xMass(iAtom) = rMass(iANr(iAtom))
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Hessian for translational coordinates in case of a RF calculation.

if (VarT) then

  do ixyz=1,3
    Invariant(ixyz) = .false.
    iTest = 2**(ixyz-1)
    do iSym=0,nIrrep-1
      if (iOper(iSym) == iTest) Invariant(ixyz) = .true.
    end do
  end do
  if (.not. all(Invariant)) then

    Fact = One
    if (.not. VarR) Fact = 2.0e-2_wp

    TMass = sum(xMass)
    do iAtom=1,mTtAtm
      f1 = xMass(iAtom)/TMass
      do jAtom=1,iAtom-1
        f2 = xMass(jAtom)/TMass

        f_const = max(Trans_Const,f_const_Min_)
        gmm = Fact*f_const*f1*f2

        if (.not. Invariant(1)) Hess(LHR(1,iAtom,1,jAtom)) = Hess(LHR(1,iAtom,1,jAtom))+gmm
        if (.not. Invariant(2)) Hess(LHR(2,iAtom,2,jAtom)) = Hess(LHR(2,iAtom,2,jAtom))+gmm
        if (.not. Invariant(3)) Hess(LHR(3,iAtom,3,jAtom)) = Hess(LHR(3,iAtom,3,jAtom))+gmm

      end do

      f_const = max(Trans_Const,f_const_Min_)
      gmm = Fact*f_const*f1*f1

      if (.not. Invariant(1)) Hess(LHR(1,iAtom,1,iAtom)) = Hess(LHR(1,iAtom,1,iAtom))+gmm
      if (.not. Invariant(1)) Hess(LHR(2,iAtom,2,iAtom)) = Hess(LHR(2,iAtom,2,iAtom))+gmm
      if (.not. Invariant(1)) Hess(LHR(3,iAtom,3,iAtom)) = Hess(LHR(3,iAtom,3,iAtom))+gmm

    end do
#   ifdef _DEBUGPRINT_
    call TriPrt(' In LNM: Hessian after Translation','(6f12.7)',Hess,n3)
    call DiagMtrx_T(Hess,n3,iNeg)
#   endif
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Hessian for rotational coordinates

if (VarR) then

  do ixyz=1,3
    Invariant(ixyz) = .false.
    if (ixyz == 1) then
      iTest = 6
    else if (ixyz == 2) then
      iTest = 5
    else
      iTest = 3
    end if
    do iSym=0,nIrrep-1
      if (iOper(iSym) == iTest) Invariant(ixyz) = .true.
    end do
  end do
  if (.not. all(Invariant)) then

    !if (mTtAtm <= 2) then
    !  call WarningMessage(2,'Error in ddV')
    !  write(u6,*)
    !  write(u6,*) ' Warning!'
    !  write(u6,*) ' Rotational internal coordinates not implemented for fewer than 3 atoms!'
    !  write(u6,*) ' Add dummy atoms to your input and try again!'
    !  write(u6,*)
    !  call Quit(_RC_GENERAL_ERROR_)
    !end if
    call mma_allocate(Grad,3,3,mTtAtm,Label='Grad')
    call mma_allocate(CurrXYZ,3,mTtAtm,Label='CurrXYZ')
    do iAtom=1,mTtAtm
      if (iANr(iAtom) <= 0) xMass(iAtom) = 1.0e-10_wp
    end do
    nOrder = 1
    Trans(:) = Zero
    RotVec(:) = Zero
    CurrXYZ(:,:) = Cart(:,:)
    call RotDer(mTtAtm,xMass,CurrXYZ,Cart,Trans,RotAng,RotVec,RotMat,nOrder,Grad,dum)
    call mma_deallocate(CurrXYZ)

    do iAtom=1,mTtAtm
      dO1_dx1 = Grad(1,1,iAtom)
      dO2_dx1 = Grad(2,1,iAtom)
      dO3_dx1 = Grad(3,1,iAtom)
      dO1_dy1 = Grad(1,2,iAtom)
      dO2_dy1 = Grad(2,2,iAtom)
      dO3_dy1 = Grad(3,2,iAtom)
      dO1_dz1 = Grad(1,3,iAtom)
      dO2_dz1 = Grad(2,3,iAtom)
      dO3_dz1 = Grad(3,3,iAtom)
      if (Invariant(1)) then
        dO1_dx1 = Zero
        dO1_dy1 = Zero
        dO1_dz1 = Zero
      end if
      if (Invariant(2)) then
        dO2_dx1 = Zero
        dO2_dy1 = Zero
        dO2_dz1 = Zero
      end if
      if (Invariant(3)) then
        dO3_dx1 = Zero
        dO3_dy1 = Zero
        dO3_dz1 = Zero
      end if
      do jAtom=1,iAtom-1
        dO1_dx2 = Grad(1,1,jAtom)
        dO2_dx2 = Grad(2,1,jAtom)
        dO3_dx2 = Grad(3,1,jAtom)
        dO1_dy2 = Grad(1,2,jAtom)
        dO2_dy2 = Grad(2,2,jAtom)
        dO3_dy2 = Grad(3,2,jAtom)
        dO1_dz2 = Grad(1,3,jAtom)
        dO2_dz2 = Grad(2,3,jAtom)
        dO3_dz2 = Grad(3,3,jAtom)
        if (Invariant(1)) then
          dO1_dx2 = Zero
          dO1_dy2 = Zero
          dO1_dz2 = Zero
        end if
        if (Invariant(2)) then
          dO2_dx2 = Zero
          dO2_dy2 = Zero
          dO2_dz2 = Zero
        end if
        if (Invariant(3)) then
          dO3_dx2 = Zero
          dO3_dy2 = Zero
          dO3_dz2 = Zero
        end if

        f_const = max(Rot_Const,f_const_Min_)
        Hess(LHR(1,iAtom,1,jAtom)) = Hess(LHR(1,iAtom,1,jAtom))+f_const*(dO1_dx1*dO1_dx2+dO2_dx1*dO2_dx2+dO3_dx1*dO3_dx2)
        Hess(LHR(1,iAtom,2,jAtom)) = Hess(LHR(1,iAtom,2,jAtom))+f_const*(dO1_dx1*dO1_dy2+dO2_dx1*dO2_dy2+dO3_dx1*dO3_dy2)
        Hess(LHR(1,iAtom,3,jAtom)) = Hess(LHR(1,iAtom,3,jAtom))+f_const*(dO1_dx1*dO1_dz2+dO2_dx1*dO2_dz2+dO3_dx1*dO3_dz2)
        Hess(LHR(2,iAtom,1,jAtom)) = Hess(LHR(2,iAtom,1,jAtom))+f_const*(dO1_dy1*dO1_dx2+dO2_dy1*dO2_dx2+dO3_dy1*dO3_dx2)
        Hess(LHR(2,iAtom,2,jAtom)) = Hess(LHR(2,iAtom,2,jAtom))+f_const*(dO1_dy1*dO1_dy2+dO2_dy1*dO2_dy2+dO3_dy1*dO3_dy2)
        Hess(LHR(2,iAtom,3,jAtom)) = Hess(LHR(2,iAtom,3,jAtom))+f_const*(dO1_dy1*dO1_dz2+dO2_dy1*dO2_dz2+dO3_dy1*dO3_dz2)
        Hess(LHR(3,iAtom,1,jAtom)) = Hess(LHR(3,iAtom,1,jAtom))+f_const*(dO1_dz1*dO1_dx2+dO2_dz1*dO2_dx2+dO3_dz1*dO3_dx2)
        Hess(LHR(3,iAtom,2,jAtom)) = Hess(LHR(3,iAtom,2,jAtom))+f_const*(dO1_dz1*dO1_dy2+dO2_dz1*dO2_dy2+dO3_dz1*dO3_dy2)
        Hess(LHR(3,iAtom,3,jAtom)) = Hess(LHR(3,iAtom,3,jAtom))+f_const*(dO1_dz1*dO1_dz2+dO2_dz1*dO2_dz2+dO3_dz1*dO3_dz2)

      end do

      Hess(LHR(1,iAtom,1,iAtom)) = Hess(LHR(1,iAtom,1,iAtom))+f_const*(dO1_dx1*dO1_dx1+dO2_dx1*dO2_dx1+dO3_dx1*dO3_dx1)
      Hess(LHR(2,iAtom,1,iAtom)) = Hess(LHR(2,iAtom,1,iAtom))+f_const*(dO1_dy1*dO1_dx1+dO2_dy1*dO2_dx1+dO3_dy1*dO3_dx1)
      Hess(LHR(2,iAtom,2,iAtom)) = Hess(LHR(2,iAtom,2,iAtom))+f_const*(dO1_dy1*dO1_dy1+dO2_dy1*dO2_dy1+dO3_dy1*dO3_dy1)
      Hess(LHR(3,iAtom,1,iAtom)) = Hess(LHR(3,iAtom,1,iAtom))+f_const*(dO1_dz1*dO1_dx1+dO2_dz1*dO2_dx1+dO3_dz1*dO3_dx1)
      Hess(LHR(3,iAtom,2,iAtom)) = Hess(LHR(3,iAtom,2,iAtom))+f_const*(dO1_dz1*dO1_dy1+dO2_dz1*dO2_dy1+dO3_dz1*dO3_dy1)
      Hess(LHR(3,iAtom,3,iAtom)) = Hess(LHR(3,iAtom,3,iAtom))+f_const*(dO1_dz1*dO1_dz1+dO2_dz1*dO2_dz1+dO3_dz1*dO3_dz1)

    end do
    call mma_deallocate(Grad)
#   ifdef _DEBUGPRINT_
    call TriPrt(' In LNM: Hessian after Rotation','(12f12.7)',Hess,n3)
    call DiagMtrx_T(Hess,n3,iNeg)
#   endif
  end if
end if
call mma_deallocate(xMass)
!                                                                      *
!***********************************************************************
!                                                                      *
! Hessian for tension

do iBond=1,nBonds
  kAtom = iTabBonds(1,iBond)
  lAtom = iTabBonds(2,iBond)
  iBondType = iTabBonds(3,iBond)
  !if (iBondType > Magic_Bond) cycle
  kr = iTabRow(iANr(kAtom))
  lr = iTabRow(iANr(lAtom))
  Help = (kr > 3) .or. (lr > 3)
  xkl = Cart(1,kAtom)-Cart(1,lAtom)
  ykl = Cart(2,kAtom)-Cart(2,lAtom)
  zkl = Cart(3,kAtom)-Cart(3,lAtom)
  rkl2 = xkl**2+ykl**2+zkl**2

  if (ddV_Schlegel .or. Help) then
    Rab = sqrt(rkl2)
    RabCov = CovRad(iANr(kAtom))+CovRad(iANr(lAtom))
    if (((kr == 1) .and. (lr == 1)) .or. Help) then
      gmm = Fact*A_StrH(1)*exp(-A_StrH(2)*(Rab-RabCov))
    else
      ij = iTri(kr,lr)
      gmm = Fact*A_Str/(Rab-B_Str(ij))**3
    end if
  else
    r0 = rAv(kr,lr)
    alpha = aAv(kr,lr)
    gmm = rkr*exp(alpha*(r0**2-rkl2))
    if (btest(iOptC,10)) then
      r0_vdW = r_ref_vdW(kr,lr)
      g_vdW = rkr_vdW*exp(-alpha_vdW*(r0_vdW-sqrt(rkl2))**2)
    else
      g_vdW = Zero
    end if
    gmm = gmm+g_vdW
  end if

  f_const = max(gmm,f_const_Min_)
# ifdef _DEBUGPRINT_
  nqR = nqR+1
  write(u6,*) 'ddV: bonds: kAtom,lAtom=',kAtom,LAtom
  write(u6,*) '          : Bondtype=',iBondType
  !write(u6,*) gmm/rkr, f_const, g_vdW
  write(u6,*) f_const
# endif
  Hxx = f_const*xkl*xkl/rkl2
  Hxy = f_const*xkl*ykl/rkl2
  Hxz = f_const*xkl*zkl/rkl2
  Hyy = f_const*ykl*ykl/rkl2
  Hyz = f_const*ykl*zkl/rkl2
  Hzz = f_const*zkl*zkl/rkl2

  Hess(LHR(1,kAtom,1,kAtom)) = Hess(LHR(1,kAtom,1,kAtom))+Hxx
  Hess(LHR(2,kAtom,1,kAtom)) = Hess(LHR(2,kAtom,1,kAtom))+Hxy
  Hess(LHR(2,kAtom,2,kAtom)) = Hess(LHR(2,kAtom,2,kAtom))+Hyy
  Hess(LHR(3,kAtom,1,kAtom)) = Hess(LHR(3,kAtom,1,kAtom))+Hxz
  Hess(LHR(3,kAtom,2,kAtom)) = Hess(LHR(3,kAtom,2,kAtom))+Hyz
  Hess(LHR(3,kAtom,3,kAtom)) = Hess(LHR(3,kAtom,3,kAtom))+Hzz

  Hess(LHR(1,kAtom,1,lAtom)) = Hess(LHR(1,kAtom,1,lAtom))-Hxx
  Hess(LHR(1,kAtom,2,lAtom)) = Hess(LHR(1,kAtom,2,lAtom))-Hxy
  Hess(LHR(1,kAtom,3,lAtom)) = Hess(LHR(1,kAtom,3,lAtom))-Hxz
  Hess(LHR(2,kAtom,1,lAtom)) = Hess(LHR(2,kAtom,1,lAtom))-Hxy
  Hess(LHR(2,kAtom,2,lAtom)) = Hess(LHR(2,kAtom,2,lAtom))-Hyy
  Hess(LHR(2,kAtom,3,lAtom)) = Hess(LHR(2,kAtom,3,lAtom))-Hyz
  Hess(LHR(3,kAtom,1,lAtom)) = Hess(LHR(3,kAtom,1,lAtom))-Hxz
  Hess(LHR(3,kAtom,2,lAtom)) = Hess(LHR(3,kAtom,2,lAtom))-Hyz
  Hess(LHR(3,kAtom,3,lAtom)) = Hess(LHR(3,kAtom,3,lAtom))-Hzz

  Hess(LHR(1,lAtom,1,lAtom)) = Hess(LHR(1,lAtom,1,lAtom))+Hxx
  Hess(LHR(2,lAtom,1,lAtom)) = Hess(LHR(2,lAtom,1,lAtom))+Hxy
  Hess(LHR(2,lAtom,2,lAtom)) = Hess(LHR(2,lAtom,2,lAtom))+Hyy
  Hess(LHR(3,lAtom,1,lAtom)) = Hess(LHR(3,lAtom,1,lAtom))+Hxz
  Hess(LHR(3,lAtom,2,lAtom)) = Hess(LHR(3,lAtom,2,lAtom))+Hyz
  Hess(LHR(3,lAtom,3,lAtom)) = Hess(LHR(3,lAtom,3,lAtom))+Hzz

end do
#ifdef _DEBUGPRINT_
call TriPrt(' In LNM: Hessian after tension','(12f12.7)',Hess,n3)
call DiagMtrx_T(Hess,n3,iNeg)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Hessian for bending

if (nBonds >= 2) then
  do mAtom=1,mTtAtm
    mr = iTabRow(iANr(mAtom))

    nNeighbor = iTabAtoms(1,0,mAtom)
    if (nNeighbor < 2) cycle
    do iNeighbor=1,nNeighbor
      iAtom = iTabAtoms(1,iNeighbor,mAtom)
      iBond = iTabAtoms(2,iNeighbor,mAtom)
      iBondType = iTabBonds(3,iBond)
      if (iBondType > Magic_Bond) cycle
      ir = iTabRow(iANr(iAtom))

      xmi = Cart(1,iAtom)-Cart(1,mAtom)
      ymi = Cart(2,iAtom)-Cart(2,mAtom)
      zmi = Cart(3,iAtom)-Cart(3,mAtom)
      rmi2 = xmi**2+ymi**2+zmi**2
      rmi = sqrt(rmi2)

      do jNeighbor=1,iNeighbor-1
        jAtom = iTabAtoms(1,jNeighbor,mAtom)
        jBond = iTabAtoms(2,jNeighbor,mAtom)
        jBondType = iTabBonds(3,jBond)
        if (jBondType > Magic_Bond) cycle
        jr = iTabRow(iANr(jAtom))
        Help = (mr > 3) .or. (ir > 3) .or. (jr > 3)

        xmj = Cart(1,jAtom)-Cart(1,mAtom)
        ymj = Cart(2,jAtom)-Cart(2,mAtom)
        zmj = Cart(3,jAtom)-Cart(3,mAtom)
        rmj2 = xmj**2+ymj**2+zmj**2
        rmj = sqrt(rmj2)

        ! Test if zero angle

        Test_zero = xmi*xmj+ymi*ymj+zmi*zmj
        Test_zero = Test_zero/(rmi*rmj)
        if (abs(Test_zero-One) < 1.0e-12_wp) cycle

        xij = Cart(1,jAtom)-Cart(1,iAtom)
        yij = Cart(2,jAtom)-Cart(2,iAtom)
        zij = Cart(3,jAtom)-Cart(3,iAtom)
        rij2 = xij**2+yij**2+zij**2
        rrij = sqrt(rij2)

        if (ddV_Schlegel .or. Help) then
          Rab = rmi
          RabCov = CovRad(iANr(iAtom))+CovRad(iANr(mAtom))
          Rbc = rmj
          RbcCov = CovRad(iANr(jAtom))+CovRad(iANr(mAtom))
          if ((ir == 1) .or. (jr == 1)) then
            gij = Fact*A_Bend(1)
          else
            gij = Fact*A_Bend(2)
          end if
        else
          r0mi = rAv(mr,ir)
          ami = aAv(mr,ir)
          r0mj = rAv(mr,jr)
          amj = aAv(mr,jr)
          gim = exp(ami*(r0mi**2-rmi2))
          gjm = exp(amj*(r0mj**2-rmj2))
          if (btest(iOptC,10)) then
            r0_vdW_im = r_ref_vdW(ir,mr)
            g_vdW_im = exp(-alpha_vdW*(r0_vdW_im-sqrt(rmi2))**2)
            r0_vdW_jm = r_ref_vdW(jr,mr)
            g_vdW_jm = exp(-alpha_vdW*(r0_vdW_jm-sqrt(rmj2))**2)
          else
            g_vdW_im = Zero
            g_vdW_jm = Zero
          end if
          g_vdW_im = g_vdW_im*rkr_vdW/rkr
          g_vdW_jm = g_vdW_jm*rkr_vdW/rkr
          gim = gim+Half*g_vdW_im
          gjm = gjm+Half*g_vdW_jm
          gij = rkf*gim*gjm
        end if
        rL2 = (ymi*zmj-zmi*ymj)**2+(zmi*xmj-xmi*zmj)**2+(xmi*ymj-ymi*xmj)**2
        !hjw modified
        if (rL2 < 1.0e-14_wp) then
          rL = Zero
        else
          rL = sqrt(rL2)
        end if
        gij = max(gij,f_const_Min_)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'iAtom,mAtom,jAtom=',iAtom,mAtom,jAtom
        write(u6,*) 'gij=',gij
        write(u6,*) 'rmj=',rmj
        write(u6,*) 'rmi=',rmi
        write(u6,*) 'rrij=',rrij
#       endif

        if ((rmj > rZero) .and. (rmi > rZero) .and. (rrij > rZero)) then
#         ifdef _DEBUGPRINT_
          nqB = nqB+1
#         endif
          SinPhi = rL/(rmj*rmi)
          rmidotrmj = xmi*xmj+ymi*ymj+zmi*zmj
          CosPhi = rmidotrmj/(rmj*rmi)

          ! Non-linear case

          Thr_Line = sin(25.0_wp*deg2rad)
          if (mTtAtm == 3) Thr_Line = rZero
          if (SinPhi > Thr_Line) then
            si(1) = (xmi/rmi*cosphi-xmj/rmj)/(rmi*sinphi)
            si(2) = (ymi/rmi*cosphi-ymj/rmj)/(rmi*sinphi)
            si(3) = (zmi/rmi*cosphi-zmj/rmj)/(rmi*sinphi)
            sj(1) = (cosphi*xmj/rmj-xmi/rmi)/(rmj*sinphi)
            sj(2) = (cosphi*ymj/rmj-ymi/rmi)/(rmj*sinphi)
            sj(3) = (cosphi*zmj/rmj-zmi/rmi)/(rmj*sinphi)
            sm(1) = -si(1)-sj(1)
            sm(2) = -si(2)-sj(2)
            sm(3) = -si(3)-sj(3)
            do icoor=1,3
              do jCoor=1,3
                if (mAtom > iAtom) then
                  Hess(LHR(icoor,mAtom,jcoor,iAtom)) = Hess(LHR(icoor,mAtom,jcoor,iAtom))+gij*sm(icoor)*si(jcoor)
                else
                  Hess(LHR(icoor,iAtom,jcoor,mAtom)) = Hess(LHR(icoor,iAtom,jcoor,mAtom))+gij*si(icoor)*sm(jcoor)
                end if
                if (mAtom > jAtom) then
                  Hess(LHR(icoor,mAtom,jcoor,jAtom)) = Hess(LHR(icoor,mAtom,jcoor,jAtom))+gij*sm(icoor)*sj(jcoor)
                else
                  Hess(LHR(icoor,jAtom,jcoor,mAtom)) = Hess(LHR(icoor,jAtom,jcoor,mAtom))+gij*sj(icoor)*sm(jcoor)
                end if
                if (iAtom > jAtom) then
                  Hess(LHR(icoor,iAtom,jcoor,jAtom)) = Hess(LHR(icoor,iAtom,jcoor,jAtom))+gij*si(icoor)*sj(jcoor)
                else
                  Hess(LHR(icoor,jAtom,jcoor,iAtom)) = Hess(LHR(icoor,jAtom,jcoor,iAtom))+gij*sj(icoor)*si(jcoor)
                end if
              end do
            end do
            do icoor=1,3
              do jCoor=1,icoor
                Hess(LHR(icoor,iAtom,jcoor,iAtom)) = Hess(LHR(icoor,iAtom,jcoor,iAtom))+gij*si(icoor)*si(jcoor)
                Hess(LHR(icoor,mAtom,jcoor,mAtom)) = Hess(LHR(icoor,mAtom,jcoor,mAtom))+gij*sm(icoor)*sm(jcoor)
                Hess(LHR(icoor,jAtom,jcoor,jAtom)) = Hess(LHR(icoor,jAtom,jcoor,jAtom))+gij*sj(icoor)*sj(jcoor)

              end do
            end do
          else

            ! Linear case

            if ((abs(ymi) > rZero) .or. (abs(xmi) > rZero)) then
              x(1) = -ymi
              y(1) = xmi
              z(1) = Zero
              x(2) = -xmi*zmi
              y(2) = -ymi*zmi
              z(2) = xmi*xmi+ymi*ymi
            else
              x(1) = One
              y(1) = Zero
              z(1) = Zero
              x(2) = Zero
              y(2) = One
              z(2) = Zero
            end if
#           ifdef _DEBUGPRINT_
            nqB = nqB+2
#           endif
            do i=1,2
              r1 = sqrt(x(i)**2+y(i)**2+z(i)**2)
              cosThetax = x(i)/r1
              cosThetay = y(i)/r1
              cosThetaz = z(i)/r1
              si(1) = -cosThetax/rmi
              si(2) = -cosThetay/rmi
              si(3) = -cosThetaz/rmi
              sj(1) = -cosThetax/rmj
              sj(2) = -cosThetay/rmj
              sj(3) = -cosThetaz/rmj
              sm(1) = -(si(1)+sj(1))
              sm(2) = -(si(2)+sj(2))
              sm(3) = -(si(3)+sj(3))

              do icoor=1,3
                do jCoor=1,3
                  if (mAtom > iAtom) then
                    Hess(LHR(icoor,mAtom,jcoor,iAtom)) = Hess(LHR(icoor,mAtom,jcoor,iAtom))+gij*sm(icoor)*si(jcoor)
                  else
                    Hess(LHR(icoor,iAtom,jcoor,mAtom)) = Hess(LHR(icoor,iAtom,jcoor,mAtom))+gij*si(icoor)*sm(jcoor)
                  end if
                  if (mAtom > jAtom) then
                    Hess(LHR(icoor,mAtom,jcoor,jAtom)) = Hess(LHR(icoor,mAtom,jcoor,jAtom))+gij*sm(icoor)*sj(jcoor)
                  else
                    Hess(LHR(icoor,jAtom,jcoor,mAtom)) = Hess(LHR(icoor,jAtom,jcoor,mAtom))+gij*sj(icoor)*sm(jcoor)
                  end if
                  if (iAtom > jAtom) then
                    Hess(LHR(icoor,iAtom,jcoor,jAtom)) = Hess(LHR(icoor,iAtom,jcoor,jAtom))+gij*si(icoor)*sj(jcoor)
                  else
                    Hess(LHR(icoor,jAtom,jcoor,iAtom)) = Hess(LHR(icoor,jAtom,jcoor,iAtom))+gij*sj(icoor)*si(jcoor)
                  end if
                end do
              end do
              do icoor=1,3
                do jCoor=1,icoor
                  Hess(LHR(icoor,iAtom,jcoor,iAtom)) = Hess(LHR(icoor,iAtom,jcoor,iAtom))+gij*si(icoor)*si(jcoor)
                  Hess(LHR(icoor,mAtom,jcoor,mAtom)) = Hess(LHR(icoor,mAtom,jcoor,mAtom))+gij*sm(icoor)*sm(jcoor)
                  Hess(LHR(icoor,jAtom,jcoor,jAtom)) = Hess(LHR(icoor,jAtom,jcoor,jAtom))+gij*sj(icoor)*sj(jcoor)
                end do
              end do
            end do
          end if
        end if

      end do
    end do
  end do
# ifdef _DEBUGPRINT_
  call TriPrt(' In LNM: Hessian after bending','(12f12.7)',Hess,n3)
  call DiagMtrx_T(Hess,n3,iNeg)
# endif
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Hessian for torsion

if (nBonds >= 3) then
  do iBond=1,nBonds
    jAtom = iTabBonds(1,iBond)
    kAtom = iTabBonds(2,iBond)
    iBondType = iTabBonds(3,iBond)
    Fact = One
    if (iBondType > Magic_Bond) Fact = Two
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) '*',jAtom,kAtom,' *'
    write(u6,*)
    write(u6,*) 'BondType=',iBondType
#   endif

    ! Allow center bond to be a "magic" bond

    !if (iBondType == vdW_Bond) cycle

    jr = iTabRow(iANr(jAtom))
    kr = iTabRow(iANr(kAtom))

    xyz(:,2) = Cart(:,jAtom)
    xyz(:,3) = Cart(:,kAtom)

    nNeighbor_j = iTabAtoms(1,0,jAtom)
    if (nNeighbor_j < 2) cycle
    nNeighbor_k = iTabAtoms(1,0,kAtom)
    if (nNeighbor_k < 2) cycle

    do jNeighbor=1,nNeighbor_j
      iAtom = iTabAtoms(1,jNeighbor,jAtom)
#     ifdef _DEBUGPRINT_
      !write(u6,*)
      !write(u6,*) iAtom,jAtom,kAtom,' *'
      !write(u6,*)
#     endif
      jBond = iTabAtoms(2,jNeighbor,jAtom)
      if (iBond == jBond) cycle
      jBondType = iTabBonds(3,jBond)
      !if (jBondType == vdW_Bond) cycle
      if (jBondType > Magic_Bond) cycle
      ir = iTabRow(iANr(iAtom))

      xyz(:,1) = Cart(:,iAtom)

      do kNeighbor=1,nNeighbor_k
        lAtom = iTabAtoms(1,kNeighbor,kAtom)
        kBond = iTabAtoms(2,kNeighbor,kAtom)
        if (iBond == kBond) cycle
        if (lAtom == iAtom) cycle
        kBondType = iTabBonds(3,kBond)
        !if (kBondType == vdW_Bond) cycle
        if (kBondType > Magic_Bond) cycle
        lr = iTabRow(iANr(lAtom))
        Help = (kr > 3) .or. (ir > 3) .or. (jr > 3) .or. (lr > 3)

        xyz(:,4) = Cart(:,lAtom)

        rij(1) = Cart(1,iAtom)-Cart(1,jAtom)
        rij(2) = Cart(2,iAtom)-Cart(2,jAtom)
        rij(3) = Cart(3,iAtom)-Cart(3,jAtom)
        rij2 = rij(1)**2+rij(2)**2+rij(3)**2

        rjk(1) = Cart(1,jAtom)-Cart(1,kAtom)
        rjk(2) = Cart(2,jAtom)-Cart(2,kAtom)
        rjk(3) = Cart(3,jAtom)-Cart(3,kAtom)
        rjk2 = rjk(1)**2+rjk(2)**2+rjk(3)**2

        rkl(1) = Cart(1,kAtom)-Cart(1,lAtom)
        rkl(2) = Cart(2,kAtom)-Cart(2,lAtom)
        rkl(3) = Cart(3,kAtom)-Cart(3,lAtom)
        rkl2 = rkl(1)**2+rkl(2)**2+rkl(3)**2

        ! Allow only angles in the range of 35-145
        A35 = 35.0_wp*deg2rad
        CosFi_Max = cos(A35)
        CosFi2 = (rij(1)*rjk(1)+rij(2)*rjk(2)+rij(3)*rjk(3))/sqrt(rij2*rjk2)
        if (abs(CosFi2) > CosFi_Max) cycle
        CosFi3 = (rkl(1)*rjk(1)+rkl(2)*rjk(2)+rkl(3)*rjk(3))/sqrt(rkl2*rjk2)
        if (abs(CosFi3) > CosFi_Max) cycle
#       ifdef _DEBUGPRINT_
        write(u6,*) 'CosFi2,CosFi3=',CosFi2,CosFi3
        write(u6,*) 'rij=',rij,rij2
        write(u6,*) 'rjk=',rjk,rjk2
        write(u6,*) 'rkl=',rkl,rkl2
#       endif

        if (ddV_Schlegel .or. Help) then
          Rab = sqrt(rij2)
          RabCov = (CovRadT(iANr(iAtom))+CovRadT(iANr(jAtom)))/Angstrom
          Rbc = sqrt(rjk2)/Fact
          RbcCov = (CovRadT(iANr(jAtom))+CovRadT(iANr(kAtom)))/Angstrom
          Diff = RbcCov-Rbc
          if (Diff < Zero) Diff = Zero
          tij = Fact*A_Trsn(1)+A_Trsn(2)*Diff
        else
          rij0 = rAv(ir,jr)**2
          aij = aAv(ir,jr)
          rjk0 = rAv(jr,kr)**2
          ajk = aAv(jr,kr)
          rkl0 = rAv(kr,lr)**2
          akl = aAv(kr,lr)
          ! Magic bond fix
          rjk2 = rjk2/Fact**2

          g_ij = exp(aij*(rij0-rij2))
          g_jk = exp(ajk*(rjk0-rjk2))
          g_kl = exp(akl*(rkl0-rkl2))
          if (btest(iOptC,10)) then
            r0_vdW_ij = r_ref_vdW(ir,jr)
            g_vdW_ij = exp(-alpha_vdW*(r0_vdW_ij-sqrt(rij2))**2)
            r0_vdW_jk = r_ref_vdW(jr,kr)
            g_vdW_jk = exp(-alpha_vdW*(r0_vdW_jk-sqrt(rjk2))**2)
            r0_vdW_kl = r_ref_vdW(kr,lr)
            g_vdW_kl = exp(-alpha_vdW*(r0_vdW_kl-sqrt(rkl2))**2)
          else
            g_vdW_ij = Zero
            g_vdW_jk = Zero
            g_vdW_kl = Zero
          end if
          g_vdW_ij = g_vdW_ij*rkr_vdW/rkr
          g_vdW_jk = g_vdW_jk*rkr_vdW/rkr
          g_vdW_kl = g_vdW_kl*rkr_vdW/rkr
          g_ij = g_ij+Half*g_vdW_ij
          g_jk = g_jk+Half*g_vdW_jk
          g_kl = g_kl+Half*g_vdW_kl
          tij = rkt*g_ij*g_jk*g_kl
        end if
        tij = max(tij,f_const_Min_)
        if (Torsion_Check(iAtom,jAtom,kAtom,lAtom,xyz,iTabAtoms,nMax,mTtAtm)) tij = max(tij,Ten*f_const_Min_)
#       ifdef _DEBUGPRINT_
        nqT = nqT+1
        write(u6,*)
        write(u6,*) iAtom,jAtom,kAtom,lAtom
        write(u6,*) tij
#       endif

        call Trsn(xyz,4,Tau,C,.false.,.false.,'        ',Dum,.false.)
        si(:) = C(:,1)
        sj(:) = C(:,2)
        sk(:) = C(:,3)
        sl(:) = C(:,4)
#       ifdef _DEBUGPRINT_
        !call RecPrt('C',' ',C,3,4)
#       endif

        ! Off diagonal block

        do icoor=1,3
          do jCoor=1,3
            Hess(LHR(icoor,iAtom,jcoor,jAtom)) = Hess(LHR(icoor,iAtom,jcoor,jAtom))+tij*si(icoor)*sj(jcoor)
            Hess(LHR(icoor,iAtom,jcoor,kAtom)) = Hess(LHR(icoor,iAtom,jcoor,kAtom))+tij*si(icoor)*sk(jcoor)
            Hess(LHR(icoor,iAtom,jcoor,lAtom)) = Hess(LHR(icoor,iAtom,jcoor,lAtom))+tij*si(icoor)*sl(jcoor)
            Hess(LHR(icoor,jAtom,jcoor,kAtom)) = Hess(LHR(icoor,jAtom,jcoor,kAtom))+tij*sj(icoor)*sk(jcoor)
            Hess(LHR(icoor,jAtom,jcoor,lAtom)) = Hess(LHR(icoor,jAtom,jcoor,lAtom))+tij*sj(icoor)*sl(jcoor)
            Hess(LHR(icoor,kAtom,jcoor,lAtom)) = Hess(LHR(icoor,kAtom,jcoor,lAtom))+tij*sk(icoor)*sl(jcoor)

          end do
        end do

        ! Diagonal block

        do icoor=1,3
          do jCoor=1,icoor
            Hess(LHR(icoor,iAtom,jcoor,iAtom)) = Hess(LHR(icoor,iAtom,jcoor,iAtom))+tij*si(icoor)*si(jcoor)
            Hess(LHR(icoor,jAtom,jcoor,jAtom)) = Hess(LHR(icoor,jAtom,jcoor,jAtom))+tij*sj(icoor)*sj(jcoor)
            Hess(LHR(icoor,kAtom,jcoor,kAtom)) = Hess(LHR(icoor,kAtom,jcoor,kAtom))+tij*sk(icoor)*sk(jcoor)
            Hess(LHR(icoor,lAtom,jcoor,lAtom)) = Hess(LHR(icoor,lAtom,jcoor,lAtom))+tij*sl(icoor)*sl(jcoor)

          end do
        end do

      end do          ! iNeigbor_k
    end do            ! iNeighbor_j
  end do              ! iBonds
# ifdef _DEBUGPRINT_
  call TriPrt(' In LNM: Hessian after torsion','(12f12.7)',Hess,n3)
  call DiagMtrx_T(Hess,n3,iNeg)
# endif
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!                                  k
!                                 /
! Hessian for out-of-plane   j - i
!                                 \
!                                  l

!if (.false.) then
if (nBonds >= 3) then

  do iAtom=1,mTtAtm

    nNeighbor_i = iTabAtoms(1,0,iAtom)
    !write(u6,*) 'iAtom,nNeighbor_i=',iAtom,nNeighbor_i
    if (nNeighbor_i < 3) cycle
    ir = iTabRow(iANr(iAtom))
    xyz(:,4) = Cart(:,iAtom)

    do iNb0=1,nNeighbor_i
      jAtom = iTabAtoms(1,iNb0,iAtom)
      !write(u6,*) 'jAtom=',jAtom
      jr = iTabRow(iANr(jAtom))
      iBond = iTabAtoms(2,iNb0,iAtom)
      iBondType = iTabBonds(3,iBond)
      !write(u6,*) 'iBondType=',iBondType
      nCoBond_j = nCoBond(jAtom,mTtAtm,nMax,iTabBonds,nBonds,iTabAtoms)
      if (nCoBond_j > 1) cycle
      !if (iBondType == vdW_Bond) cycle
      if (iBondType > Magic_Bond) cycle
      xyz(:,1) = Cart(:,jAtom)

      do iNb1=1,nNeighbor_i
        kAtom = iTabAtoms(1,iNb1,iAtom)
        !write(u6,*) 'kAtom=',kAtom
        kBond = iTabAtoms(2,iNb1,iAtom)
        if (kAtom == jAtom) cycle
        kBondType = iTabBonds(3,kBond)
        !write(u6,*) 'kBondType=',kBondType
        !if (kBondType == vdW_Bond) cycle
        if (kBondType > Magic_Bond) cycle
        kr = iTabRow(iANr(kAtom))

        xyz(:,2) = Cart(:,kAtom)

        do iNb2=1,nNeighbor_i
          lAtom = iTabAtoms(1,iNb2,iAtom)
          !write(u6,*) 'lAtom=',lAtom
          lBond = iTabAtoms(2,iNb2,iAtom)

          if (lAtom == jAtom) cycle
          if (lAtom <= kAtom) cycle
          lBondType = iTabBonds(3,lBond)
          !write(u6,*) 'lBondType=',lBondType
          !if (lBondType == vdW_Bond) cycle
          if (lBondType > Magic_Bond) cycle
          lr = iTabRow(iANr(lAtom))
          Help = (kr > 3) .or. (ir > 3) .or. (jr > 3) .or. (lr > 3)

          !Write(u6,*) 'i,j,k,l=',iAtom,jAtom,kAtom,lAtom
          !Write(u6,*) 'Help=',Help

          xyz(:,3) = Cart(:,lAtom)

          rij(:) = Cart(:,iAtom)-Cart(:,jAtom)

          rik(:) = Cart(:,iAtom)-Cart(:,kAtom)

          ril(:) = Cart(:,iAtom)-Cart(:,lAtom)

          rij2 = rij(1)**2+rij(2)**2+rij(3)**2
          rik2 = rik(1)**2+rik(2)**2+rik(3)**2
          ril2 = ril(1)**2+ril(2)**2+ril(3)**2

          ThrFi1 = cos(90.0_wp*deg2rad)
          ThrFi2 = cos(150.0_wp*deg2rad)
          CosFi2 = (rij(1)*rik(1)+rij(2)*rik(2)+rij(3)*rik(3))/sqrt(rij2*rik2)
          if ((CosFi2 > ThrFi1) .or. (CosFi2 < ThrFi2)) cycle

          CosFi3 = (rij(1)*ril(1)+rij(2)*ril(2)+rij(3)*ril(3))/sqrt(rij2*ril2)
          if ((CosFi3 > ThrFi1) .or. (CosFi3 < ThrFi2)) cycle

          CosFi4 = (rik(1)*ril(1)+rik(2)*ril(2)+rik(3)*ril(3))/sqrt(rik2*ril2)
          if ((CosFi4 > ThrFi1) .or. (CosFi4 < ThrFi2)) cycle
#         ifdef _DEBUGPRINT_
          write(u6,*) 'CosFi2,CosFi3,CosFi4=',CosFi2,CosFi3,CosFi4
#         endif

          if (ddV_Schlegel .or. Help) then

            ! I do not have a clue to how this will really work!

            tij = f_const_Min_
          else
            rij0 = rAv(ir,jr)**2
            aij = aAv(ir,jr)
            rik0 = rAv(ir,kr)**2
            aik = aAv(ir,kr)
            ril0 = rAv(ir,lr)**2
            ail = aAv(ir,lr)
            beta = rko*exp((aij*rij0+aik*rik0+ail*ril0))
            tij = beta*exp(-(aij*rij2+aik*rik2+ail*ril2))
          end if
          !tij = max(tij,f_const_Min_)

          call OutofP(xyz,4,Tau,C,.false.,.false.,'        ',Dum,.false.)
          if (abs(Tau) > 25.0_wp*deg2rad) cycle
#         ifdef _DEBUGPRINT_
          nqO = nqO+1
#         endif

          si(:) = C(:,4)
          sj(:) = C(:,1)
          sk(:) = C(:,2)
          sl(:) = C(:,3)
#         ifdef _DEBUGPRINT_
          write(u6,*) 'iAtoms=',iAtom,jAtom,kAtom,lAtom
          write(u6,*) 'tij,Tau=',tij,Tau
          !call RecPrt('si',' ',si,1,3)
          !call RecPrt('sj',' ',sj,1,3)
          !call RecPrt('sk',' ',sk,1,3)
          !call RecPrt('sl',' ',sl,1,3)
#         endif

          ! Off diagonal block

          do icoor=1,3
            do jCoor=1,3
              Hess(LHR(icoor,iAtom,jcoor,jAtom)) = Hess(LHR(icoor,iAtom,jcoor,jAtom))+tij*si(icoor)*sj(jcoor)
              Hess(LHR(icoor,iAtom,jcoor,kAtom)) = Hess(LHR(icoor,iAtom,jcoor,kAtom))+tij*si(icoor)*sk(jcoor)
              Hess(LHR(icoor,iAtom,jcoor,lAtom)) = Hess(LHR(icoor,iAtom,jcoor,lAtom))+tij*si(icoor)*sl(jcoor)
              Hess(LHR(icoor,jAtom,jcoor,kAtom)) = Hess(LHR(icoor,jAtom,jcoor,kAtom))+tij*sj(icoor)*sk(jcoor)
              Hess(LHR(icoor,jAtom,jcoor,lAtom)) = Hess(LHR(icoor,jAtom,jcoor,lAtom))+tij*sj(icoor)*sl(jcoor)
              Hess(LHR(icoor,kAtom,jcoor,lAtom)) = Hess(LHR(icoor,kAtom,jcoor,lAtom))+tij*sk(icoor)*sl(jcoor)

            end do
          end do

          ! Diagonal block

          do icoor=1,3
            do jCoor=1,icoor
              Hess(LHR(icoor,iAtom,jcoor,iAtom)) = Hess(LHR(icoor,iAtom,jcoor,iAtom))+tij*si(icoor)*si(jcoor)
              Hess(LHR(icoor,jAtom,jcoor,jAtom)) = Hess(LHR(icoor,jAtom,jcoor,jAtom))+tij*sj(icoor)*sj(jcoor)
              Hess(LHR(icoor,kAtom,jcoor,kAtom)) = Hess(LHR(icoor,kAtom,jcoor,kAtom))+tij*sk(icoor)*sk(jcoor)
              Hess(LHR(icoor,lAtom,jcoor,lAtom)) = Hess(LHR(icoor,lAtom,jcoor,lAtom))+tij*sl(icoor)*sl(jcoor)

            end do
          end do

        end do        ! iNb2
      end do          ! iNb1
    end do            ! iCase
  end do              ! iBond
# ifdef _DEBUGPRINT_
  call TriPrt(' In LNM: Hessian after out-of-plane','(12f12.7)',Hess,n3)
  call DiagMtrx_T(Hess,n3,iNeg)
# endif
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*) 'ddV: nqR, nqB, nqT, nqO=',nqR,nqB,nqT,nqO
#endif

return

end subroutine ddV_
