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
! Copyright (C) 2021,2026, Yoshio Nishimoto                            *
!***********************************************************************

subroutine DEPSAOff(nConf,nState,nAshT,nOLag,NBSQT,CLag,OLag,DEPSA,DEPSA_can,CDensF,FIFA)

use caspt2_global, only: IDTCEX, LUCIEX
use general_data, only: NASH
use caspt2_module, only: ISCF, NAES, NBAS, NDEL, NFRO, NISH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nConf, nState, nAshT, nOLag, NBSQT
real(kind=wp), intent(in) :: CLag(nConf,nState), OLag(nOLag), DEPSA_can(nAshT,nAshT), CDensF(nAshT,nAshT), FIFA(NBSQT)
real(kind=wp), intent(inout) :: DEPSA(nAshT,nAshT)
integer(kind=iwp) :: ID, iMO, iState, iSym, iT, iTabs, iTtot, iU, iUabs, iUtot, nAshI, nBasI, nFroI, nIshI, nOrbI
real(kind=wp) :: EigT, EigU, Tmp
real(kind=wp), allocatable :: CDens(:,:), CI(:), G1(:,:), SGM(:)
real(kind=wp), parameter :: ThrDeg = 1.0e-8_wp

! See J. Chem. Phys. 2023, 158, 174112. for details, in particular,
! Section II C 4 "Non-invariance with respect to active MOs"
! Off-diagonal active density, Eq. (70):
! D_{tu} += -(1/2)*[(OLag_tu - OLag_ut) + X_tu ] / (f_t-f_u)

! The original code first solved the linear equation: |CLag_I> = H^{CC}|Q_I> = 2w_I (H-E_I)|Q_I> for Q,
! where H^{CC} is CI Hessian and |Q> is the solution.
! Then, contracted with the CI-active rotation Hessian:
! X_{tu} = Q * H_{tu}^{Ca} = 2 sum_I w_I <Q_I | [E_{tu}^-,H] |I>, where E_{tu}^- = E_{tu} - E_{ut}
! The active orbitals and CI vectors rotated together in the opposite direction keep the reference states, so
! H_{tu}^{Ca} = -H^{CC}*r_{tu}, with r_{tu,I} = E_{tu}^-|I>
! Then X_{tu} can be directly computed: X_{tu} = -Q*H^{CC}*r_{tu} = -CLag*r_{tu} = -sum_I <CLag_I|E_{tu}^-|I>
! CDens(t,u) = sum_I <CLag_I|E_ut|I> with the CI vectors in the quasi-canonical basis, so X_tu = CDens(t,u)-CDens(u,t).
! The rotation of a single CSF is outside the CI space.

call mma_allocate(CDens,nAshT,nAshT,Label='CDens')
! CDensF is the part of the frozen orbitals, which depend on the SCF density D.
! With F = DRDMSCF, |CLagF_I> = w_I sum_xy F(x,y) Exy |PsiI>, and
! X_tu = -sum_I <CLagF_I|Etu-Eut|PsiI> = -sum_I w_I sum_xy F(x,y) <PsiI|Exy*(Etu-Eut)|PsiI>
! F is symmetric, so
!   sum_xy F(x,y)*<PsiI|Exy*(Etu-Eut)|PsiI>
!    = 0.5*sum_xy F(x,y)*<PsiI|[Exy,Etu-Eut]|PsiI>
!    = 0.5*sum_xy F(x,y)*<PsiI|Exu*del(ty)-Ety*del(xu)-Ext*del(uy)+Euy*del(xt)|PsiI>
!    = 0.5*<PsiI|Fxt*Exu-Fuy*Ety-Fxu*Ext+Fty*Euy|PsiI>
! and X_tu = -0.5*(Fxt*Dxu-Fuy*Dty-Fxu*Dxt+Fty*Duy) = (DF-FD)_tu.
! CDensF = D*F, so CDens(t,u)-CDens(u,t) gives frozen orbital contributions
CDens(:,:) = CDensF(:,:)
if (ISCF == 0) then
  call mma_allocate(CI,nConf,Label='CI')
  call mma_allocate(SGM,nConf,Label='SGM')
  call mma_allocate(G1,nAshT,nAshT,Label='G1')
  do iState=1,nState
    ID = IDTCEX(iState)
    call DDaFile(LUCIEX,2,CI,nConf,ID)
    call DENS1T_RPT2(CI,CLag(:,iState),SGM,G1,nAshT)
    CDens(:,:) = CDens(:,:)+G1(:,:)
  end do
  call mma_deallocate(CI)
  call mma_deallocate(SGM)
  call mma_deallocate(G1)
end if

! The orbital part is much easier; similar to the frozen core approximation.
! Corresponds to the first term in Eq. (70)

iMO = 1
do iSym=1,nSym
  nAshI = nAsh(iSym)
  nOrbI = nBas(iSym)-nDel(iSym)
  if (nAshI /= 0) then
    nFroI = nFro(iSym)
    nIshI = nIsh(iSym)
    nBasI = nBas(iSym)
    do iT=1,nAshI
      iTabs = iT+nAES(iSym)
      iTtot = iT+nFroI+nIshI
      EigT = FIFA(iMO+iTtot-1+nBasI*(iTtot-1))
      do iU=1,iT-1
        iUabs = iU+nAES(iSym)
        iUtot = iU+nFroI+nIshI
        EigU = FIFA(iMO+iUtot-1+nBasI*(iUtot-1))
        if (abs(EigT-EigU) < ThrDeg) then
          Tmp = DEPSA_can(iTabs,iUabs) ! avoid singular
        else
          Tmp = OLag(iMO+iTtot-1+nOrbI*(iUtot-1))-OLag(iMO+iUtot-1+nOrbI*(iTtot-1))+CDens(iTabs,iUabs)-CDens(iUabs,iTabs)
          Tmp = -Tmp/(EigT-EigU)*Half
        end if

        DEPSA(iTabs,iUabs) = DEPSA(iTabs,iUabs)+Tmp
        DEPSA(iUabs,iTabs) = DEPSA(iUabs,iTabs)+Tmp
      end do
    end do
  end if
  iMO = iMO+nOrbI*nOrbI
end do

call mma_deallocate(CDens)

end subroutine DEPSAOff
