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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine sg_two_pdm_full(SGS,CIS,EXS,Psi,nCSFs,PsiSym,P,NLEV)

use Symmetry_Info, only: Mul
use sguga, only: CIStruct, EXStruct, sg_epq_psi, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) :: nCSFs, PsiSym, nLev
real(kind=wp), intent(in) :: Psi(nCSFs)
real(kind=wp), intent(inout) :: P(nLev,nLev,nLev,nLev)
integer(kind=iwp) :: ijSym, iOrb, iSym, jOrb, jSym, klSym, kOrb, kSym, lOrb, lSym, MaxDim, mCSFs
real(kind=wp) :: D_ij, P_klij
real(kind=wp), pointer :: Eij_Psi(:), Elk_Psi(:)
real(kind=wp), allocatable, target :: Eij_Psi_X(:), Elk_Psi_X(:)
real(kind=wp), parameter :: CPQ = One

MaxDim = maxval(CIS%nCSF(:))
call mma_allocate(Eij_Psi_X,MaxDim,Label='Eij_Psi_X')
call mma_allocate(Elk_Psi_X,MaxDim,Label='Elk_Psi_X')

P(:,:,:,:) = Zero
D_ij = Zero

do iOrb=1,nLev
  iSym = SGS%ISM(iOrb)
  do jOrb=1,nLev
    jSym = SGS%ISM(jOrb)
    ijSym = Mul(iSym,jSym)
    mCSFs = CIS%nCSF(Mul(PsiSym,ijSym))

    Eij_Psi(1:mCSFs) => Eij_Psi_X
    Elk_Psi(1:mCSFs) => Elk_Psi_X

    ! Operate with E_ij on |Psi> and produce E_ij|Psi>
    ! Note, iOrb>=jOrb

    ! Compute Dij to be distributed below

    Eij_Psi(:) = Zero
    call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,PsiSym,Psi,Eij_Psi)
    if (iSym == jSym) D_ij = dot_product(Psi,Eij_Psi)

    ! Note, in the case of a RASSCF the resulting sigma vector is incomplete. For the case of E_ij
    ! the complete expansion of E_ij|Psi> would require CSFs that are outside the set of CSFs defining the
    ! RAS expansion. However, if we cap this sigma vector with vectors that have zero coefficient for these
    ! external CSFs we experience no errors.

    do kOrb=1,nLev
      kSym = SGS%ISM(kOrb)
      lOrb = kOrb
      lSym = SGS%ISM(lOrb)
      klSym = Mul(kSym,lSym)
      if (ijSym /= klSym) cycle

      ! Add the -d_il E_jk term, here in the form of a -d_kl E_ij term of Piklj
      P(iOrb,kOrb,lOrb,jOrb) = P(iOrb,kOrb,lOrb,jOrb)-D_ij
    end do

    ! kOrb>=lOrb
    do kOrb=1,nLev
      kSym = SGS%ISM(kOrb)
      do lOrb=1,nLev
        lSym = SGS%ISM(lOrb)
        klSym = Mul(kSym,lSym)
        if (ijSym /= klSym) cycle

        ! Operate with E_lk on |Psi> and produce E_lk|Psi>. Note, since k>=l E_lk|Psi>
        ! is completely definded by the CSFs of the RASSCF space.

        Elk_Psi(1:mCSFs) = Zero
        ! <Psi|E_kl is computed as E_lk|Psi>
        call SG_Epq_Psi(SGS,CIS,EXS,lOrb,kOrb,CPQ,PsiSym,Psi,Elk_Psi)
        ! cap with <Psi|E_kl on E_ij|Psi>, contribution to Pklij
        P_klij = dot_product(Elk_Psi,Eij_Psi)

        P(kOrb,lOrb,iOrb,jOrb) = P(kOrb,lOrb,iOrb,jOrb)+P_klij

      end do
    end do
    nullify(Eij_Psi,Elk_Psi)

  end do
end do

call mma_deallocate(Elk_Psi_X)
call mma_deallocate(Eij_Psi_X)

P(:,:,:,:) = Half*P(:,:,:,:)
call RecPrt('P full',' ',P,nLev**2,nLev**2)

end subroutine sg_two_pdm_full
