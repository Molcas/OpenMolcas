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
subroutine sg_two_pdm(SGS,CIS,EXS,Psi,nCSFs,PsiSym,P,PA,nP)

use Index_functions, only: iTri
use Symmetry_Info, only: Mul
use sguga, only: CIStruct, EXStruct, sg_epq_psi, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) :: nCSFs, PsiSym, nP
real(kind=wp), intent(in) :: Psi(nCSFs)
real(kind=wp), intent(inout) :: P(nP), PA(nP)
integer(kind=iwp) :: ijOrb, ijSym, ikljOrb, ikOrb, iOrb, iSym, jOrb, jSym, klijOrb, klOrb, klSym, kOrb, kSym, ljOrb, lOrb, &
                     lOrb_Min, lSym, MaxDim, mCSFs
real(kind=wp) :: D_ij, P_klij
real(kind=wp), pointer :: Eij_Psi(:), Eji_Psi(:), Elk_Psi(:)
real(kind=wp), allocatable, target :: Eij_Psi_X(:), Elk_Psi_X(:), Eji_Psi_X(:)
real(kind=wp), parameter :: CPQ = One

MaxDim = maxval(CIS%nCSF(:))
call mma_allocate(Eij_Psi_X,MaxDim,Label='Eij_Psi_X')
call mma_allocate(Eji_Psi_X,MaxDim,Label='Eji_Psi_X')
call mma_allocate(Elk_Psi_X,MaxDim,Label='Elk_Psi_X')

P(:) = Zero
PA(:) = Zero
D_ij = Zero

! Compute P_lk,ij = <Psi|e_lkij|Psi> + <Psi|e_lkji|Psi>
! Where e_lkij = E_lk E_ij - delta_ki E_lj
! Hence,
!
! P_lk,ij = sum_m <Psi|E_lk|m>(<m|Eij|Psi>+<m|E_ji|Psi>) - delta_ki D_lj -delta_kj D_li
!
! In the symmetrization set we use the fact that for none zero elements we have that
! <m|E_ij|m'>=<m'|E_ji|m>
! This is trivial in the case of no symmetry. However, in the case of symmetry we have to be careful.
!
! Let say that |m> belongs to the set of CSFs that are included in |Psi>, while |m'> doesn't. Then
! there is no symmetrization to be achieved. We need to keep track on this in Mk_Eij_Psi.
!
! The size of the outer loop is nOrb*(nOrb+1)/2, i.e. nOrb=10 => 55 tasks to parallelize over.

do iOrb=1,SGS%nLev
  iSym = SGS%ISM(iOrb)
  do jOrb=1,iOrb
    ijOrb = iTri(iOrb,jOrb)
    jSym = SGS%ISM(jOrb)
    ijSym = Mul(iSym,jSym)
    mCSFs = CIS%nCSF(Mul(PsiSym,ijSym))

    Eij_Psi(1:mCSFs) => Eij_Psi_X
    Eji_Psi(1:mCSFs) => Eji_Psi_X
    Elk_Psi(1:mCSFs) => Elk_Psi_X

    Eij_Psi(:) = Zero
    Eji_Psi(:) = Zero
    ! Operate with E_ij on |Psi> and produce E_ij|Psi>
    call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,PsiSym,Psi,Eij_Psi)
    ! Compute Dij to be distributed below
    if (iSym == jSym) D_ij = dot_product(Psi,Eij_Psi)
    if (iOrb /= jOrb) call SG_Epq_Psi(SGS,CIS,EXS,jOrb,iOrb,CPQ,PsiSym,Psi,Eji_Psi)

    ! Note, in the case of a RASSCF the resulting sigma vector is incomplete. For the case of E_ij
    ! the complete expansion of E_ij|Psi> would require CSFs that are outside the set of CSFs defining the
    ! RAS expansion. However, if we cap this sigma vector with vectors that have zero coefficient for these
    ! external CSFs we experience no errors.

    do kOrb=1,iOrb
      kSym = SGS%ISM(kOrb)
      lOrb = kOrb
      lSym = SGS%ISM(lOrb)
      klSym = Mul(kSym,lSym)
      if (ijSym /= klSym) cycle

      ! Add the -d_il E_jk term, here in the form of a -d_kl E_ij term of Piklj
      ikOrb = iTri(iOrb,kOrb)
      ljOrb = iTri(lOrb,jOrb)
      ikljOrb = iTri(ikOrb,ljOrb)
      P(ikljOrb) = P(ikljOrb)-D_ij
      if ((iOrb /= kOrb) .and. (lOrb > jOrb)) PA(ikljOrb) = PA(ikljOrb)-D_ij
      if ((iOrb /= kOrb) .and. (lOrb < jOrb)) PA(ikljOrb) = PA(ikljOrb)+D_ij
    end do

    ! kOrb>=lOrb
    do kOrb=iOrb,SGS%nLev
      kSym = SGS%ISM(kOrb)
      lOrb_Min = 1
      if (kOrb == iOrb) lOrb_Min = jOrb
      do lOrb=lOrb_Min,kOrb
        klOrb = iTri(kOrb,lOrb)
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

        klijOrb = iTri(klOrb,ijOrb)
        P(klijOrb) = P(klijOrb)+P_klij
        if ((iOrb /= jOrb) .and. (kOrb /= lOrb)) PA(klijOrb) = PA(klijOrb)+P_klij

        if (iOrb /= jOrb) then
          P_klij = dot_product(Elk_Psi,Eji_Psi)
          P(klijOrb) = P(klijOrb)+P_klij
          if (kOrb /= lOrb) PA(klijOrb) = PA(klijOrb)-P_klij
        end if

      end do
    end do
    nullify(Eij_Psi,Eji_Psi,Elk_Psi)

  end do
end do

call mma_deallocate(Elk_Psi_X)
call mma_deallocate(Eij_Psi_X)
call mma_deallocate(Eji_Psi_X)

P = Half*P
PA = Half*PA

end subroutine sg_two_pdm
