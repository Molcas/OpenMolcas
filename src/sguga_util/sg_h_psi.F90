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
! Copyright (C) 2025, Roland Lindh                                     *
!***********************************************************************

subroutine sg_h_psi(SGS,CIS,EXS,Psi,nCSFs,PsiSym,Sigma,TUVX_Tri,nTUVX_Tri,TU_Tri,nTU_Tri)

use Index_functions, only: iTri
use symmetry_info, only: MUL
use sguga, only: CIStruct, EXStruct, sg_epq_psi, SGStruct
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
type(SGStruct), intent(in) :: SGS
type(CIStruct), intent(in) :: CIS
type(EXStruct), intent(inout) :: EXS
integer(kind=iwp), intent(in) :: nCSFs, PsiSym, nTUVX_Tri, nTU_Tri
real(kind=wp), intent(in) :: Psi(nCSFs), TUVX_Tri(nTUVX_Tri), TU_Tri(nTU_Tri)
real(kind=wp), intent(out) :: Sigma(nCSFs)

integer(kind=iwp), parameter :: incx = 1, incy = 1, nBuff = 10
real(kind=wp), parameter :: Alpha = One, Beta = One
integer(kind=iwp) :: i, iBuff, ijOrb, ijSym, ikljOrb, ikOrb, iOrb, iSym, jOrb, jSym, klijOrb, klOrb, klSym, kOrb, kSym, ljOrb, &
                     lOrb, lOrb_Max, lSym, MaxDim, mCSFs, nOrb, SigmaSym
real(kind=wp) :: CPQ, Fact, OneInt, TUVX(nBuff), TwoInt
real(kind=wp), pointer :: Eij_Psi(:)
real(kind=wp), allocatable, target :: Ekl_Eij_Psi(:,:), Eij_Psi_X(:)

! SGUGA driven algorithm for H|Psi>, DGEMV version
! i>=j Symmetrize: E_ij + E_ji, k>=l Symmetrize E_kl + E_lj, ij>=kl

Sigma(:) = Zero

MaxDim = maxval(CIS%nCSF(:))
call mma_allocate(Eij_Psi_X,MaxDim,Label='Eij_Psi_X')
call mma_allocate(Ekl_Eij_Psi,nCSFs,nBuff,Label='Ekl_Eij_Psi')

OneInt = Zero  ! Cardholder variable
TwoInt = Zero  ! Cardholder variable

! The E_kl E_ij part

CPQ = One
nOrb = SGS%nLev

do iOrb=1,nOrb
  iSym = SGS%ISM(iOrb)
  do jOrb=1,iOrb
    ijOrb = iTri(iOrb,jOrb)

    jSym = SGS%ISM(jOrb)
    ijSym = MUL(iSym,jSym)
    SigmaSym = MUL(PsiSym,ijSym)
    mCSFs = CIS%nCSF(SigmaSym)

    Eij_Psi(1:mCSFs) => Eij_Psi_X

    ! Operate with E_ij on |Psi> and produce E_ij|Psi>
    Eij_Psi(:) = Zero
    call SG_Epq_Psi(SGS,CIS,EXS,iOrb,jOrb,CPQ,PsiSym,Psi,Eij_Psi)
    if (iOrb /= jOrb) call SG_Epq_Psi(SGS,CIS,EXS,jOrb,iOrb,CPQ,PsiSym,Psi,Eij_Psi)

    if (ijSym == 0) then
      ! ijSym==0, do the E_ij part of the sigma vector
      TwoInt = Zero
      do kOrb=1,nOrb
        lOrb = kOrb
        ikOrb = iTri(iOrb,kOrb)
        ljOrb = iTri(lOrb,jOrb)
        ikljOrb = iTri(ikOrb,ljOrb)
        TwoInt = TwoInt+TUVX_Tri(ikljOrb)
      end do
      OneInt = TU_Tri(ijOrb)-Half*TwoInt
      !call DaXpY_(nCSFs,OneInt,Eij_Psi(:),1,Sigma(:),1)
      Sigma(:) = Sigma(:)+OneInt*Eij_Psi(:)
    end if

    iBuff = 0
    do kOrb=1,iOrb
      kSym = SGS%ISM(kOrb)
      lOrb_Max = kOrb
      if (kOrb == iOrb) lOrb_Max = jOrb
      do lOrb=1,lOrb_Max
        klOrb = iTri(kOrb,lOrb)

        lSym = SGS%ISM(lOrb)
        klSym = MUL(kSym,lSym)

        if (ijSym /= klSym) cycle

        Fact = One
        if (klOrb /= ijOrb) Fact = Two

        ! Operate with E_kl on E_ij_|Psi> and produce E_kl_E_ij_|Psi>
        iBuff = iBuff+1
        Ekl_Eij_Psi(1:nCSFs,iBuff) = Zero

        call SG_Epq_Psi(SGS,CIS,EXS,kOrb,lOrb,CPQ,SigmaSym,Eij_Psi,Ekl_Eij_Psi)
        if (kOrb /= lOrb) call SG_Epq_Psi(SGS,CIS,EXS,lOrb,kOrb,CPQ,SigmaSym,Eij_Psi,Ekl_Eij_Psi)

        klijOrb = iTri(klOrb,ijOrb)
        TUVX(iBuff) = Fact*Half*TUVX_Tri(klijOrb)

        if (iBuff == nBuff) then
          !call DGEMV_('N',nCSFs,nBuff,Alpha,Ekl_Eij_Psi(1:nCSFs,1:nBuff),nCSFs,TUVX(1:nBuff),incx,Beta,Sigma,incy)
          !Sigma(:) = Sigma(:)+Mat_Mul(Ekl_Eij_Psi(1:nCSFs,1:nBuff),TUVX(1:nBuff))
          do i=1,nBuff
            Sigma(:) = Sigma(:)+Ekl_Eij_Psi(:,i)*TUVX(i)
          end do
          iBuff = 0
          Ekl_Eij_Psi(:,:) = Zero
        end if

      end do
    end do

    if (iBuff /= 0) then
      !call DGEMV_('N',nCSFs,iBuff,Alpha,Ekl_Eij_Psi(1:nCSFs,1:iBuff),nCSFs,TUVX(1:iBuff),incx,Beta,Sigma,incy)
      !Sigma(:) = Sigma(:)+Mat_Mul(Ekl_Eij_Psi(1:nCSFs,1:nBuff),TUVX(1:nBuff))
      do i=1,iBuff
        Sigma(:) = Sigma(:)+Ekl_Eij_Psi(:,i)*TUVX(i)
      end do
      iBuff = 0
    end if

    nullify(Eij_Psi)

  end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call mma_deallocate(Ekl_Eij_Psi)
call mma_deallocate(Eij_Psi_X)

end subroutine sg_h_psi
