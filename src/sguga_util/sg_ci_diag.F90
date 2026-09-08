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

subroutine sg_hdiag_ref(SGS,CIS,EXS,nCSFs,CSFSym, &
                        TUVX_Tri,nTUVX_Tri, &
                        TU_Tri,nTU_Tri, &
                        Diag)

  use sguga, only: SGStruct, CIStruct, EXStruct, sg_Epq_Psi, sg_EpqEqp_diag
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One
  use Definitions, only: wp, iwp, u6

  implicit none

  type(SGStruct), intent(in)    :: SGS
  type(CIStruct), intent(in)    :: CIS
  type(EXStruct), intent(inout) :: EXS

  integer(kind=iwp), intent(in) :: nCSFs
  integer(kind=iwp), intent(in) :: CSFSym

  integer(kind=iwp), intent(in) :: nTUVX_Tri
  integer(kind=iwp), intent(in) :: nTU_Tri

  real(kind=wp), intent(in) :: TUVX_Tri(nTUVX_Tri)
  real(kind=wp), intent(in) :: TU_Tri(nTU_Tri)

  real(kind=wp), intent(out) :: Diag(nCSFs)

  real(kind=wp), allocatable :: Psi(:)
  real(kind=wp), allocatable :: Sigma(:)

  integer(kind=iwp) :: iCSF

#define _TEST_
#ifdef _TEST_
  logical(kind=iwp) :: ReUse_Save
  integer(kind=iwp) :: p, q
  real(kind=wp), allocatable :: Epq(:), Eqp(:)
#endif

  call mma_allocate(Psi,nCSFs,Label='Psi')
  call mma_allocate(Sigma,nCSFs,Label='Sigma')

#ifdef _TEST_
! turn off reuse of SGTMP
ReUse_Save=EXS%Reuse_SGTMP
EXS%Reuse_SGTMP=.False.
call mma_allocate(Epq,nCSFs,Label='Epq')
call mma_allocate(Eqp,nCSFs,Label='Eqp')

do p=1,SGS%nLev
do q=1,SGS%nLev
Write (6,*) 'p,q=',p,q

do iCSF=1,nCSFs
  Psi(:) = Zero
  Psi(iCSF) = One
  Sigma(:) = Zero
  call SG_Epq_Psi(SGS,CIS,EXS, &
                  p,q, &
                  One, &
                  CSFSym, &
                  Psi, &
                  Sigma)
! Call RecPrt('Sigma',' ',Sigma,1,nCSFs)
  Epq(iCSF)=Dot_product(Sigma,Sigma)
end do

Call RecPrt('EqpEpq(exact)',' ',Epq,1,nCSFs)

Psi(:) = One
Sigma(:) = Zero
If (p==q) Then
  call SG_Epq_Psi(SGS,CIS,EXS, &
                  p,q, &
                  One, &
                  CSFSym, &
                  Psi, &
                  Sigma)
  Sigma(:) = Sigma(:)**2
Else
  call SG_EpqEqp_Diag(SGS,CIS,EXS, &
                  p,q, &
                  One, &
                  CSFSym, &
                  Psi, &
                  Sigma)
End If
Call RecPrt('EqpEpq(Direct)',' ',Sigma,1,nCSFs)

end do
end do

call mma_deallocate(Eqp)
call mma_deallocate(Epq)

  EXS%Reuse_SGTMP=ReUse_Save
Call Abend()
#endif

  do iCSF=1,nCSFs

     Psi(:)   = Zero
     Sigma(:) = Zero

     Psi(iCSF) = One

     call sg_h_psi(SGS,CIS,EXS, &
                   Psi,nCSFs,CSFSym, &
                   Sigma, &
                   TUVX_Tri,nTUVX_Tri, &
                   TU_Tri,nTU_Tri)

     Diag(iCSF) = Sigma(iCSF)

  end do
   call RecPrt('Diag',' ',Diag,1,nCSFs)

  call mma_deallocate(Sigma)
  call mma_deallocate(Psi)
  Call Abend()

end subroutine sg_hdiag_ref

subroutine sg_hdiag_fast(SGS,CIS,EXS,nCSFs,CSFSym, &
                         TUVX_Tri,nTUVX_Tri, &
                         TU_Tri,nTU_Tri, &
                         Diag)

  use sguga, only: SGStruct, CIStruct, EXStruct, SG_Epq_Psi
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Half
  use Definitions, only: wp, iwp
  use Index_Functions, only: iTri

  implicit none

  type(SGStruct), intent(in)    :: SGS
  type(CIStruct), intent(in)    :: CIS
  type(EXStruct), intent(inout) :: EXS

  integer(kind=iwp), intent(in) :: nCSFs
  integer(kind=iwp), intent(in) :: CSFSym

  integer(kind=iwp), intent(in) :: nTUVX_Tri
  integer(kind=iwp), intent(in) :: nTU_Tri

  real(kind=wp), intent(in)  :: TUVX_Tri(nTUVX_Tri)
  real(kind=wp), intent(in)  :: TU_Tri(nTU_Tri)

  real(kind=wp), intent(out) :: Diag(nCSFs)

  real(kind=wp), allocatable :: Psi(:)
  real(kind=wp), allocatable :: Epp(:)
  real(kind=wp), allocatable :: EqqEpp(:)
  real(kind=wp), allocatable :: EqpEpq(:)

  integer(kind=iwp) :: p,q
  integer(kind=iwp) :: ppTri
  integer(kind=iwp) :: pqTri
  integer(kind=iwp) :: ppqqTri
  integer(kind=iwp) :: pqqpTri

  real(kind=wp) :: OneInt
  real(kind=wp) :: Coul
  real(kind=wp) :: Exch

  Diag(:) = Zero

  call mma_allocate(Psi,nCSFs,Label='Psi')
  call mma_allocate(Epp,nCSFs,Label='Epp')
  call mma_allocate(EqqEpp,nCSFs,Label='EqqEpp')
  call mma_allocate(EqpEpq,nCSFs,Label='EqpEpq')

  !
  ! Trick:
  !
  ! Psi(:)=1 allows SG_Epq_Psi to return the diagonal
  ! expectation values for ALL CSFs simultaneously.
  !
  Psi(:) = One

  !
  ! One-electron contribution
  !
  do p=1,SGS%nLev

     Epp(:) = Zero

     call SG_Epq_Psi(SGS,CIS,EXS, &
                     p,p, &
                     One, &
                     CSFSym, &
                     Psi, &
                     Epp)

     ppTri = iTri(p,p)

     OneInt = TU_Tri(ppTri)

     Diag(:) = Diag(:) + OneInt*Epp(:)

  end do

  !
  ! Two-electron Coulomb / exchange contribution
  !
  do p=1,SGS%nLev

     Epp(:) = Zero

     call SG_Epq_Psi(SGS,CIS,EXS, &
                     p,p, &
                     One, &
                     CSFSym, &
                     Psi, &
                     Epp)

     do q=1,SGS%nLev

        !
        ! <Eqq Epp>
        !
        EqqEpp(:) = Zero

        call SG_Epq_Psi(SGS,CIS,EXS, &
                        q,q, &
                        One, &
                        CSFSym, &
                        Epp, &
                        EqqEpp)

        !
        ! <Eqp Epq>
        !
        EqpEpq(:) = Zero

        call SG_Epq_Psi(SGS,CIS,EXS, &
                        p,q, &
                        One, &
                        CSFSym, &
                        Psi, &
                        EqpEpq)

        call SG_Epq_Psi(SGS,CIS,EXS, &
                        q,p, &
                        One, &
                        CSFSym, &
                        EqpEpq, &
                        EqpEpq)

        pqTri   = iTri(max(p,q),min(p,q))
        ppqqTri = iTri(iTri(p,p),iTri(q,q))
        pqqpTri = iTri(iTri(max(p,q),min(p,q)), &
                       iTri(max(p,q),min(p,q)))

        Coul = Half*TUVX_Tri(ppqqTri)
        Exch = Half*TUVX_Tri(pqqpTri)

        Diag(:) = Diag(:) &
                + Coul*EqqEpp(:) &
                - Exch*EqpEpq(:)

     end do

  end do

  call mma_deallocate(EqpEpq)
  call mma_deallocate(EqqEpp)
  call mma_deallocate(Epp)
  call mma_deallocate(Psi)

end subroutine sg_hdiag_fast
