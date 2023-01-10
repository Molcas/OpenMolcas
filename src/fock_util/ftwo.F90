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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!               1998, Roland Lindh                                     *
!***********************************************************************

subroutine Ftwo(icase,ExFac,iSym,kSym,iBas,kBas,off_sqMat,off_ltMat,D1I,FI,D1A,FA,PQRS)
!***********************************************************************
!                                                                      *
!     Assemble Fock matrices FI and FA (in AO-basis)                   *
!                                                                      *
!     calling arguments:                                               *
!     icase   : input, integer                                         *
!               symmetry case number; (II!II)=1, (II!KK)=2, (IK!IK)=3  *
!     ExFac   : input, real*8                                          *
!               coefficient of "exact exchange"                        *
!     iSym    : input, integer                                         *
!               symmetry species iSym                                  *
!     kSym    : input, integer                                         *
!               symmetry species kSym                                  *
!     iBas    : input, integer                                         *
!               1-st, fixed index of ERIs                              *
!     kBas    : input, integer                                         *
!               2-nd, fixed index of ERIs                              *
!     off_sqMat : input, array of integer                              *
!               offset of one-electron integrals (squared format)      *
!     off_ltMat : input, array of integer                              *
!               offset of one-electron integrals (lower triangular)    *
!     D1I     : input, array of real*8                                 *
!               one-body density matrix (frozen+inactive, AO-basis)    *
!     FI      : output, array of real*8                                *
!               Fock matrix (frozen+inactive, AO-basis)                *
!     D1A     : input, array of real*8                                 *
!               one-body density matrix (active, AO-basis)             *
!     FA      : output, array of real*8                                *
!               Fock matrix (active, AO-basis)                         *
!     PQRS    : input, array of real*8                                 *
!               two-electron integrals (AO-basis)                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!     Modified to only need the unique symmetry blocks, R. Lindh,      *
!     March '98.                                                       *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: icase, iSym, kSym, iBas, kBas, off_sqMat(*), off_ltMat(*)
real(kind=wp), intent(in) :: Exfac, D1I(*), D1A(*), PQRS(*)
real(kind=wp), intent(inout) :: FI(*), FA(*)
#include "rasdim.fh"
#include "general.fh"
integer(kind=iwp) :: iOff, iOff_ij, iOff_kl, jBas, k, kl, kOff, l
real(kind=wp) :: D1A_ij, D1I_ij
real(kind=wp), external :: dDot_

!                                                                      *
!***********************************************************************
!                                                                      *
! symmetry case (II|II)
! Both Coulomb and Exchange contributions
!                                                                      *
!***********************************************************************
!                                                                      *
if (icase == 1) then

  ! Coulombic contribution

  iOff = off_ltMat(iSym)+nTri_Elem(iBas-1)+kBas
  kOff = off_sqMat(kSym)+1
  FI(iOff) = FI(iOff)+dDot_(nBas(kSym)*nBas(kSym),D1I(kOff),1,PQRS,1)
  FA(iOff) = FA(iOff)+dDot_(nBas(kSym)*nBas(kSym),D1A(kOff),1,PQRS,1)

  ! Exchange contribution

  if (ExFac /= Zero) then
    iOff = off_sqMat(iSym)+(iBas-1)*nBas(iSym)+1
    kOff = off_ltMat(kSym)+nTri_Elem(kBas-1)+1
    !call DGeMX(kBas,nBas(iSym),-Half*ExFac,PQRS,nBas(iSym),D1I(iOff),1,FI(kOff),1)
    !call DGeMX(kBas,nBas(iSym),-Half*ExFac,PQRS,nBas(iSym),D1A(iOff),1,FA(kOff),1)
    call DGEMV_('N',kBas,nBas(iSym),-Half*ExFac,PQRS,nBas(iSym),D1I(iOff),1,One,FI(kOff),1)
    call DGEMV_('N',kBas,nBas(iSym),-Half*ExFac,PQRS,nBas(iSym),D1A(iOff),1,One,FA(kOff),1)
    if (iBas /= kBas) then
      iOff = off_ltMat(iSym)+nTri_Elem(iBas-1)+1
      kOff = off_sqMat(kSym)+(kBas-1)*nBas(kSym)+1
      !call DGeMX(iBas,nBas(iSym),-Half*ExFac,PQRS,nBas(iSym),D1I(kOff),1,FI(iOff),1)
      !call DGeMX(iBas,nBas(iSym),-Half*ExFac,PQRS,nBas(iSym),D1A(kOff),1,FA(iOff),1)
      call DGEMV_('N',iBas,nBas(iSym),-Half*ExFac,PQRS,nBas(iSym),D1I(kOff),1,One,FI(iOff),1)
      call DGEMV_('N',iBas,nBas(iSym),-Half*ExFac,PQRS,nBas(iSym),D1A(kOff),1,One,FA(iOff),1)
    end if
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! symmetry case (II|KK) and (KK|II)
! Only Coulomb contributions
! This is limited to iBas >= jBas
! Code modified to only require uniqe symmetry blocks of
! two-electron integrals. R. Lindh, March '98.
!                                                                      *
!***********************************************************************
!                                                                      *
if ((icase == 2) .and. (iSym > kSym)) then
  jBas = kBas ! To avoid confusion.
  iOff_ij = off_ltMat(iSym)+nTri_Elem(iBas-1)+jBas
  kOff = off_sqMat(kSym)+1

  FI(iOff_ij) = FI(iOff_ij)+DDot_(nBas(kSym)**2,D1I(kOff),1,PQRS,1)
  FA(iOff_ij) = FA(iOff_ij)+DDot_(nBas(kSym)**2,D1A(kOff),1,PQRS,1)

  D1I_ij = D1I(off_sqMat(iSym)+(jBas-1)*nBas(iSym)+iBas)
  D1A_ij = D1A(off_sqMat(iSym)+(jBas-1)*nBas(iSym)+iBas)
  if (iBas /= jBas) then
    D1I_ij = D1I_ij*Two
    D1A_ij = D1A_ij*Two
  end if
  do k=1,nBas(kSym)
    do l=1,k
      kl = (l-1)*nBas(kSym)+k
      iOff_kl = off_ltMat(kSym)+nTri_Elem(k-1)+l
      FI(iOff_kl) = FI(iOff_kl)+D1I_ij*PQRS(kl)
      FA(iOff_kl) = FA(iOff_kl)+D1A_ij*PQRS(kl)
    end do
  end do

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! symmetry case (IK!IK)
! Only Exchange contributions
!                                                                      *
!***********************************************************************
!                                                                      *
if ((icase == 3) .and. (ExFac /= Zero)) then
  iOff = off_sqMat(iSym)+(iBas-1)*nBas(iSym)+1
  kOff = off_ltMat(kSym)+nTri_Elem(kBas-1)+1
  !call DGeMX(kBas,nBas(iSym),-Half*ExFac,PQRS,nBas(kSym),D1I(iOff),1,FI(kOff),1)
  !call DGeMX(kBas,nBas(iSym),-Half*ExFac,PQRS,nBas(kSym),D1A(iOff),1,FA(kOff),1)
  call DGEMV_('N',kBas,nBas(iSym),-Half*ExFac,PQRS,nBas(kSym),D1I(iOff),1,One,FI(kOff),1)
  call DGEMV_('N',kBas,nBas(iSym),-Half*ExFac,PQRS,nBas(kSym),D1A(iOff),1,One,FA(kOff),1)
  iOff = off_ltMat(iSym)+nTri_Elem(iBas-1)+1
  kOff = off_sqMat(kSym)+(kBas-1)*nBas(kSym)+1
  !call DGeMTX(nBas(kSym),iBas,-Half*ExFac,PQRS,nBas(kSym),D1I(kOff),1,FI(iOff),1)
  !call DGeMTX(nBas(kSym),iBas,-Half*ExFac,PQRS,nBas(kSym),D1A(kOff),1,FA(iOff),1)
  call DGEMV_('T',nBas(kSym),iBas,-Half*ExFac,PQRS,nBas(kSym),D1I(kOff),1,One,FI(iOff),1)
  call DGEMV_('T',nBas(kSym),iBas,-Half*ExFac,PQRS,nBas(kSym),D1A(kOff),1,One,FA(iOff),1)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Ftwo
