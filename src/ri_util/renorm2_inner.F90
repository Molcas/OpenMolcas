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
! Copyright (C) 2008, Roland Lindh                                     *
!***********************************************************************

subroutine ReNorm2_Inner(iCnttp)
!***********************************************************************
!                                                                      *
!    Objective: Orthonormalize parts of the auxiliary basis set.       *
!                                                                      *
! Called from: Mk_RICD_Shells                                          *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theor. Chemi., Lund Univ., Sweden.*
!                                                                      *
!             Modified to transform the auxiliary basis to a true      *
!             Cholesky basis set while on TACC 2008 conference in      *
!             Songjiang District, Shanghai, China, 23-27 Sept. 2008.   *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use RI_procedures, only: Drv2El_Atomic_NoSym
use SOAO_Info, only: iAOtSO, nSOInf
use Real_Spherical, only: Sphere
use Basis_Info, only: dbsc, iCnttp_Dummy, Shells
use Sizes_of_Seward, only: S
use RICD_Info, only: Thrshld_CD
use Integral_interfaces, only: int_wrout
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iCnttp
#include "itmax.fh"
integer(kind=iwp) :: i, iAng, iAO, iBas, iCase, iCmp, iDisk, ij, ijF, ijS, ijS_req, ijT, iSeed, iShll, iShll_, iSO, j, jBas, jCmp, &
                     jiS, Keep_Shell, Lu_A, Lu_Q, m, nBasisi, nCmp, nExpi, nSO, nTest, nTInt_c
real(kind=wp) :: Thr_CB, ThrAO
logical(kind=iwp) :: In_Core
real(kind=wp), allocatable :: ADiag(:), Not_Used(:), QVec(:,:), TInt_c(:), TInt_d(:), Tmp(:,:)
procedure(int_wrout) :: Integral_RI_2
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! Let us now Gram-Schmidt orthonormalize the auxiliary basis for
! better numerics and balance.

! Update kOffAO and lOffAO to include the auxiliary basis too.

call Setup_OffAO()

! Set up transformation matrix from Cartesian to real spherical harmonics.

call Sphere(S%iAngMx)

call Flip_Flop(.false.) ! Contracted mode.

Thr_CB = max(1.0e-14_wp,Thrshld_CD*1.0e-10_wp)
ThrAO = Zero

!do iCnttp = 1, nCnttp
! Skip the dummy shell
if (iCnttp == iCnttp_dummy) return
! skip non-auxiliary basis sets
if (.not. dbsc(iCnttp)%Aux) return
!                                                                      *
!***********************************************************************
!                                                                      *
! Define some parameters to facilitate the atomic calculation

S%nShlls = dbsc(iCnttp)%nVal
nTest = dbsc(iCnttp)%nVal-1

! Define AOtSO

iAO = 0
iSO = 0
nSO = 0
do iAng=0,nTest
  iShll_ = dbsc(iCnttp)%iVal+iAng
  nCmp = nTri_Elem1(iAng)
  if (Shells(iShll_)%Prjct) nCmp = 2*iAng+1
  iSO = 0
  if (Shells(iShll_)%nBasis_C*Shells(iShll_)%nExp == 0) cycle
  do iCmp=1,nCmp
    iAO = iAO+1
    if (iAO > nSOInf) then
      write(u6,*) 'renorm2_inner: iAO>nSOInf'
      write(u6,*) 'iAO=',iAO
      write(u6,*) 'nSOInf=',nSOInf
      call Abend()
    end if
    iAOtSO(iAO,0) = iSO+1
    iSO = iSO+Shells(iShll_)%nBasis
  end do
  nSO = nSO+iSO
end do

ijS_req = 0
Keep_Shell = iTabMx
do iAng=0,nTest
  iShll = dbsc(iCnttp)%iVal+iAng
  nExpi = Shells(iShll)%nExp
  nBasisi = Shells(iShll)%nBasis
  if (nExpi*nBasisi == 0) cycle

  nCmp = nTri_Elem1(iAng)
  if (Shells(iShll)%Prjct) nCmp = 2*iAng+1

  ijS_req = ijS_req+1

  call Drv2El_Atomic_NoSym(Integral_RI_2,ThrAO,iCnttp,iCnttp,TInt_c,nTInt_c,In_Core,Not_Used,Lu_A,ijS_req,Keep_Shell)
# ifdef _DEBUGPRINT_
  call TriPrt('TInt_c',' ',TInt_c,nTInt_c)
# endif

  if (.not. In_Core) then
    call WarningMessage(2,'Error in ReNorm')
    write(u6,*) 'Out-of-core acCD not implemented!'
    call Abend()
  end if

  ! Produce the reduced set, in-place reduction.

  call mma_allocate(TInt_d,nTInt_c**2,Label='TInt_d')
  ijT = 0
  do iBas=1,nTInt_c
    do jBas=1,iBas
      ijT = ijT+1
      ijS = (jBas-1)*nTInt_c+iBas
      jiS = (iBas-1)*nTInt_c+jBas
      TInt_d(ijS) = TInt_c(ijT)
      TInt_d(jiS) = TInt_c(ijT)
    end do
  end do
  call mma_deallocate(TInt_c)
# ifdef _DEBUGPRINT_
  call RecPrt('TInt_d',' ',TInt_d,nTInt_c,nTInt_c)
# endif

  ij = 0
  iCmp = 1
  jCmp = 1
  do jBas=1,nBasisi
    j = (jCmp-1)*nBasisi+jBas
    do iBas=1,nBasisi
      i = (iCmp-1)*nBasisi+iBas
      ijF = (j-1)*nBasisi*nCmp+i
      ij = ij+1

      TInt_d(ij) = TInt_d(ijF)

    end do
  end do
# ifdef _DEBUGPRINT_
  call RecPrt('TInt_d(r)','(5G20.10)',TInt_d,nBasisi,nBasisi)
# endif

  call mma_allocate(ADiag,nBasisi,Label=' ADiag')

  iSeed = 77
  Lu_A = IsFreeUnit(iSeed)
  call DaName_MF_WA(Lu_A,'AMat09')

  iDisk = 0
  call dDaFile(Lu_A,1,TInt_d,nBasisi**2,iDisk)

  iSeed = iSeed+1
  Lu_Q = IsFreeUnit(iSeed)
  call DaName_MF_WA(Lu_Q,'QMat09')

  call dcopy_(nBasisi,TInt_d,nBasisi+1,ADiag,1)

  call CD_AInv_Inner(nBasisi,m,ADiag,Lu_A,Lu_Q,Thr_CB)

  call mma_deallocate(ADiag)
  call mma_deallocate(TInt_d)

  ! Transform the contraction coefficients according to the Cholesky vectors.

  call mma_allocate(Tmp,nExpi,nBasisi,Label='Tmp')
  call mma_allocate(QVec,nBasisi,nBasisi,Label='QVec')
  QVec(:,:) = Zero

  iDisk = 0
  call dDaFile(Lu_Q,2,QVec,nBasisi*m,iDisk)
  call DaEras(Lu_Q)
# ifdef _DEBUGPRINT_
  call RecPrt('QVec',' ',QVec,nBasisi,m)
# endif

  do iCase=1,2
    Tmp(:,:) = Shells(iShll)%Cff_c(:,:,iCase)
#   ifdef _DEBUGPRINT_
    call RecPrt('Coeff(old)',' ',Shells(iShll)%Cff_c(:,:,iCase),nExpi,nBasisi)
#   endif
    call DGEMM_('N','N',nExpi,nBasisi,nBasisi,One,Tmp,nExpi,QVec,nBasisi,Zero,Shells(iShll)%Cff_c(:,:,iCase),nExpi)
#   ifdef _DEBUGPRINT_
    call RecPrt('Coeff(new)',' ',Shells(iShll)%Cff_c(:,:,iCase),nExpi,nBasisi)
#   endif
  end do

  call mma_deallocate(QVec)
  call mma_deallocate(Tmp)

end do

!end do
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine ReNorm2_Inner
