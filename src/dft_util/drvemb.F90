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
! Copyright (C) 2010,2012,2017, Francesco Aquilante                    *
!               2015,2017, Alexander Zech                              *
!***********************************************************************

subroutine DrvEMB(nh1,KSDFT,Do_Grad,Grad,nGrad,DFTFOCK)
!***********************************************************************
!***********************************************************************
!** Orbital-Free Embedding calculation                               ***
!**                                                                  ***
!** Method:                                                          ***
!**     T. A. Wesolowski, A. Warshel, J. Phys. Chem. 97 (1993) 8050. ***
!**                                                                  ***
!** NDSD potential:                                                  ***
!**     J.-M. Garcia Lastra, J. W. Kaminski, T. A. Wesolowski,       ***
!**                               J. Chem. Phys.  129 (2008) 074107. ***
!**                                                                  ***
!** Embedding multi-determinantal wfs:                               ***
!**     T. A. Wesolowski, Phys. Rev.A. 77 (2008) 012504.             ***
!**                                                                  ***
!**                                                                  ***
!** Embedding Hartree-Fock wf:                                       ***
!**     F. Aquilante, T. A. Wesolowski                               ***
!**                       J. Chem. Phys. 135 (2011) 084120.          ***
!**                                                                  ***
!**                                                                  ***
!** Author: F. Aquilante, Geneva July 2010                           ***
!**                                                                  ***
!**                       (last update: Feb 2012)                    ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************

use OFembed, only: dFMD, Energy_NAD, Func_A, Func_AB, Func_B, NDSD, OFE_first, V_emb, V_Nuc_AB, V_Nuc_BA, Xsigma
use nq_Info, only: Dens_I
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nh1, nGrad
character(len=*), intent(inout) :: KSDFT
logical(kind=iwp), intent(in) :: Do_Grad
real(kind=wp), intent(inout) :: Grad(nGrad)
character(len=4), intent(in) :: DFTFOCK
integer(kind=iwp) :: i, iSpin, j, kSpin, nD, nFckDim
real(kind=wp) :: d_Alpha, d_Beta, DSpn, DTot, Ec_A, Fact, Fact_, Fakt_, Func_A_TF, Func_B_TF, tmp, xElAB, Vxc_ref(2)
logical(kind=iwp) :: is_rhoA_on_file
real(kind=wp), allocatable :: D_DS(:,:), F_DFT(:,:), Fcorr(:,:), TmpA(:)
real(kind=wp), external :: dDot_, Xlambda
#ifdef _NOT_USED_
integer(kind=iwp) :: nDens
real(kind=wp) :: Func_AB_TF, TF_NAD, V_emb_x, V_emb_x_ref, Xint_Ts_A, Xint_Ts_AB, Xint_Ts_NAD, Xnorm, Ynorm
real(kind=wp), allocatable :: D1ao_x(:), D1ao_y(:), Vemb(:)
#endif

is_rhoA_on_file = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
call Setup_iSD()
if (Do_Grad) Grad(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _NOT_USED_
! Section to calculate Nonelectr. V_emb with current density
! Temporarily turned off (clean output)
if (.not. OFE_first) then
  call mma_allocate(D1ao_y,nh1)
  call NameRun('AUXRFIL') ! switch RUNFILE name
  call mma_allocate(Vemb,nh1,label='Vemb')
  call Get_dArray('dExcdRa',Vemb,nh1)
  call mma_allocate(TmpA,nh1,Label='TmpA')
  call Get_dArray('Nuc Potential',TmpA,nh1)
  ! Subtract V_nuc_B
  Vemb(:) = Vemb-TmpA
  ! Calculate nonelectr. V_emb with current Density
  Ynorm = dDot_(nh1,D1ao_y,1,D1ao_y,1)
  V_emb_x = dDot_(nh1,Vemb,1,D1ao_y,1)
  write(u6,'(A,F19.10,4X,A,F10.5)') 'Nonelectr. Vemb w. current density: ',V_emb_x,'Y_Norm = ',Ynorm
  call mma_deallocate(D1ao_y)
  ! Get rho_A_ref
  call NameRun('PRERFIL')
  call mma_allocate(D1ao_x,nDens,Label='D1ao_x')
  call get_dArray('D1ao',D1ao_x,nDens)
  Xnorm = dDot_(nh1,D1ao_x,1,D1ao_x,1)
  V_emb_x_ref = dDot_(nh1,Vemb,1,D1ao_x,1)
  write(u6,'(A,F19.10,4X,A,F10.5)') 'Nonelectr. Vemb w.    ref. density: ',V_emb_x_ref,'X_Norm = ',Xnorm
  call VEMB_Exc_states(Vemb,nh1,KSDFT,Func_B)
  call mma_deallocate(TmpA)
  call mma_deallocate(D1ao_x)
  call mma_dealloacte(Vemb)
  call NameRun('#Pop')   ! switch back to AUXRFIL
  call NameRun('#Pop')   ! switch back to RUNFILE
end if
! Section End
#endif
call f_Inquire('PRERFIL',is_rhoA_on_file) ! rho_A from file
if (is_rhoA_on_file .and. (.not. OFE_first)) return ! Vemb on disk

!***********************************************************************
!                                                                      *
!     Setup of density matrices for subsys B (environment)             *
!                                                                      *
!***********************************************************************
call NameRun('AUXRFIL') ! switch RUNFILE name
!                                                                      *
!***********************************************************************
!                                                                      *
nD = 4
call mma_allocate(F_DFT,nh1,nD,Label='F_DFT')
call mma_allocate(D_DS,nh1,nD,Label='D_DS')
Vxc_ref(1) = Zero
Vxc_ref(2) = Zero

! Get the density matrix of the environment (rho_B)

call Get_iScalar('Multiplicity',kSpin)
call Get_dArray_chk('D1ao',D_DS(:,1),nh1)
!call RecPrt('D1ao',' ',D_DS(:,1),nh1,1)

! Get the spin density matrix of the environment

if (kSpin /= 1) then
  call Get_dArray_chk('D1sao',D_DS(:,2),nh1)
  !call RecPrt('D1Sao',' ',D_DS(:,2),nh1,1)
end if

! Compute alpha and beta density matrices of the environment

nFckDim = 2
if (kSpin == 1) then
  D_DS(:,1) = Half*D_DS(:,1)
  D_DS(:,2) = D_DS(:,1)
  nFckDim = 1
else
  do i=1,nh1
    DTot = D_DS(i,1)
    DSpn = D_DS(i,2)
    d_Alpha = Half*(DTot+DSpn)
    d_Beta = Half*(DTot-DSpn)
    D_DS(i,1) = d_Alpha
    D_DS(i,2) = d_Beta
  end do
  !call RecPrt('Da',' ',D_DS(:,1),nh1,1)
  !call RecPrt('Db',' ',D_DS(:,2),nh1,1)
end if

!if (OFE_first) then
!---AZECH 10/2015
! kinetic part of E_xct, Subsys B
Func_B_TF = Zero
call wrap_DrvNQ('TF_only',F_DFT(:,1:nFckDim),nFckDim,Func_B_TF,D_DS(:,1:nFckDim),nh1,nFckDim,Do_Grad,Grad,nGrad,DFTFOCK)

if (OFE_first) then

  call wrap_DrvNQ(KSDFT,F_DFT(:,1:nFckDim),nFckDim,Func_B,D_DS(:,1:nFckDim),nh1,nFckDim,Do_Grad,Grad,nGrad,DFTFOCK)

  if (KSDFT(1:4) == 'NDSD') then
    call mma_allocate(NDSD,nh1,nFckDim,label='NDSD')
    NDSD(:,:) = F_DFT(:,1:nFckDim)
    KSDFT(1:4) = 'LDTF' !set to Thomas-Fermi for subsequent calls
  end if

end if

!***********************************************************************
!                                                                      *
!     Setup of density matrices for subsys A                           *
!                                                                      *
!***********************************************************************
call NameRun('#Pop')     ! switch back RUNFILE name

if (is_rhoA_on_file) call NameRun('PRERFIL')
! Get the density matrix for rho_A

call Get_dArray_chk('D1ao',D_DS(:,3),nh1)
!call RecPrt('D1ao',' ',D_DS(:,3),nh1,1)

call Get_iScalar('Multiplicity',iSpin)
if ((iSpin == 1) .and. (kSpin /= 1) .and. OFE_first) then
  call WarningMessage(0,'Non-singlet environment perturbation on singlet state!;'// &
                      'Spin-components of the OFE potential will be averaged.')
end if

! Get the spin density matrix of A

if (iSpin /= 1) then
  call Get_dArray_chk('D1sao',D_DS(:,4),nh1)
  !call RecPrt('D1Sao',' ',D_DS(:,4),nh1,1)
end if

! Compute alpha and beta density matrices of subsystem A

nFckDim = 2
if (iSpin == 1) then
  D_DS(:,3) = Half*D_DS(:,3)
  D_DS(:,4) = D_DS(:,3)
  if (kSpin == 1) nFckDim = 1
else
  do i=1,nh1
    DTot = D_DS(i,3)
    DSpn = D_DS(i,4)
    d_Alpha = Half*(DTot+DSpn)
    d_Beta = Half*(DTot-DSpn)
    D_DS(i,3) = d_Alpha
    D_DS(i,4) = d_Beta
  end do
  !call RecPrt('Da',' ',D_DS(:,3),nh1,1)
  !call RecPrt('Db',' ',D_DS(:,4),nh1,1)
end if

!---AZECH 10/2015
! kinetic part of E_xct, Subsys A
call wrap_DrvNQ('TF_only',F_DFT(:,3:nFckDim+2),nFckDim,Func_A_TF,D_DS(:,3:nFckDim+2),nh1,nFckDim,Do_Grad,Grad,nGrad,DFTFOCK)

call wrap_DrvNQ(KSDFT,F_DFT(:,3:nFckDim+2),nFckDim,Func_A,D_DS(:,3:nFckDim+2),nh1,nFckDim,Do_Grad,Grad,nGrad,DFTFOCK)

! Fraction of correlation potential from A (cases: HF or Trunc. CI)

if (dFMD > Zero) then

  call mma_allocate(Fcorr,nh1,nFckDim,Label='Fcorr')

  call cwrap_DrvNQ(KSDFT,nFckDim,Ec_A,D_DS(:,3:nFckDim+2),nh1,nFckDim,Do_Grad,Grad,nGrad,DFTFOCK,Fcorr(:,1:nFckDim))
end if

!***********************************************************************
!                                                                      *
!     Calculation on the supermolecule                                 *
!                                                                      *
!***********************************************************************
nFckDim = 2
if ((iSpin == 1) .and. (kSpin == 1)) then
  nFckDim = 1
  D_DS(:,1) = D_DS(:,1)+D_DS(:,3)
else
  D_DS(:,1) = D_DS(:,1)+D_DS(:,3)
  D_DS(:,2) = D_DS(:,2)+D_DS(:,4)
end if
#ifdef _NOT_USED_
!---AZECH 10/2015
! kinetic part of E_xct, Subsys A+B
! temporarily turned off to clean output
if (.false.) then
  Func_AB_TF = Zero
  call wrap_DrvNQ('TF_only',F_DFT(:,1:nFckDim),nFckDim,Func_AB_TF,D_DS(:,1:nFckDim),nh1,nFckDim,Do_Grad,Grad,nGrad,DFTFOCK)
  TF_NAD = Func_AB_TF-Func_A_TF-Func_B_TF
  write(u6,*) 'kinetic part of E_xc,T (Thomas-Fermi ONLY)'
  write(u6,'(A,F19.10)') 'Ts(A+B): ',Func_AB_TF
  write(u6,'(A,F19.10)') 'Ts(A):   ',Func_A_TF
  write(u6,'(A,F19.10)') 'Ts(B):   ',Func_B_TF
  write(u6,'(A,F19.10)') '-------------------'
  write(u6,'(A,F19.10)') 'Ts_NAD:  ',TF_NAD
  ! calculate v_T, Subsys A+B
  Xint_Ts_AB = dDot_(nh1,F_DFT(:,1),1,D_DS(:,3),1)
  Xint_Ts_NAD = Xint_Ts_AB-Xint_Ts_A
  ! scale by 2 because wrapper only handles spin-densities
  Xint_Ts_NAD = Two*Xint_Ts_NAD
  write(u6,*) 'integrated v_Ts_NAD (Thomas-Fermi) with rhoA current'
  write(u6,'(A,F19.10)') 'Ts(A+B)_integral: ',Xint_Ts_AB
  write(u6,'(A,F19.10)') 'Ts(A)_integral:   ',Xint_Ts_A
  write(u6,'(A,F19.10)') '-------------------'
  write(u6,'(A,F19.10)') 'Ts_NAD_integral:  ',Xint_Ts_NAD
end if
#endif

call wrap_DrvNQ(KSDFT,F_DFT(:,1:nFckDim),nFckDim,Func_AB,D_DS(:,1:nFckDim),nh1,nFckDim,Do_Grad,Grad,nGrad,DFTFOCK)

Energy_NAD = Func_AB-Func_A-Func_B

!---AZECH 10/2015
! exchange-correlation part of E_xct, Subsys A+B
! temporarily turned off to clean output
!write(u6,*) 'E_xc_NAD (determined with Thomas-Fermi)'
!Func_xc_NAD = Energy_NAD-TF_NAD
!write(u6,'(A,F19.10)') 'E_xc_NAD: ',Func_xc_NAD

if (dFMD > Zero) then
  xElAB = Dens_I
  Fakt_ = -Xlambda(abs(Energy_NAD)/xElAB,Xsigma)
  F_DFT(:,3:nFckDim+2) = F_DFT(:,3:nFckDim+2)+Fakt_*Fcorr(:,:)
  call mma_deallocate(Fcorr)
# ifdef _DEBUGPRINT_
  write(u6,*) ' lambda(E_nad) = ',dFMD*Fakt_
# endif
end if

!                                                                      *
!***********************************************************************
!                                                                      *
! Non Additive (NAD) potential: F(AB)-F(A)
do i=1,nFckDim
  F_DFT(:,i) = F_DFT(:,i)-F_DFT(:,2+i)
end do

! NDSD potential for T_nad: add the (B)-dependent term
if (allocated(NDSD)) then
  j = 1
  do i=1,nFckDim
    F_DFT(:,i) = F_DFT(:,i)+NDSD(:,j)
    if (kSpin /= 1) j = j+1
  end do
end if

! Add the Nuc Attr potential (from subsystem B) and then
! put out the DFT Fock matrices from the (NAD) embedding potential
! on the runfile (AUXRFIL). Note that the classical Coulomb
! interaction potential from subsystem B is computed in the std
! Fock matrix builders

if (is_rhoA_on_file) call NameRun('#Pop')
call NameRun('AUXRFIL')  ! switch RUNFILE name

call mma_allocate(TmpA,nh1,Label='TmpA')
call Get_dArray('Nuc Potential',TmpA,nh1)

Fact = Two ! because Dmat has been scaled by half
if (kSpin /= 1) Fact = One
Fact_ = Fact

V_emb = Fact*dDot_(nh1,F_DFT(:,1),1,D_DS(:,3),1)
V_Nuc_AB = Fact*dDot_(nh1,TmpA,1,D_DS(:,3),1)
if (kSpin /= 1) then
  V_emb = V_emb+Fact*dDot_(nh1,F_DFT(:,2),1,D_DS(:,4),1)
  V_Nuc_AB = V_Nuc_AB+Fact*dDot_(nh1,TmpA,1,D_DS(:,4),1)
end if

! Averaging the spin-components of F(AB) iff non-spol(A)//spol(B)
if ((iSpin == 1) .and. (kSpin /= 1)) then
  do i=1,nh1
    tmp = Half*(F_DFT(i,1)+F_DFT(i,2))
    F_DFT(i,1) = tmp
  end do
  nFckDim = 1  ! reset stuff as if A+B had been spin compensated
  Fact = Two
end if

do i=1,nFckDim
  F_DFT(:,i) = F_DFT(:,i)+TmpA
  Vxc_ref(i) = Fact*dDot_(nh1,F_DFT(:,i),1,D_DS(:,i+2),1)
end do

if (dFMD > Zero) call Put_dScalar('KSDFT energy',Ec_A)
call Put_dArray('Vxc_ref ',Vxc_ref,2)

call Put_dArray('dExcdRa',F_DFT(:,1:nFckDim),nh1*nFckDim)
call NameRun('#Pop')  ! switch back RUNFILE name

call Get_dArray('Nuc Potential',TmpA,nh1)
V_Nuc_BA = Fact_*(dDot_(nh1,TmpA,1,D_DS(:,1),1)-dDot_(nh1,TmpA,1,D_DS(:,3),1))
if (kSpin /= 1) then
  V_Nuc_BA = V_Nuc_BA+Fact_*(dDot_(nh1,TmpA,1,D_DS(:,2),1)-dDot_(nh1,TmpA,1,D_DS(:,4),1))
end if

call mma_deallocate(TmpA)

#ifdef _DEBUGPRINT_
if (nFckDim == 1) then
  do i=1,nh1
    write(u6,'(i4,f22.16)') i,F_DFT(i,1)
  end do
else
  do i=1,nh1
    write(u6,'(i4,3f22.16)') i,F_DFT(i,1),F_DFT(i,2),Half*(F_DFT(i,1)+F_DFT(i,2))
  end do
end if
write(u6,'(a,f22.16)') ' NAD DFT Energy :',Energy_NAD
#endif

call mma_deallocate(F_DFT)
call mma_deallocate(D_DS)
call Free_iSD()

return

end subroutine DrvEMB
