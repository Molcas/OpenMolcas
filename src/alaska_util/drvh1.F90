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
! Copyright (C) 1991,1995, Roland Lindh                                *
!               Markus P. Fuelscher                                    *
!***********************************************************************

subroutine Drvh1(Grad,Temp,nGrad)
!***********************************************************************
!                                                                      *
! Object: driver for computation of gradient with respect to the one-  *
!         electron hamiltonian and the overlap matrix. The former will *
!         be contracted with the "variational" first order density     *
!         matrix and the latter will be contracted with the generalized*
!         Fock matrix.                                                 *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             October '91                                              *
!                                                                      *
! modified by M.P. Fuelscher                                           *
! make use of the RELAX file                                           *
!                                                                      *
!             Modified to Sperical Well, External Field, and SRO       *
!             gradients, April '95. R. Lindh                           *
!             Modified to Self Consistent Reaction Fields, May '95     *
!***********************************************************************

use PCM_arrays, only: PCM_SQ
use External_Centers, only: nWel, XF, Wel_Info
use Basis_Info, only: nCnttp, dbsc, nBas
use Symmetry_Info, only: nIrrep
#ifdef _NEXTFFIELD_
use finfld, only: force
#endif
use Index_Functions, only: nTri_Elem1
use Grd_interface, only: grd_kernel, grd_mem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(out) :: Temp(nGrad)
integer(kind=iwp) :: i, iComp, iCOSMO, ii, iIrrep, iMltpl, iPrint, iRout, iWel, ix, iy, nComp, nCompf, nDens, nFock, nOrdOp, nOrdOpf
real(kind=wp) :: Fact, TCpu1, TCpu2, TWall1, TWall2
character(len=80) :: Label
character(len=8) :: Method
logical(kind=iwp) :: DiffOp, lECP, lFAIEMP, lPP
integer(kind=iwp), allocatable :: lOper(:), lOperf(:)
real(kind=wp), allocatable :: Coor(:,:), Coorf(:,:), D_Var(:), Fock(:)
procedure(grd_kernel) :: COSGrd, FragPGrd, KneGrd, M1Grd, M2Grd, NAGrd, OvrGrd, PCMGrd, PPGrd, PrjGrd, RFGrd, SROGrd, WelGrd, XFdGrd
procedure(grd_mem) :: FragPMmG, KneMmG, M1MmG, M2MmG, NAMmG, OvrMmG, PCMMmG, PPMmG, PrjMmG, RFMmg, SROMmG, WelMmg, XFdMmg
#ifdef _NEXTFFIELD_
!AOM<
integer(kind=iwp) :: ncmp, nextfld
character(len=30) :: fldname
procedure(grd_kernel) :: MltGrd
procedure(grd_mem) :: MltMmG
!AOM>
#endif
#include "Molcas.fh"
#include "print.fh"
#include "disp.fh"
#include "wldata.fh"
#include "rctfld.fh"

! Prologue
iRout = 131
iPrint = nPrint(iRout)
call CWTime(TCpu1,TWall1)
call StatusLine(' Alaska:',' Computing 1-electron gradients')

! Allocate memory for density and Fock matrices

nFock = 0
nDens = 0
do iIrrep=0,nIrrep-1
  nFock = nFock+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  nDens = nDens+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do
lECP = .false.
lPP = .false.
lFAIEMP = .false.
do i=1,nCnttp
  lECP = lECP .or. dbsc(i)%ECP
  lPP = lPP .or. (dbsc(i)%nPP /= 0)
  lFAIEMP = LFAIEMP .or. dbsc(i)%Frag
end do

! Get the method label
!write(u6,*) ' Read Method label'
call Get_cArray('Relax Method',Method,8)

! Read the variational 1st order density matrix
! density matrix in AO/SO basis
!write(u6,*) ' Read density matrix'

call mma_allocate(D_Var,nDens,Label='D_Var')
call Get_D1ao_Var(D_Var,nDens)
if (iPrint >= 99) then
  write(u6,*) 'variational 1st order density matrix'
  ii = 1
  do iIrrep=0,nIrrep-1
    write(Label,*) 'symmetry block',iIrrep
    call TriPrt(Label,' ',D_Var(ii),nBas(iIrrep))
    ii = ii+nBas(iIrrep)*(nBas(iIrrep)+1)/2
  end do
end if

! Read the generalized Fock matrix
! Fock matrix in AO/SO basis
!write(u6,*) ' Read Fock matrix'
if (.not. HF_Force) then
  call mma_allocate(Fock,nDens,Label='Fock')
  call Get_dArray_chk('FockOcc',Fock,nDens)
  if (iPrint >= 99) then
    write(u6,*) 'generalized Fock matrix'
    ii = 1
    do iIrrep=0,nIrrep-1
      write(Label,*) 'symmetry block',iIrrep
      call TriPrt(Label,' ',Fock(ii),nBas(iIrrep))
      ii = ii+nBas(iIrrep)*(nBas(iIrrep)+1)/2
    end do
  end if
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! nOrdOp: order/rank of the operator
! lOper: lOper of each component of the operator

nOrdOp = 0
nComp = nTri_Elem1(nOrdOp)
call mma_allocate(Coor,3,nComp,Label='Coor')
call mma_allocate(lOper,nComp,Label='lOper')
Coor(:,:) = Zero
lOper(1) = 1
!***********************************************************************
!1)                                                                    *
!     Trace the generalized Fock matrix with the gradient of the       *
!     overlap matrix.                                                  *
!                                                                      *
!***********************************************************************

if (.not. HF_Force) then
  DiffOp = .false.
  Label = ' The Renormalization Contribution'
  call OneEl_g(OvrGrd,OvrMmG,Temp,nGrad,DiffOp,Coor,Fock,nFock,lOper,nComp,nOrdOp,Label)
  Grad(:) = Grad(:)-Temp(:)
end if

!***********************************************************************
!2)                                                                    *
!     Trace the "variational" first order density matrix with the      *
!     gradient of the kinetic energy integrals.                        *
!                                                                      *
!***********************************************************************

if (.not. HF_Force) then
  DiffOp = .false.
  Label = ' The Kinetic Energy Contribution'
  call OneEl_g(KneGrd,KneMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
  Grad(:) = Grad(:)+Temp(:)
  ! This is either not implemented or disabled, so commenting it out
# ifdef _NEXTFFIELD_
  !AOM<
  ! Check finite field operators...
  nextfld = 0
  do
    call GetNextFfield(nextfld,fldname,nOrdOpf,ncmp,force,size(force))
    if (nextfld == 0) exit
    if (nOrdOpf == 0) then
      if (fldname(1:2) /= 'OV') then
        call WarningMessage(2,'Error in Drvh1')
        write(u6,*) 'Finite field gradients only for MLTP'
        call Quit_OnUserError()
      end if
    end if
    nCompf = nTri_Elem1(nOrdOpf)
    Label = fldname
    if (ncompf /= ncmp) then
      call WarningMessage(2,'Error in Drvh1')
      write(u6,*) 'Wrong number of components in FF grad'
      call Quit_OnUserError()
    end if
    call mma_allocate(lOperf,nCompf,Label='lOperf')
    lOperf(:) = 1
    DiffOp = .false.
    if (nOrdOpf > 0) DiffOp = .true.
    call OneEl_g(MltGrd,MltMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOperf,nCompf,nOrdOpf,Label)
    call MltGrdNuc(Temp,nGrad,nOrdOpf)
    call mma_deallocate(lOperf)
    Grad(:) = Grad(:)-Temp(:)
  end do
  !AOM>
# endif
end if

!***********************************************************************
!3)                                                                    *
!     Trace the "variational" first order density matrix with the      *
!     gradient of the nuclear attraction integrals.                    *
!                                                                      *
!***********************************************************************

DiffOp = .true.
Label = ' The Nuclear Attraction Contribution'
call OneEl_g(NAGrd,NAMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
Grad(:) = Grad(:)+Temp(:)
if (HF_Force) then
  if (lECP .or. (nWel /= 0) .or. allocated(XF) .or. lRF) then
    call WarningMessage(2,'Error in Drvh1')
    write(u6,*) 'HF forces not implemented yet for this case!'
    call Quit_OnUserError()
  end if
end if

!***********************************************************************
!4)                                                                    *
!     Trace the "variational" first order density matrix with the      *
!     gradient of the ECP integrals.                                   *
!                                                                      *
!***********************************************************************

if (.not. HF_Force) then
  if (lECP) then
    DiffOp = .true.
    Label = ' The Projection Operator Contribution'
    call OneEl_g(PrjGrd,PrjMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
    Grad(:) = Grad(:)+Temp(:)

    Label = ' The M1 Operator Contribution'
    call OneEl_g(M1Grd,M1MmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
    Grad(:) = Grad(:)+Temp(:)

    Label = ' The M2 Operator Contribution'
    call OneEl_g(M2Grd,M2MmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
    Grad(:) = Grad(:)+Temp(:)

    Label = ' The SR Operator Contribution'
    call OneEl_g(SROGrd,SROMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
    Grad(:) = Grad(:)+Temp(:)
  end if
  if (lPP) then
    Label = ' The Pseudo Potential Contribution'
    call OneEl_g(PPGrd,PPMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
    Grad(:) = Grad(:)+Temp(:)
  end if
end if

!***********************************************************************
!5)                                                                    *
!     Trace the "variational" first order density matrix with the      *
!     gradient of the Spherical well integrals.                        *
!                                                                      *
!***********************************************************************

if (.not. HF_Force) then
  DiffOp = .true.
  do iWel=1,nWel
    r0 = Wel_Info(1,iWel)
    ExpB = Wel_Info(2,iWel)
    Label = ' The Spherical Well Contribution'
    call OneEl_g(WelGrd,WelMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
    Fact = Wel_Info(3,iWel)
    Grad(:) = Grad(:)+Fact*Temp(:)
  end do
end if

!***********************************************************************
!6)                                                                    *
!     Trace the "variational" first order density matrix with the      *
!     gradient of the external field integrals.                        *
!                                                                      *
!***********************************************************************

if (.not. HF_Force) then
  if (allocated(XF)) then
    DiffOp = .true.
    Label = ' The External Field Contribution'
    call OneEl_g(XFdGrd,XFdMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
    Grad(:) = Grad(:)+Temp(:)
  end if
end if

!***********************************************************************
!7)                                                                    *
!     Trace the "variational" first order density matrix with the      *
!     gradient of the reaction field integrals.                        *
!                                                                      *
!***********************************************************************

if (.not. HF_Force) then
  if (lRF .and. (.not. lLangevin) .and. (.not. PCM)) then

    ! The Kirkwood model

    nOrdOpf = lMax
    nCompf = (lMax+1)*(lMax+2)*(lMax+3)/6
    call mma_allocate(lOperf,nCompf,Label='lOperf')

    ! Store permutation symmetry of components of the EF

    iComp = 1
    do iMltpl=0,lMax
      do ix=iMltpl,0,-1
        !if (Mod(ix,2) == 0) then
        !  iSymX = 1
        !else
        !  ixyz = 1
        !  iSymX = 2**IrrFnc(ixyz)
        !end if
        do iy=iMltpl-ix,0,-1
          !if (Mod(iy,2) == 0) then
          !  iSymY = 1
          !else
          !  ixyz = 2
          !  iSymY = 2**IrrFnc(ixyz)
          !end if
          !iz = iMltpl-ix-iy
          !if (Mod(iz,2) == 0) then
          !  iSymZ = 1
          !else
          !  ixyz = 4
          !  iSymZ = 2**IrrFnc(ixyz)
          !end if
          !lOperf(iComp) = MltLbl(iSymX,MltLbl(iSymY,iSymZ))
          ! Compute only total symmetric contributions
          lOperf(iComp) = 1
          iComp = iComp+1
        end do
      end do
    end do
    call mma_allocate(Coorf,3,nCompf,Label='Coorf')
    Coorf(:,:) = Zero
    DiffOp = .true.
    Label = ' The Electronic Reaction Field Contribution'
    call OneEl_g(RFGrd,RFMmG,Temp,nGrad,DiffOp,Coorf,D_Var,nDens,lOperf,nCompf,nOrdOpf,Label)
    Grad(:) = Grad(:)+Temp(:)

    call mma_deallocate(lOperf)
    call mma_deallocate(Coorf)

  else if (lRF .and. PCM) then
    iCOSMO = 0
    ! The PCM / COSMO model

    if (iCOSMO <= 0) then
      iPrint = 15
      PCM_SQ(:,:) = PCM_SQ(:,:)/real(nIrrep,kind=wp)
    end if
    lOper(1) = 1
    DiffOp = .true.
    if (iCOSMO > 0) then
      call fzero(Temp,ngrad)
      Label = ' The Electronic Reaction Field Contribution (COSMO)'
      call OneEl_g(COSGrd,PCMMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
      if (iPrint >= 15) then
        Label = ' Reaction Field (COSMO) Contribution'
        call PrGrad(Label,Temp,nGrad,ChDisp)
      end if
    else
      Label = ' The Electronic Reaction Field Contribution (PCM)'
      call OneEl_g(PCMGrd,PCMMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
      if (iPrint >= 15) then
        Label = ' Reaction Field (PCM) Contribution'
        call PrGrad(Label,Temp,nGrad,ChDisp)
      end if
    end if

    Grad(:) = Grad(:)+Temp(:)
    if (iCOSMO == 0) PCM_SQ(:,:) = PCM_SQ(:,:)*real(nIrrep,kind=wp)

  end if
end if

!***********************************************************************
!8)                                                                    *
!     Trace the "variational" first order density matrix with the      *
!     gradient of the FAIEMP integrals.                                *
!                                                                      *
!***********************************************************************

if (.not. HF_Force) then
  if (lFAIEMP) then
    DiffOp = .true.
    Label = ' The FAIEMP Projection Operator Contribution'
    call OneEl_g(FragPGrd,FragPMmG,Temp,nGrad,DiffOp,Coor,D_Var,nDens,lOper,nComp,nOrdOp,Label)
    Grad(:) = Grad(:)+Temp(:)
    call DrvG_FAIEMP(Grad,Temp,nGrad)
  end if
end if

!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(lOper)
call mma_deallocate(Coor)
!                                                                      *
!***********************************************************************
!                                                                      *
! Epilogue, end

if (.not. HF_Force) call mma_deallocate(Fock)
call mma_deallocate(D_Var)

call CWTime(TCpu2,TWall2)
return

end subroutine Drvh1
