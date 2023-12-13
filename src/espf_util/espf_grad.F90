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

subroutine espf_grad(natom,nGrdPt,nAtQM,Ext,Grid,B,DB,IsMM,GradCl,DoTinker,DoGromacs)
! Gradient due to the external potential

use espf_global, only: MxExtPotComp
use Index_Functions, only: nTri_Elem
use Disp, only: ChDisp, lDisp
use NAC, only: isNAC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Angstrom, auTokcalmol
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: natom, nGrdPt, nAtQM, IsMM(natom)
real(kind=wp), intent(in) :: Ext(MxExtPotComp,natom), Grid(3,nGrdPt), B(nGrdPt), DB(nGrdPt,3,nAtQM), GradCl(3,natom)
logical(kind=iwp), intent(in) :: DoTinker, DoGromacs
integer(kind=iwp) :: iAddPot, iAt, iAtom, iBla, iCur, iNumb, iPL, iPnt, iQM, iStep, ITkQMMM, iXYZ, jAtom, jPnt, jXYZ, nAtMM, ncmp, &
                     nData, nGrad
real(kind=wp) :: FX(4), opnuc(1)
logical(kind=iwp) :: Exists, isNAC_tmp, lMMGrd, lMMHess
character(len=180) :: Line
real(kind=wp), allocatable :: D2(:), Grad(:,:), GrdI(:,:), HC(:), MMGrd(:,:,:), Temp(:), XC(:)
integer(kind=iwp), external :: iPL_espf, IsFreeUnit, LHR
character(len=180), external :: Get_Ln

iPL = iPL_espf()

call mma_allocate(Grad,3,nAtom,Label='Grad')
nGrad = 3*nAtom
call Get_dArray_chk('GRAD',Grad,nGrad)
lMMHess = .false.
if (iPL >= 3) call PrGrad(' Molecular gradients, entering ESPF',Grad,lDisp(0),ChDisp)

! Recover MM gradient and hessian, if any, in QMMM file

call F_Inquire('QMMM',Exists)
nAtMM = 0
if (((Exists .and. DoTinker) .or. DoGromacs) .and. (.not. isNAC)) then
  call mma_allocate(MMGrd,3,natom,2,label='MMGrd')
  call Qpg_dArray('MM Grad',lMMGrd,nData)
  if (lMMGrd) then
    call Get_dArray('MM Grad',MMGrd,6*natom)
    MMGrd(:,:,1) = MMGrd(:,:,2)
    MMGrd(:,:,2) = Zero
  else
    MMGrd(:,:,:) = Zero
  end if
end if

if (Exists .and. DoTinker .and. (.not. isNAC)) then
  call mma_allocate(HC,nTri_Elem(3*natom),label='Hess')
  HC(:) = Zero
end if

if (Exists .and. DoTinker .and. (.not. isNAC)) then
  ITkQMMM = IsFreeUnit(15)
  call Molcas_Open(ITkQMMM,'QMMM')
  Line = ' '
  do while (index(Line,'TheEnd ') == 0)
    Line = Get_Ln(ITkQMMM)
    if (index(Line,'NMM') /= 0) then
      call Get_I1(2,nAtMM)
    else if (index(Line,'MMGradient') /= 0) then
      call Get_I1(2,iAtom)
      call Get_F(3,FX,3)
      MMGrd(:,iAtom,2) = FX(1:3)*Angstrom/auTokcalmol
      Grad(:,iAtom) = Grad(:,iAtom)+MMGrd(:,iAtom,2)
    else if (index(Line,'MMHDiag') /= 0) then
      lMMHess = .true.
      call Get_I1(2,iAtom)
      call Get_F(3,FX,3)
      do iXYZ=1,3
        HC(LHR(iXYZ,iAtom,iXYZ,iAtom)) = FX(iXYZ)*Angstrom**2/auTokcalmol
      end do
    else if (index(Line,'MMHOff') /= 0) then
      lMMHess = .true.
      call Get_I1(2,iAtom)
      call Get_I1(3,iXYZ)
      call Get_I1(4,iNumb)
      !write(u6,*) 'HOff read ',iAtom,iXYZ,iNumb
      iCur = 0
      jAtom = iAtom
      jXYZ = iXYZ
      do
        iStep = min(4,iNumb-iCur)
        Line = Get_Ln(ITkQMMM)
        call Get_F(1,FX,iStep)
        !write(u6,'(A,4f10.5)') 'HOff read ',(FX(j),j=1,iStep)
        do iBla=1,iStep
          jXYZ = jXYZ+1
          if (jXYZ == 4) then
            jXYZ = 1
            jAtom = jAtom+1
            if (jAtom > natom) call Abend()
          end if
          HC(LHR(iXYZ,iAtom,jXYZ,jAtom)) = FX(iBla)*Angstrom**2/auTokcalmol
        end do
        iCur = iCur+4
        if (iCur >= iNumb) exit
      end do
    end if
  end do
  close(ITkQMMM)
end if

if (DoGromacs .and. (.not. isNAC)) then
  MMGrd(:,:,2) = GradCl
  Grad(:,:) = Grad(:,:)+MMGrd(:,:,2)
end if

if (((Exists .and. DoTinker) .or. DoGromacs) .and. (.not. isNAC)) then
  if (iPL >= 4) then
    call RecPrt('Old MM Grad:',' ',MMGrd(:,:,1),3,natom)
    call RecPrt('New MM Grad:',' ',MMGrd(:,:,2),3,natom)
  end if
  if (nAtMM > 0) call Put_dArray('MM Grad',MMGrd,6*natom)
  call mma_deallocate(MMGrd)
end if

if (Exists .and. DoTinker .and. (.not. isNAC)) then
  if (lMMHess .and. iPL >= 4) call TriPrt(' In ESPF_grad: MM Hessian','(12f12.7)',HC,3*natom)
  if (lMMHess) call Put_dArray('MMHessian',HC,nTri_Elem(3*natom))
  call mma_deallocate(HC)
end if

if (((Exists .and. DoTinker) .or. DoGromacs) .and. (.not. isNAC)) then
  call Put_iScalar('No of Internal coordinates',3*natom)
  if (iPL >= 3) call PrGrad(' Molecular gradients, after MM',Grad,lDisp(0),ChDisp)
end if

! External field acting on nuclear charges

if (isNac) then
  write(u6,*) 'ESPF: Skipping nuclear-external field contribution'
else
  call mma_allocate(XC,nAtom,label='XCharge')
  call Get_dArray('Effective nuclear Charge',XC,nAtom)
  do iAt=1,nAtom
    Grad(:,iAt) = Grad(:,iAt)+XC(iAt)*Ext(2:4,iAt)
  end do
  call mma_deallocate(XC)
  if (iPL >= 3) call PrGrad(' Molecular grad, after nuc ESPF',Grad,lDisp(0),ChDisp)
end if

! Here I need the integral derivatives, weighted by B and contracted
! with the density matrix: B * Sum_mu,nu P_mu,nu*d/dq(<mu|1/R_grid|nu>)
!
! Basically, it should work like the following lines ... but it can't now.
!
!  opnuc = Zero
!  ncmp = 3
!  iAddPot = -1
!  call mma_allocate(D1,nGrdPt,3*natom,label='D1')
!  call DrvPot(Grid,opnuc,ncmp,D1,nGrdPt,iAddPot)
!  call mma_deallocate(D1)
!
! Now I try another thing, following the way derivatives with respect to point
! charges are computed in alaska. I just copied 2 files from alaska: drvh1 and
! xfdgrd then renamed them drvespf and bdvgrd.

call mma_allocate(Temp,3*natom,label='Temp')
call mma_allocate(GrdI,4,nGrdPt,label='GridInfo')
! Need to save isNAC here because Prepare calling inisew calling init_seward which resets isNAC .to False.
isNAC_tmp = isNAC
call Prepare(nGrdPt,Grid,B,GrdI)
call Drvespf(Grad,Temp,3*natom,GrdI)
call mma_deallocate(GrdI)
if (iPL >= 3) call PrGrad(' Molecular gradients, after P*B*dV',Grad,lDisp(0),ChDisp)
call mma_deallocate(Temp)

! Here I need the integrals contracted with the density matrix and weighted
! by the derivatives of B: dB/dq * Sum_mu,nu P_mu,nu*(<mu|1/R_grid|nu>)

opnuc = Zero
ncmp = 1
iAddPot = -1
call mma_allocate(D2,nGrdPt,label='D2')
call DrvPot(Grid,opnuc,ncmp,D2,nGrdPt,iAddPot)
if (iPL >= 4) then
  write(u6,'(/,A,/)') ' PV = '
  do iPnt=1,nGrdPt
    write(u6,*) D2(iPnt)
  end do
end if
iQM = 0
do iAt=1,natom
  if (IsMM(iAt) == 1) cycle
  iQM = iQM+1
  do jPnt=1,NGrdPt
    Grad(:,iAt) = Grad(:,iAt)+DB(jPnt,:,iQM)*D2(jPnt)
  end do
end do
isNAC = isNAC_tmp

! Apply Morokuma's scheme if needed

if ((Exists .and. DoTinker) .or. DoGromacs) call LA_Morok(natom,Grad,1)

! Finally

call Put_dArray('GRAD',Grad,3*natom)
call mma_deallocate(D2)
if (iPL >= 2) call PrGrad(' Molecular gradients, after ESPF',Grad,lDisp(0),ChDisp)
call Add_Info('Grad',Grad,3*natom,6)
call mma_deallocate(Grad)

return

end subroutine espf_grad
