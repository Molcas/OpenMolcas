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

subroutine espf_grad(natom,nGrdPt,ipExt,ipGrid,ipB,ipDB,ipIsMM,ipGradCl,DoTinker,DoGromacs)
! Gradient due to the external potential

implicit real*8(A-H,O-Z)
#include "espf.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "nac.fh"
logical Exist, lMMHess, lMMGrd, DoTinker, DoGromacs, isNAC_tmp
character*180 Line
character*180 Get_Ln
external Get_Ln
dimension FX(4)
dimension opnuc(1)
real*8, allocatable :: Grad(:,:)

iPL = iPL_espf()

call mma_allocate(Grad,3,nAtom,Label='Grad')
nGrad = 3*nAtom
call Get_dArray_chk('GRAD',Grad,nGrad)
lMMHess = .false.
if (iPL >= 3) call PrGrad(' Molecular gradients, entering ESPF',Grad,lDisp(0),ChDisp)

! Recover MM gradient and hessian, if any, in QMMM file

call F_Inquire('QMMM',Exist)
natMM = 0
if (((Exist .and. DoTinker) .or. DoGromacs) .and. (.not. isNAC)) then
  call Allocate_Work(ipMMGrd,6*natom)
  call Qpg_dArray('MM Grad',lMMGrd,nData)
  if (lMMGrd) then
    call Get_dArray('MM Grad',Work(ipMMGrd),6*natom)
    call dcopy_(3*natom,Work(ipMMGrd+3*natom),1,Work(ipMMGrd),1)
    call dcopy_(3*natom,[Zero],0,Work(ipMMGrd+3*natom),1)
  else
    call dcopy_(6*natom,[Zero],0,Work(ipMMGrd),1)
  end if
end if

if (Exist .and. DoTinker .and. (.not. isNAC)) then
  call GetMem('Hess','Allo','Real',ipHC,3*natom*(3*natom+1)/2)
  call dcopy_((3*natom*(3*natom+1)/2),[Zero],0,Work(ipHC),1)
end if

if (Exist .and. DoTinker .and. (.not. isNAC)) then
  ITkQMMM = IsFreeUnit(15)
  call Molcas_Open(ITkQMMM,'QMMM')
  Line = ' '
  do while (index(Line,'TheEnd ') == 0)
    Line = Get_Ln(ITkQMMM)
    if (index(Line,'NMM') /= 0) then
      call Get_I1(2,natMM)
    else if (index(Line,'MMGradient') /= 0) then
      call Get_I1(2,iAtom)
      call Get_F(3,FX,3)
      do iXYZ=0,2
        iOff = (iAtom-1)*3+iXYZ
        Work(ipMMGrd+3*natom+iOff) = FX(iXYZ+1)*Angstrom*ToHartree
        Grad(iXYZ+1,iAtom) = Grad(iXYZ+1,iAtom)+Work(ipMMGrd+3*natom+iOff)
      end do
    else if (index(Line,'MMHDiag') /= 0) then
      lMMHess = .true.
      call Get_I1(2,iAtom)
      call Get_F(3,FX,3)
      do iXYZ=1,3
        Work(ipHC+LHR(iXYZ,iAtom,iXYZ,iAtom)-1) = FX(iXYZ)*Angstrom*Angstrom*ToHartree
      end do
    else if (index(Line,'MMHOff') /= 0) then
      lMMHess = .true.
      call Get_I1(2,iAtom)
      call Get_I1(3,iXYZ)
      call Get_I1(4,iNumb)
      !write(6,*) 'HOff read ',iAtom,iXYZ,iNumb
      iCur = 0
      jAtom = iAtom
      jXYZ = iXYZ
      do
        iStep = min(4,iNumb-iCur)
        Line = Get_Ln(ITkQMMM)
        call Get_F(1,FX,iStep)
        !write(6,'(A,4f10.5)') 'HOff read ',(FX(j),j=1,iStep)
        do iBla=1,iStep
          jXYZ = jXYZ+1
          if (jXYZ == 4) then
            jXYZ = 1
            jAtom = jAtom+1
            if (jAtom > natom) call Abend()
          end if
          Work(ipHC+LHR(iXYZ,iAtom,jXYZ,jAtom)-1) = FX(iBla)*Angstrom*Angstrom*ToHartree
        end do
        iCur = iCur+4
        if (iCur >= iNumb) exit
      end do
    end if
  end do
  close(ITkQMMM)
end if

if (DoGromacs .and. (.not. isNAC)) then
  call dcopy_(3*natom,Work(ipGradCl),1,Work(ipMMGrd+3*natom),1)
  do iAtom=1,nAtom
    do ixyz=1,3
      i = ixyz+3*(iAtom-1)-1
      Grad(ixyz,iAtom) = Grad(ixyz,iAtom)+Work(ipMMGrd+3*natom+i)
    end do
  end do
end if

if (((Exist .and. DoTinker) .or. DoGromacs) .and. (.not. isNAC)) then
  if (iPL >= 4) then
    call RecPrt('Old MM Grad:',' ',Work(ipMMGrd),3,natom)
    call RecPrt('New MM Grad:',' ',Work(ipMMGrd+3*natom),3,natom)
  end if
  if (natMM > 0) call Put_dArray('MM Grad',Work(ipMMGrd),6*natom)
  call Free_Work(ipMMGrd)
end if

if (Exist .and. DoTinker .and. (.not. isNAC)) then
  if (lMMHess .and. iPL >= 4) call TriPrt(' In ESPF_grad: MM Hessian','(12f12.7)',Work(ipHC),3*natom)
  if (lMMHess) call Put_dArray('MMHessian',Work(ipHC),3*natom*(3*natom+1)/2)
  call GetMem('Hess','Free','Real',ipHC,3*natom*(3*natom+1)/2)
end if

if (((Exist .and. DoTinker) .or. DoGromacs) .and. (.not. isNAC)) then
  call Put_iScalar('No of Internal coordinates',3*natom)
  if (iPL >= 3) call PrGrad(' Molecular gradients, after MM',Grad,lDisp(0),ChDisp)
end if

! External field acting on nuclear charges

if (isNac) then
  write(6,*) 'ESPF: Skipping nuclear-external field contribution'
else
  call GetMem('XCharge','Allo','Real',ipXC,nAtom)
  call Get_dArray('Effective nuclear Charge',Work(ipXC),nAtom)
  do iAt=1,nAtom
    iCurXC = ipXC+iAt-1
    iCurE = ipExt+(iAt-1)*MxExtPotComp
    Grad(1,iAt) = Grad(1,iAt)+Work(iCurXC)*Work(iCurE+1)
    Grad(2,iAt) = Grad(2,iAt)+Work(iCurXC)*Work(iCurE+2)
    Grad(3,iAt) = Grad(3,iAt)+Work(iCurXC)*Work(iCurE+3)
  end do
  call GetMem('XCharge','Free','Real',ipXC,natom)
  if (iPL >= 3) call PrGrad(' Molecular grad, after nuc ESPF',Grad,lDisp(0),ChDisp)
end if

! Here I need the integral derivatives, weighted by B and contracted
! with the density matrix: B * Sum_mu,nu P_mu,nu*d/dq(<mu|1/R_grid|nu>)
!
! Basically, it should work like the following lines ... but it can't now.
!
!  opnuc = Dum
!  ncmp = 3
!  iAddPot = -1
!  call GetMem('dESPF1','Allo','Real',ipD1,nGrdPt*natom*3)
!  call DrvPot(Work(ipGrid),opnuc,ncmp,Work(ipD1),nGrdPt,iAddPot)
!  call GetMem('dESPF1','Free','Real',ipD1,nGrdPt*natom*3)
!
! Now I try another thing, following the way derivatives with respect to point
! charges are computed in alaska. I just copied 2 files from alaska: drvh1 and
! xfdgrd then renamed them drvespf and bdvgrd.

call GetMem('Temp','Allo','Real',ipTemp,3*natom)
call GetMem('GridInfo','Allo','Real',ipGrdI,4*nGrdPt)
! Need to save isNAC here because Prepare calling inisew calling init_seward which resets isNAC .to False.
isNAC_tmp = isNAC
call Prepare(nGrdPt,ipGrid,ipB,ipGrdI)
call Drvespf(Grad,Work(ipTemp),3*natom,Work(ipGrdI))
call GetMem('GridInfo','Free','Real',ipGrdI,4*nGrdPt)
if (iPL >= 3) call PrGrad(' Molecular gradients, after P*B*dV',Grad,lDisp(0),ChDisp)
call GetMem('Temp','Free','Real',ipTemp,3*natom)

! Here I need the integrals contracted with the density matrix and weighted
! by the derivatives of B: dB/dq * Sum_mu,nu P_mu,nu*(<mu|1/R_grid|nu>)

opnuc = Dum
ncmp = 1
iAddPot = -1
call GetMem('dESPF2','Allo','Real',ipD2,nGrdPt)
call DrvPot(Work(ipGrid),opnuc,ncmp,Work(ipD2),nGrdPt,iAddPot)
if (iPL >= 4) then
  write(6,'(/,A,/)') ' PV = '
  do iPnt=1,nGrdPt
    write(6,*) Work(ipD2+iPnt-1)
  end do
end if
iQM = 0
do iAt=1,natom
  if (iWork(ipIsMM+iAt-1) == 1) cycle
  iQM = iQM+1
  do jPnt=1,NGrdPt
    iCurDB1 = ipDB+(jPnt-1)+((iQM-1)*3+0)*nGrdPt
    iCurDB2 = ipDB+(jPnt-1)+((iQM-1)*3+1)*nGrdPt
    iCurDB3 = ipDB+(jPnt-1)+((iQM-1)*3+2)*nGrdPt
    iCurI = ipD2+jPnt-1
    Grad(1,iAt) = Grad(1,iAt)+Work(iCurDB1)*Work(iCurI)
    Grad(2,iAt) = Grad(2,iAt)+Work(iCurDB2)*Work(iCurI)
    Grad(3,iAt) = Grad(3,iAt)+Work(iCurDB3)*Work(iCurI)
  end do
end do
isNAC = isNAC_tmp

! Apply Morokuma's scheme if needed

if ((Exist .and. DoTinker) .or. DoGromacs) call LA_Morok(natom,Grad,1)

! Finally

call Put_dArray('GRAD',Grad,3*natom)
call GetMem('dESPF2','Free','Real',ipD2,nGrdPt)
if (iPL >= 2) call PrGrad(' Molecular gradients, after ESPF',Grad,lDisp(0),ChDisp)
call Add_Info('Grad',Grad,3*natom,6)
call mma_deallocate(Grad)

return

end subroutine espf_grad
