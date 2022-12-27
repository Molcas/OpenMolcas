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

subroutine AverMEP(Kword,Eint,Poli,iCi,SumElcPot,NCountField,PertElcInt,iQ_Atoms,nBas,nOcc,natyp,nntyp)

use qmstat_global, only: AvElcPot, ChaNuc, FieldNuc, iPrint, MxMltp, nMlt, outxyz, PertNElcInt
use Index_Functions, only: iTri, nTri3_Elem, nTri_Elem
use Data_Structures, only: Alloc1DArray_Type, Allocate_DT, Deallocate_DT
use OneDat, only: sNoNuc, sNoOri, sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Three, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
character(len=4), intent(inout) :: Kword
integer(kind=iwp), intent(in) :: iCi, NCountField, iQ_Atoms, nBas, nntyp, nOcc(nntyp), natyp(nntyp)
real(kind=wp), intent(in) :: Eint(iCi,10), Poli(iCi,10)
real(kind=wp), intent(inout) :: SumElcPot(iCi,10)
real(kind=wp), intent(out) :: PertElcInt(nTri_Elem(nBas))
integer(kind=iwp) :: i, i1, i2, iB1, iB2, iComp, iiDum(1), iLuField, indMME, iOpt, irc, iSmLbl, iTyp, j, kaunta, Lu_One, nSize, nTyp
real(kind=wp) :: Tra
logical(kind=iwp) :: Exists
character(len=8) :: Label
integer(kind=iwp), allocatable :: Dum(:,:), iCent(:)
real(kind=wp), allocatable :: AvTemp(:), ForceNuc(:,:), H0(:), H1(:)
type(Alloc1DArray_Type), allocatable :: MME(:)
integer(kind=iwp), external :: IsFreeUnit
#include "warnings.h"

call UpCase(Kword)
!**************
! This subroutine includes three different options. All have to do with the
! calculation of a Mean Electrostatic Potential, Field and Field gradients,
! and to evalute the perturbation of them in the One electron Hamiltoniam
! this perturbation is added to the SEWARD One Electron File.
! Option 1: Add the components of Potential, etc
! Option 2: Obtain the Average
! Option 3: Calculate the electrostatic perturbation energy integrals
!           and add them to the one-electron file.
! Calculations involve up to the field gradients because the charge density
! is expanded to the quadrupoles. If the expansion is bigger the number 10
! must be changed, but also all the eqscf and eqras subroutines shoud be
! changed. In the last option the array nMlt is used instead of the
! number so if a smaller expantion is used, non problem,
! since this index takes care of that.
!*****************

select case (Kword(1:4))

  case default !('ADD ')
    do i=1,iCi
      do j=1,10 !Charges(1),Dipoles(3),Quadrupoles(6)
        SumElcPot(i,j) = SumElcPot(i,j)+Eint(i,j)+Poli(i,j)
      end do
    end do
    if (iPrint >= 9) then
      write(u6,*) 'Total Sum Potential'
      do i=1,iCi
        write(u6,*) SumElcPot(i,:)
      end do
    end if

  case ('AVER')
    AvElcPot(:,:) = SumElcPot/real(NCountField,kind=wp)

    ! Charges (1),Dipoles(3),Quadrupoles(6)

    ! The order of Field gradients is changed in order to follow the same order than Molcas
    ! This change is due to the different order of quadrupoles in QmStat and Molcas.
    ! QmStat:xx,xy,yy,xz,yz,zz Molcas:xx,xy,xz,yy,yz,zz
    call mma_allocate(AvTemp,iCi,label='AvTemp')
    AvTemp(:) = AvElcPot(:,8)
    AvElcPot(:,8) = AvElcPot(:,7)
    AvElcPot(:,8) = AvTemp
    call mma_deallocate(AvTemp)

    !*******************************
    ! This multiplication comes because the off-diagonal
    ! quadrupoles must be multiplied by two since we use
    ! a triangular form to compute the Interaction
    ! Energy with the Electric Field Gradient.
    ! Since it is easier to multiply the Average potential
    ! than the quadrupole for each pair of bases, we perform
    ! the multiplication here
    !**********************
    AvElcPot(:,6) = Two*AvElcPot(:,6)
    AvElcPot(:,7) = Two*AvElcPot(:,7)
    AvElcPot(:,9) = Two*AvElcPot(:,9)
    !**********************
    if (iPrint >= 9) then
      write(u6,*) 'Total Averg Potential'
      do i=1,iCi
        write(u6,*) AvElcPot(i,:)
      end do
    end if

  case ('PERT')
    ! First we read the multipoles expansion for each pair of bases.
    ! The index iCent(i) will give us to which center belongs each pair of bases.

    call mma_allocate(outxyz,3,nTri_Elem(iQ_Atoms),label='outxyz')
    call mma_allocate(iCent,nTri_Elem(nBas),label='iCent')
    call mma_allocate(Dum,nBas,nBas,label='Dummy')
    call Allocate_DT(MME,[1,nTri3_Elem(MxMltp)],label='MME')
    call MultiNew(iQ_Atoms,nBas,nOcc,natyp,nntyp,MME,iCent,Dum,nMlt,outxyz,.false.)
    call mma_deallocate(Dum)

    !*********************
    ! Calculate the forces for the nuclei these forces will compensate partially
    ! the forces due to the electrons, They will be printed and added to the
    ! RUNFILE in the optimization procedure after Alaska module
    !********************
    ! This model does not work for calculating the forces in the nuclei
    ! with a Slater representation since you have to calculate the field
    ! in a set of point charges and not in distributed charges as the field
    ! is calculated when used Slater representation also there is a more
    ! dark and complicated problem about the string interaction keeping
    ! together the distributed electronic charge and the point nuclear
    !  charge under different forces.
    !*********************
    call mma_allocate(ForceNuc,3,iQ_Atoms,label='ForceNuc')
    do i=1,iQ_Atoms
      ForceNuc(:,i) = ChaNuc(i)*AvElcPot(i,2:4)
    end do
    iLuField = IsFreeUnit(63)
    call OpnFl(FieldNuc,iLuField,Exists)
    write(u6,*) 'FieldNuc',FieldNuc
    do i=1,iQ_Atoms
      write(iLuField,*) ForceNuc(:,i)
    end do
    close(iLuField)

    if (iPrint >= 9) then
      write(u6,*) 'Nuclei charge and Forces'
      do i=1,iQ_Atoms
        write(u6,*) ChaNuc(i),ForceNuc(:,i)
      end do
    end if
    call mma_deallocate(ForceNuc)
    !********************

    nTyp = 0
    do i=1,nMlt
      nTyp = nTyp+nTri_Elem(i)
    end do
    PertElcInt(:) = Zero

    ! Put quadrupoles in Buckingham form.

    do i1=1,nBas
      do i2=1,i1
        indMME = iTri(i1,i2)
        do j=5,10
          MME(j)%A(indMME) = MME(j)%A(indMME)*OneHalf
        end do
        Tra = (MME(5)%A(indMME)+MME(8)%A(indMME)+MME(10)%A(indMME))/Three
        MME(5)%A(indMME) = MME(5)%A(indMME)-Tra
        MME(8)%A(indMME) = MME(8)%A(indMME)-Tra
        MME(10)%A(indMME) = MME(10)%A(indMME)-Tra
      end do
    end do

    irc = -1
    Lu_One = IsFreeUnit(49)
    iOpt = 0
    call OpnOne(irc,iOpt,'ONEINT',Lu_One)
    if (irc /= 0) then
      write(u6,*)
      write(u6,*) 'ERROR! Could not open one-electron integral file.'
      call Quit(_RC_IO_ERROR_READ_)
    end if

    ! We read the size of the unperturbed Hamiltonian 'OneHam 0' in OneInt.

    irc = -1
    iOpt = ibset(0,sOpSiz)
    iSmLbl = 1
    nSize = 0
    Label = 'OneHam 0'
    iComp = 1
    call iRdOne(irc,iOpt,Label,iComp,iiDum,iSmLbl)
    nSize = iiDum(1)
    if (irc /= 0) then
      write(u6,*)
      write(u6,*) 'ERROR! Failed to read number of one-electron integrals.'
      call Quit(_RC_IO_ERROR_READ_)
    end if
    if (nSize == 0) then
      write(u6,*)
      write(u6,*) 'ERROR! Problem reading size of unperturbed Hamiltonian in OneInt'
      call Quit(_RC_IO_ERROR_READ_)
    end if

    ! Memory allocation for the unperturbed Hamiltonian
    call mma_allocate(H0,nSize,label='MAver')
    irc = -1
    iOpt = ibset(ibset(0,sNoOri),sNoNuc)
    iSmLbl = 0

    ! Read the unperturbed Hamiltonian
    call RdOne(irc,iOpt,Label,iComp,H0,iSmLbl) !Collect non perturbed integrals
    call mma_allocate(H1,nSize,label='MAver1')
    if (iPrint >= 9) call TriPrt('Non Perturb One-e',' ',H0,nBas)

    ! We perform the multiplication for each pair of bases in a triangular form.
    ! The perturbation is added to the unperturbed Hamiltonian 'H0'.

    kaunta = 0
    do iB1=1,nBas
      do iB2=1,iB1
        kaunta = kaunta+1
        indMME = iTri(iB1,iB2)
        do iTyp=1,nTyp
          PertElcInt(indMME) = PertElcInt(indMME)+AvElcPot(iCent(kaunta),iTyp)*MME(iTyp)%A(indMME)
        end do
        H1(kaunta) = H0(kaunta)+PertElcInt(indMME)
      end do
    end do
    call mma_deallocate(iCent)
    call Deallocate_DT(MME)

    if (iPrint >= 9) call TriPrt('H0+Elec One-e',' ',H1,nBas)

    ! The non-Electrostatic perturbation is added. The PertNElcInt array comes
    ! through the module qmstat_global

    if (iPrint >= 10) call TriPrt('PertNElcInt-e',' ',PertNElcInt,nBas)

    H1(:) = H1+PertNElcInt

    if (iPrint >= 9) call TriPrt('H0+Elec+nonEl One-e',' ',H1,nBas)

    ! The perturbed Hamiltonian 'H1' is writen in OneInt.
    irc = -1
    iOpt = 0
    iSmLbl = 1
    call WrOne(irc,iOpt,'OneHam  ',1,H1,iSmLbl) !Write perturbed integrals
    if (iPrint >= 9) call TriPrt('Perturb One-e',' ',H1,nBas)

    if (iPrint >= 10) call TriPrt('Non Perturb One-e AGAIN',' ',H0,nBas)

    call ClsOne(irc,Lu_One)

    call mma_deallocate(H0)
    call mma_deallocate(H1)

end select

return

end subroutine AverMEP
