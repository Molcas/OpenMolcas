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

subroutine AverMEP(Kword,Eint,Poli,ici,SumElcPot,NCountField,PertElcInt,iQ_Atoms,nBas,nOcc,natyp,nntyp)

use Index_Functions, only: nTri3_Elem
use Constants, only: Zero, One, Two, OneHalf
use Definitions, only: wp, iwp, u6

implicit none
#include "maxi.fh"
#include "qminp.fh"
#include "qm1.fh"
#include "qmcom.fh"
#include "files_qmstat.fh"
#include "WrkSpc.fh"
#include "warnings.h"
character(len=4) :: Kword
real(kind=wp) :: Eint(MxQCen,10), Poli(MxQCen,10), SumElcPot(MxQCen,10), PertElcInt(MxBas*(MxBas+1)/2)
integer(kind=iwp) :: ici, NCountField, iQ_Atoms, nBas, nOcc(*), natyp(*), nntyp
integer(kind=iwp) :: i, i1, i2, iB1, iB2, iCent(MxBas**2), iDum, iH0, iH1, iiDum(1), iLuField, iMME(nTri3_Elem(MxMltp)), indMME, & !IFG
                     iOpt, irc, iSmLbl, iTriBasQ, iTyp, j, kaunta, Lu_One, nSize, nTyp
real(kind=wp) :: AvTemp, ForceNuc(MxAt,3), SumOld(MxQCen,10), Tra !IFG
logical(kind=iwp) :: Exists
character(len=20) :: MemLab, MemLab1
integer(kind=iwp), external :: IsFreeUnit

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
      do j=1,10 !Charges (1),Dipoles(3),Quadrupoles(6)
        SumOld(i,j) = SumElcPot(i,j)
        SumElcPot(i,j) = SumOld(i,j)+Eint(i,j)+Poli(i,j)
      end do
    end do
    if (iPrint >= 9) then
      write(u6,*) 'Total Sum Potential'
      do i=1,iCi
        write(u6,*) (SumElcPot(i,j),j=1,10)
      end do
    end if

  case ('AVER')
    do i=1,iCi
      do j=1,10 !Charges (1),Dipoles(3),Quadrupoles(6)
        AvElcPot(i,j) = SumElcPot(i,j)/real(NCountField,kind=wp)
      end do

      ! The order of Field gradients is changed in order to follow
      ! the same order than Molcas

      AvTemp = AvElcPot(i,8)         ! This change is due to the different order of quadrupoles in QmStat and Molcas.
      AvElcPot(i,8) = AvElcPot(i,7)  ! QmStat:xx,xy,yy,xz,yz,zz Molcas:xx,xy,xz,yy,yz,zz
      AvElcPot(i,7) = AvTemp
    end do

    !*******************************
    ! This multiplication comes because the off-diagonal
    ! quadrupoles must be multiplied by two since we use
    ! a triangular form to compute the Interaction
    ! Energy with the Electric Field Gradient.
    ! Since it is easier to multiply the Average potential
    ! than the quadrupole for each pair of bases, we perform
    ! the multiplication here
    !**********************
    AvElcPot(i,6) = Two*AvElcPot(i,6)
    AvElcPot(i,7) = Two*AvElcPot(i,7)
    AvElcPot(i,9) = Two*AvElcPot(i,9)
    !**********************
    if (iPrint >= 9) then
      write(u6,*) 'Total Averg Potential'
      do i=1,iCi
        write(u6,*) (AvElcPot(i,j),j=1,10)
      end do
    end if

  case ('PERT')
    ! First we read the multipoles expansion for each pair of bases.
    ! The index iCent(i) will give us to which center belongs each pair of bases.

    call GetMem('Dummy','Allo','Inte',iDum,nBas**2)
    call MultiNew(iQ_Atoms,nBas,nOcc,natyp,nntyp,iMME,iCent,iWork(iDum),nMlt,outxyz,SlExpQ,.false.)
    call GetMem('Dummy','Free','Inte',iDum,nBas**2)

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
    do i=1,iQ_Atoms
      do j=1,3
        ForceNuc(i,j) = ChaNuc(i)*AvElcPot(i,j+1)
      end do
    end do
    iLuField = 63
    iLuField = IsFreeUnit(iLuField)
    call OpnFl(FieldNuc,iLuField,Exists)
    write(u6,*) 'FieldNuc',FieldNuc
    do i=1,iQ_Atoms
      write(iLuField,*) (ForceNuc(i,j),j=1,3)
    end do
    close(iLuField)

    if (iPrint >= 9) then
      write(u6,*) 'Nuclei charge and Forces'
      do i=1,iQ_Atoms
        write(u6,*) ChaNuc(i),(ForceNuc(i,j),j=1,3)
      end do
    end if
    !********************

    nTyp = 0
    do i=1,nMlt
      nTyp = nTyp+i*(i+1)/2
    end do
    do i=1,(nBas*(nBas+1)/2)
      PertElcInt(i) = Zero
    end do

    ! Put quadrupoles in Buckingham form.

    do i1=1,nBas
      do i2=1,i1
        indMME = i2+i1*(i1-1)/2
        do j=5,10
          Work(iMME(j)+indMME-1) = Work(iMME(j)+indMME-1)*OneHalf
        end do
        Tra = Work(iMME(5)+indMME-1)+Work(iMME(8)+indMME-1)+Work(iMME(10)+indMME-1)
        Tra = Tra/3
        Work(iMME(5)+indMME-1) = Work(iMME(5)+indMME-1)-Tra
        Work(iMME(8)+indMME-1) = Work(iMME(8)+indMME-1)-Tra
        Work(iMME(10)+indMME-1) = Work(iMME(10)+indMME-1)-Tra
      end do
    end do

    irc = -1
    Lu_One = 49
    Lu_One = IsFreeUnit(Lu_One)
    call OpnOne(irc,0,'ONEINT',Lu_One)
    if (irc /= 0) then
      write(u6,*)
      write(u6,*) 'ERROR! Could not open one-electron integral file.'
      call Quit(_RC_IO_ERROR_READ_)
    end if

    ! We read the size of the unperturbed Hamiltonian 'OneHam 0' in OneInt.

    irc = -1
    iOpt = 1
    iSmLbl = 1
    nSize = 0
    call iRdOne(irc,iOpt,'OneHam 0',1,iiDum,iSmLbl)
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
    write(MemLab,*) 'MAver'
    call GetMem(MemLab,'Allo','Real',iH0,nSize+4)
    irc = -1
    iOpt = 0
    iSmLbl = 0

    ! Read the unperturbed Hamiltonian
    call RdOne(irc,iOpt,'OneHam 0',1,Work(iH0),iSmLbl) !Collect non perturbed integrals
    write(MemLab1,*) 'MAver1'
    call GetMem(MemLab1,'Allo','Real',iH1,nSize+4)
    if (iPrint >= 9) then
      call TriPrt('Non Perturb One-e',' ',Work(iH0),nBas)
    end if

    ! We perform the multiplication for each pair of bases in a triangular form.
    ! The perturbation is added to the unperturbed Hamiltonian 'iH0'.

    kaunta = 0
    do iB1=1,nBas
      do iB2=1,iB1
        kaunta = kaunta+1
        indMME = iB2+iB1*(iB1-1)/2
        do iTyp=1,nTyp
          PertElcInt(indMME) = PertElcInt(indMME)+AvElcPot(iCent(kaunta),iTyp)*Work(iMME(iTyp)+indMME-1)
        end do
        Work(iH1+kaunta-1) = Work(iH0+kaunta-1)+PertElcInt(indMME)
      end do
    end do

    if (iPrint >= 9) then
      call TriPrt('H0+Elec One-e',' ',Work(iH1),nBas)
    end if

    ! The non-Electrostatic perturbation is added. The PertNElcInt array comes
    ! through the include file qminp.fh.

    if (iPrint >= 10) then
      call TriPrt('PertNElcInt-e',' ',PertNElcInt,nBas)
    end if

    iTriBasQ = nBas*(nBas+1)/2
    call DaxPy_(iTriBasQ,One,PertNElcInt,1,Work(iH1),1)

    if (iPrint >= 9) then
      call TriPrt('H0+Elec+nonEl One-e',' ',Work(iH1),nBas)
    end if

    ! The perturbed Hamiltonian 'H1' is writen in OneInt.
    irc = -1
    iOpt = 0
    iSmLbl = 1
    call WrOne(irc,iOpt,'OneHam  ',1,Work(iH1),iSmLbl) !Write perturbed integrals
    if (iPrint >= 9) then
      call TriPrt('Perturb One-e',' ',Work(iH1),nBas)
    end if

    if (iPrint >= 10) then
      call TriPrt('Non Perturb One-e AGAIN',' ',Work(iH0),nBas)
    end if

    call ClsOne(irc,Lu_One)

    call GetMem(MemLab,'Free','Real',iH0,nSize+4)
    call GetMem(MemLab1,'Free','Real',iH1,nSize+4)

end select

return

end subroutine AverMEP
