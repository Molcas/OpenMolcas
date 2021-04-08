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
! Copyright (C) 1991, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Rd1Int_Motra(ipOvlp,ipHOne,ipKine)
!***********************************************************************
!                                                                      *
! Objective: Read the header of the one-electron integral file         *
!            Extract also symmetry and basis set information.          *
!            In addition read the overlap, the nuclear attraction      *
!            and kinetic integrals.                                    *
!                                                                      *
!**** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************

use motra_global, only: BsLbl, Header, iRFpert, n2max, nBas, nSym, nTot1, nTot2, PotNuc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ipOvlp, ipHOne, ipKine
#include "mxdm.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iBas, iComp, iOpt, iRc, iSyLbl, iSym, nDim, nTemp
real(kind=wp) :: ERFself
character(len=8) :: OneLbl
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: Temp(:)

!----------------------------------------------------------------------*
! Read one-electron integral file header etc.                          *
!----------------------------------------------------------------------*
call Get_cArray('Seward Title',Header,144)
!----------------------------------------------------------------------*
! Read no.of symm. species                                             *
!----------------------------------------------------------------------*
call Get_iScalar('nSym',nSym)
!----------------------------------------------------------------------*
! Read symm. oper per symm. species                                    *
!----------------------------------------------------------------------*
!call Get_iArray('Symmetry operations',iOper,nSym)
!----------------------------------------------------------------------*
! Read no. of basis functions per symm. species                        *
!----------------------------------------------------------------------*
call Get_iArray('nBas',nBas,nSym)
!----------------------------------------------------------------------*
! Read no. of basis functions per symm. species                        *
!----------------------------------------------------------------------*
nDim = 0
do iSym=1,nSym
  nDim = nDim+nBas(iSym)
end do
call mma_allocate(BsLbl,mxOrb,label='BsLbl') ! allocate much more than needed, it will be assumed by other programs
call Get_cArray('Unique Basis Names',BsLbl,len(BsLbl)*nDim)
BsLbl(nDim+1:) = ''
!----------------------------------------------------------------------*
! Read no. of unique atoms in the system                               *
!----------------------------------------------------------------------*
!call Get_iScalar('Unique atoms',nAtoms)
!----------------------------------------------------------------------*
! Read atom labels                                                     *
!----------------------------------------------------------------------*
!call Get_cArray('Unique Atom Names',AtLbl,LenIn*nAtoms)
!----------------------------------------------------------------------*
! Read coordinates of atoms                                            *
!----------------------------------------------------------------------*
!call Get_dArray('Unique Coordinates',Coor,3*nAtoms)
!----------------------------------------------------------------------*
! Read nuclear repulsion energy                                        *
!----------------------------------------------------------------------*
!call Get_PotNuc(PotNuc)
call Get_dScalar('PotNuc',PotNuc)
!----------------------------------------------------------------------*
! Allocate memory for one-electron integrals                           *
!----------------------------------------------------------------------*
nTot1 = 0
nTot2 = 0
n2max = 0
do iSym=1,nSym
  iBas = nBas(iSym)
  nTot1 = nTot1+iBAs*(iBas+1)/2
  nTot2 = nTot2+iBas*iBas
  n2max = max(n2max,iBas*iBas)
end do

call GetMem('Ovlp','Allo','Real',ipOvlp,nTot1+4)
call GetMem('Kine','Allo','Real',ipKine,nTot1+4)
call GetMem('HOne','Allo','Real',ipHOne,nTot1+4)
!----------------------------------------------------------------------*
! Read overlap integrals                                               *
!----------------------------------------------------------------------*
iRc = -1
iOpt = 6
iComp = 1
iSyLbl = 1
OneLbl = 'Mltpl  0'
call RdOne(iRc,iOpt,OneLbl,iComp,Work(ipOvlp),iSyLbl)

if (iRc /= 0) call Error()
!----------------------------------------------------------------------*
! Read core Hamiltonian                                                *
!----------------------------------------------------------------------*
iRc = -1
iOpt = 6
iComp = 1
iSyLbl = 1
OneLbl = 'OneHam  '
call RdOne(iRc,iOpt,OneLbl,iComp,Work(ipHone),iSyLbl)

if (iRc /= 0) call Error()
!----------------------------------------------------------------------*
! Read kinetic energy integrals                                        *
!----------------------------------------------------------------------*
iRc = -1
iOpt = 6
iComp = 1
iSyLbl = 1
OneLbl = 'Kinetic '
call RdOne(iRc,iOpt,OneLbl,iComp,Work(ipKine),iSyLbl)

if (iRc /= 0) call Error()
!----------------------------------------------------------------------*
! If this is a perturbative reaction field calculation then            *
! modifiy the one-electron Hamiltonian by the reaction field and       *
! the nuclear attraction by the cavity self-energy                     *
!----------------------------------------------------------------------*
if (iRFpert /= 0) then
  nTemp = 0
  do iSym=1,nSym
    nTemp = nTemp+nBas(iSym)*(nBas(iSym)+1)/2
  end do
  call mma_allocate(Temp,nTemp,label='RFFLD')
  call f_Inquire('RUNOLD',Found)
  if (Found) call NameRun('RUNOLD')
  call Get_dScalar('RF Self Energy',ERFself)
  call Get_dArray('Reaction field',Temp,nTemp)
  if (Found) call NameRun('RUNFILE')
  PotNuc = PotNuc+ERFself
  call Daxpy_(nTemp,One,Temp,1,Work(ipHone),1)
  call mma_deallocate(Temp)
end if
!----------------------------------------------------------------------*
! Normal termination                                                   *
!----------------------------------------------------------------------*
return

contains

subroutine Error()
  write(u6,*) 'Rd1Int: Error reading from ONEINT'
  write(u6,*) 'OneLbl=',OneLbl
  call Abend()
end subroutine Error

end subroutine Rd1Int_Motra
