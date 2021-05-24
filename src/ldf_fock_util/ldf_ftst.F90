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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine LDF_FTst(UsePartPermSym,Mode,tau,nD,FactC,ip_DBlocks,ip_FBlocks)
! Thomas Bondo Pedersen, January 2012.
!
! Purpose: Compute Coulomb contribution to Fock matrix
!          using either exact or LDF integrals (using the formula
!          specified by Mode), depending on
!          positivity of the LDF integrals.
!          (Debug code.)

implicit none
logical UsePartPermSym
integer Mode
real*8 tau
integer nD
real*8 FactC(nD)
integer ip_DBlocks(nD)
integer ip_FBlocks(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

character*8 SecNam
parameter(SecNam='LDF_FTst')

integer PrintLevel
parameter(PrintLevel=2)

real*8 TolNeg
parameter(TolNeg=-1.0d-8)

integer LDF_nBas_Atom
external LDF_nBas_Atom

character*5 IntegralID

integer AB, CD
integer A, B, C, D
integer nAB, nCD
integer ip_Int, l_Int
integer iD
integer ipD, ipF

real*8 r, t

integer i, j
integer AP_Atoms
AP_Atoms(i,j) = iWork(ip_AP_Atoms-1+2*(j-1)+i)

r = 0.0d0
if (UsePartPermSym) then ! use particle permutation symmetry
  do AB=1,NumberOfAtomPairs
    A = AP_Atoms(1,AB)
    B = AP_Atoms(2,AB)
    nAB = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    do CD=1,AB-1
      C = AP_Atoms(1,CD)
      D = AP_Atoms(2,CD)
      nCD = LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
      l_Int = nAB*nCD
      call GetMem('FTstInt','Allo','Real',ip_Int,l_Int)
      call LDF_getIntegralsSelectedByPSD(PrintLevel,Mode,tau,TolNeg,AB,CD,l_Int,Work(ip_Int),IntegralID)
      if (IntegralID == 'exact') then
        r = r+1.0d0
      end if
      do iD=1,nD
        ipD = iWork(ip_DBlocks(iD)-1+CD)
        ipF = iWork(ip_FBlocks(iD)-1+AB)
        call dGeMV_('N',nAB,nCD,FactC(iD),Work(ip_Int),max(nAB,1),Work(ipD),1,1.0d0,Work(ipF),1)
      end do
      do iD=1,nD
        ipD = iWork(ip_DBlocks(iD)-1+AB)
        ipF = iWork(ip_FBlocks(iD)-1+CD)
        call dGeMV_('T',nAB,nCD,FactC(iD),Work(ip_Int),max(nAB,1),Work(ipD),1,1.0d0,Work(ipF),1)
      end do
      call GetMem('FTstInt','Free','Real',ip_Int,l_Int)
    end do
    l_Int = nAB**2
    call GetMem('FTstInt','Allo','Real',ip_Int,l_Int)
    call LDF_getIntegralsSelectedByPSD(PrintLevel,Mode,tau,TolNeg,AB,CD,l_Int,Work(ip_Int),IntegralID)
    if (IntegralID == 'exact') then
      r = r+1.0d0
    end if
    do iD=1,nD
      ipD = iWork(ip_DBlocks(iD)-1+AB)
      ipF = iWork(ip_FBlocks(iD)-1+AB)
      call dGeMV_('N',nAB,nAB,FactC(iD),Work(ip_Int),max(nAB,1),Work(ipD),1,1.0d0,Work(ipF),1)
    end do
    call GetMem('FTstInt','Free','Real',ip_Int,l_Int)
  end do
else
  do AB=1,NumberOfAtomPairs
    A = AP_Atoms(1,AB)
    B = AP_Atoms(2,AB)
    nAB = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    do CD=1,NumberOfAtomPairs
      C = AP_Atoms(1,CD)
      D = AP_Atoms(2,CD)
      nCD = LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
      l_Int = nAB*nCD
      call GetMem('FTstInt','Allo','Real',ip_Int,l_Int)
      call LDF_getIntegralsSelectedByPSD(PrintLevel,Mode,tau,TolNeg,AB,CD,l_Int,Work(ip_Int),IntegralID)
      if (IntegralID == 'exact') then
        r = r+1.0d0
      end if
      do iD=1,nD
        ipD = iWork(ip_DBlocks(iD)-1+CD)
        ipF = iWork(ip_FBlocks(iD)-1+AB)
        call dGeMV_('N',nAB,nCD,FactC(iD),Work(ip_Int),nAB,Work(ipD),1,1.0d0,Work(ipF),1)
      end do
      call GetMem('FTstInt','Free','Real',ip_Int,l_Int)
    end do
  end do
end if

write(6,'(A,/,80A)') SecNam,('=',iD=1,len(SecNam))
write(6,'(3X,A,I10)') 'LDF integral mode......................',Mode
write(6,'(3X,A,L2)') 'Particle permutation symmetry used.....',UsePartPermSym
if (NumberOfAtomPairs > 0) then
  t = dble(NumberOfAtomPairs)
  if (UsePartPermSym) then
    t = t*(t+1.0d0)/2.0d0
  else
    t = t*t
  end if
  t = 1.0d2*r/t
  write(6,'(3X,A,I10,1X,A,F7.2,A)') 'Number of exact integral blocks used...',int(r),'(',t,'%)'
end if
call xFlush(6)

end subroutine LDF_FTst
