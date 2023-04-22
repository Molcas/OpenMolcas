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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: UsePartPermSym
integer(kind=iwp), intent(in) :: Mode, nD, ip_DBlocks(nD), ip_FBlocks(nD)
real(kind=wp), intent(in) :: tau, FactC(nD)
character(len=5) :: IntegralID
integer(kind=iwp) :: AB, CD, A, B, C, D, nAB, nCD, l_Int, iD, ipD, ipF
real(kind=wp) :: r, t
real(kind=wp), allocatable :: FTstInt(:)
character(len=*), parameter :: SecNam = 'LDF_FTst'
integer(kind=iwp), parameter :: PrintLevel = 2
real(kind=wp), parameter :: TolNeg = -1.0e-8_wp
integer(kind=iwp), external :: LDF_nBas_Atom
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

r = Zero
if (UsePartPermSym) then ! use particle permutation symmetry
  do AB=1,NumberOfAtomPairs
    A = iWork(ip_AP_Atoms-1+2*(AB-1)+1)
    B = iWork(ip_AP_Atoms-1+2*(AB-1)+2)
    nAB = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    do CD=1,AB-1
      C = iWork(ip_AP_Atoms-1+2*(CD-1)+1)
      D = iWork(ip_AP_Atoms-1+2*(CD-1)+2)
      nCD = LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
      l_Int = nAB*nCD
      call mma_allocate(FTstInt,l_Int,label='FTstInt')
      call LDF_getIntegralsSelectedByPSD(PrintLevel,Mode,tau,TolNeg,AB,CD,l_Int,FTstInt,IntegralID)
      if (IntegralID == 'exact') then
        r = r+One
      end if
      do iD=1,nD
        ipD = iWork(ip_DBlocks(iD)-1+CD)
        ipF = iWork(ip_FBlocks(iD)-1+AB)
        call dGeMV_('N',nAB,nCD,FactC(iD),FTstInt,max(nAB,1),Work(ipD),1,One,Work(ipF),1)
      end do
      do iD=1,nD
        ipD = iWork(ip_DBlocks(iD)-1+AB)
        ipF = iWork(ip_FBlocks(iD)-1+CD)
        call dGeMV_('T',nAB,nCD,FactC(iD),FTstInt,max(nAB,1),Work(ipD),1,One,Work(ipF),1)
      end do
      call mma_deallocate(FTstInt)
    end do
    l_Int = nAB**2
    call mma_allocate(FTstInt,l_Int,label='FTstInt')
    call LDF_getIntegralsSelectedByPSD(PrintLevel,Mode,tau,TolNeg,AB,CD,l_Int,FTstInt,IntegralID)
    if (IntegralID == 'exact') then
      r = r+One
    end if
    do iD=1,nD
      ipD = iWork(ip_DBlocks(iD)-1+AB)
      ipF = iWork(ip_FBlocks(iD)-1+AB)
      call dGeMV_('N',nAB,nAB,FactC(iD),FTstInt,max(nAB,1),Work(ipD),1,One,Work(ipF),1)
    end do
    call mma_deallocate(FTstInt)
  end do
else
  do AB=1,NumberOfAtomPairs
    A = iWork(ip_AP_Atoms-1+2*(AB-1)+1)
    B = iWork(ip_AP_Atoms-1+2*(AB-1)+2)
    nAB = LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
    do CD=1,NumberOfAtomPairs
      C = iWork(ip_AP_Atoms-1+2*(CD-1)+1)
      D = iWork(ip_AP_Atoms-1+2*(CD-1)+2)
      nCD = LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
      l_Int = nAB*nCD
      call mma_allocate(FTstInt,l_Int,label='FTstInt')
      call LDF_getIntegralsSelectedByPSD(PrintLevel,Mode,tau,TolNeg,AB,CD,l_Int,FTstInt,IntegralID)
      if (IntegralID == 'exact') then
        r = r+One
      end if
      do iD=1,nD
        ipD = iWork(ip_DBlocks(iD)-1+CD)
        ipF = iWork(ip_FBlocks(iD)-1+AB)
        call dGeMV_('N',nAB,nCD,FactC(iD),FTstInt,nAB,Work(ipD),1,One,Work(ipF),1)
      end do
      call mma_deallocate(FTstInt)
    end do
  end do
end if

write(u6,'(A,/,A)') SecNam,repeat('=',len(SecNam))
write(u6,'(3X,A,I10)') 'LDF integral mode......................',Mode
write(u6,'(3X,A,L2)') 'Particle permutation symmetry used.....',UsePartPermSym
if (NumberOfAtomPairs > 0) then
  t = real(NumberOfAtomPairs,kind=wp)
  if (UsePartPermSym) then
    t = t*(t+One)*Half
  else
    t = t*t
  end if
  t = 1.0e2_wp*r/t
  write(u6,'(3X,A,I10,1X,A,F7.2,A)') 'Number of exact integral blocks used...',int(r),'(',t,'%)'
end if
call xFlush(u6)

end subroutine LDF_FTst
