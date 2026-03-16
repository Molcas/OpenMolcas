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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************

subroutine NTOSymAnalysis(NUseSym,NUseBF,NUsedBF,NTO,NTOType,STATENAME,EigVal,UsetoReal,RealtoUse,Spin,Sym,Ind,SumEigVal)

use Symmetry_Info, only: nSym => nIrrep
use Constants, only: Zero
use rassi_data, only: NBST, NASHT, NASH, NBASF, NISH, NSSH
use Definitions, only: u6

implicit none
! input variables
integer NUseSym
real*8 :: NTO(*)
integer, dimension(NSym) :: NUseBF, NUsedBF, UsetoReal, RealtoUse
character(len=8) NTOType
character(len=1) Spin
character(len=5) STATENAME
real*8 :: EigVal(*)
integer :: Sym(*), Ind(*)
character(len=128) FILENAME, molden_name
! Loop control
integer I, J, ICount
! Variables needed for judging the symmetry of a NTO
integer INTO, IUseSym, NNTO, ISym, IOrb
real*8, dimension(NUseSym) :: SquareSum
real*8, dimension(NBST) :: EigValArray
! SquareSum=Sum over square of coefficients for certain symmetry
integer, dimension(NUseSym) :: NOrbinSym
! Total number of orbitals in IUseSym
integer, dimension(NUseSym,NASHT) :: OrbSymIndex
! OrbSymIndex gives the original orbital index for a orbital in iusesym
real*8 SumEigVal
real*8, parameter :: Threshold = Zero
! If SquareSum(IUseSym) > Threshold, then print the coefficients in IUseSym symmetry
! If there are more than one symmetry with SquareSum(IUseSym) > Threshold,
! then give a warning message and print the one with the largest SquareSum
integer NPCMO, IPCMO
real*8, dimension(:), allocatable :: PCMO
! Printing control
integer iPrintSym, OrbNum, IOrbinSym
integer LU
real*8, dimension(2) :: vDum
integer, dimension(7,8) :: v2Dum
character(len=72) Note
integer, external :: ISFREEUNIT

do IUseSym=1,NUseSym
  NOrbinSym(IUseSym) = 0
end do

NPCMO = 0
do ISym=1,NSym
  NPCMO = NPCMO+NBASF(ISym)**2
end do
allocate(PCMO(NPCMO))

do INTO=1,NASHT
  iPrintSym = 0
  do IUseSym=1,NUseSym
    SquareSum(IUseSym) = Zero
    do ICount=1,NUseBF(IUseSym)
      I = INTO
      J = ICount+NUsedBF(IUseSym)
      SquareSum(IUseSym) = SquareSum(IUseSym)+NTO(I+(J-1)*NASHT)**2
    end do
    if (SquareSum(IUseSym) > Threshold) then
      if (iPrintSym == 0) then
        iPrintSym = IUseSym
      else
        write(u6,'(a)') 'There are at least two symmetries that have a sum of coefficient**2 larger than the threshold'
        write(u6,'(5A10)') 'Threshold','Sum1','Sum2','Sym1','Sym2'
        write(u6,'(3ES10.3E2,2I10)') Threshold,SquareSum(iPrintSym),SquareSum(IUseSym),iPrintSym,IUseSym
        if (SquareSum(iPrintSym) < SquareSum(IUseSym)) iPrintSym = IUseSym
      end if
    end if
  end do
  if (iPrintsym == 0) &
    write(u6,'(a,I2,a)') 'the symmetry of orbital ',INTO,' is not found. How is this possible? Change the value of DoTest in '// &
                         'ntocalc to .true. to print out the intermediate values'
  NOrbinSym(IPrintSym) = NOrbinSym(IPrintSym)+1
  OrbSymIndex(IPrintSym,NOrbinSym(IPrintSym)) = INTO
  Sym(INTO) = UsetoReal(IPrintSym)
  Ind(INTO) = NOrbinSym(IPrintSym)+NISH(UsetoReal(IPrintSym))
end do

! generating file in a similar way to other orbital files
IOrb = 0
IPCMO = 0
SumEigVal = Zero
do ISym=1,NSym
  IUseSym = RealtoUse(ISym)
  if (IUseSym /= 0) then
    ! If there are active orbitals in this symmetry
    NNTO = NOrbinSym(IUseSym)
    ! write inactive part
    do OrbNum=1,NISH(ISym)
      IOrb = IOrb+1
      EigValArray(IOrb) = Zero
    end do
    ! Recording Printed NTO (PCMO)
    do I=1,NISH(ISym)*NBASF(ISYM)
      IPCMO = IPCMO+1
      PCMO(IPCMO) = Zero
    end do
    ! Recording Printed NTO (PCMO)
    ! write active part
    do IOrbinSym=1,NNTO
      OrbNum = IOrbinSym+NISH(ISym)
      IOrb = IOrb+1
      I = OrbSymIndex(IUseSym,IOrbinSym)
      EigValArray(IOrb) = EigVal(I)
      SumEigVal = SumEigVal+EigVal(I)
      ! Recording Printed NTO (PCMO)
      do ICount=1,NUseBF(IUSeSym)
        IPCMO = IPCMO+1
        J = ICount+NUsedBF(IUseSym)
        PCMO(IPCMO) = NTO(I+(J-1)*NASHT)
      end do
      ! Recording Printed NTO (PCMO)
    end do
    ! write virtual part
    do OrBNum=NISH(ISym)+NASH(ISym)+1,NBASF(ISym)
      IOrb = IOrb+1
      EigValArray(IOrb) = Zero
    end do
    ! Recording Printed NTO (PCMO)
    do I=1,NSSH(ISym)*NBASF(ISYM)
      IPCMO = IPCMO+1
      PCMO(IPCMO) = Zero
    end do
    ! Recording Printed NTO (PCMO)
  else
    ! If there is no active orbitals in this symmetry
    do OrbNum=1,NBASF(ISym)
      IOrb = IOrb+1
      EigValArray(IOrb) = Zero
    end do
    ! Recording Printed NTO (PCMO)
    do I=1,NBASF(ISYM)**2
      IPCMO = IPCMO+1
      PCMO(IPCMO) = Zero
    end do
    ! Recording Printed NTO (PCMO)
  end if
end do

do I=1,NBST
  EigValArray(I) = EigValArray(I)/SumEigVal
end do

LU = ISFREEUNIT(50)
Note = '*  Natural Transition Orbitals'
write(FILENAME,'(6(a))') 'NTORB.SF.',trim(adjustl(STATENAME)),'.',Spin,'.',NTOType
write(molden_name,'(6(a))') 'MD_NTO.SF.',trim(adjustl(STATENAME)),'.',Spin,'.',NTOType
call WRVEC_(FILENAME,LU,'CO',0,NSYM,NBASF,NBASF,PCMO,vDum,EigValArray,vDum,vDum,vDum,v2Dum,Note,0)
call Molden_interface(0,trim(FILENAME),trim(molden_name))

deallocate(PCMO)

end subroutine NTOSymAnalysis
