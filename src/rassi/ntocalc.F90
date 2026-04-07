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
!***************************************************
!                   Calculating NTO
!****************************************************
! Reference: J. Chem. Phys., 2014, 141, 024106
! Notation used in the following of the code follows those in the
! reference mentioned above, especially between equation 53 and
! 54 on page 024106-8.
! NASHT is the number of active orbitals, originated from this
! program.
! Umat is the U matrix, which is the eigenvector  matrix
! calculated by transition density matrix (TDM)
! multiplied by its transpose. Vmat is the V matrix, the eigen-
! vector matrix of a matrix calculated by the transpose
! multiplied by the TDM.
! Ueig is the eigenvalen matrix for the U matrix, Veig is
! that
! for the V matrix.
! ONTO is the hole NTO, calculated by multiplying MO matrix with
! the eigenvector matrix for U matrix. Note that the eigenvector
! matrix is still named as U. Similar condition is for the
! particle matirx.
!
! However, the sets of particle and hole orbitals are switched
! when I examined the results. So I put the data stored in LONTO
! as the particle NTO. (Because JOB1 is for the second JobIph
! file in the input and JOB2 is for the first.)
!
!                                  -------Jie Bao
!           in Depart. of Chemistry, University of Minnesota, USA
!                                      2018/08/09

subroutine NTOCalc(JOB1,JOB2,ISTATE,JSTATE,TRAD,TRASD,ISpin)

use fortran_strings, only: str
use Symmetry_Info, only: nIrrep
use rassi_data, only: NASH, NASHT, NBASF, NBST, NCMO, NISH, NISHT, NOSH
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: JOB1, JOB2, IState, jState, ISpin
real(kind=wp) :: TRAD(NASHT**2), TRASD(NASHT**2)
integer(kind=iwp) :: I, I_NTO, icactorb, INFO, IOrb, Iprint, isym, IUseSym, J, Jprint, LU, N_NTO, NAISHT, NDge, NScrq, NSupBas, &
                     NSUPCMO, NUseBF(nIrrep), NUsedBF(nIrrep), NUseSym, RealtoUse(nIrrep), UsetoReal(nIrrep)
real(kind=wp) :: SumEigVal, WGRONK(2)
logical(kind=iwp) :: DOTEST
character(len=128) :: FILENAME
character(len=9) :: STATENAME
character(len=8) :: NTOType
character(len=3) :: lIrrep(8)
character :: Spin(2)
integer(kind=iwp), allocatable :: Indfr(:), Indto(:), OrbAct(:), OrbBas(:), OrbSym(:), OrbUsedSym(:), Symfr(:), Symto(:)
real(kind=wp), allocatable :: CMO1(:), CMO2(:), ONTO(:), Scrq(:), SUPCMO1(:), SUPCMO2(:), TDM(:), TDMT(:), UEig(:), UMAT(:), &
                              UNTO(:), VEig(:), VMAT(:)
real(kind=wp), parameter :: PrintThres = 1.0e-5_wp
integer(kind=iwp), external :: ISFREEUNIT

! Nr. of basis functions used prior to this symmetry (NUsedBF)
! and used in this symmetry (NUseBF) nIrrep >= NusedSym
!IOrb is the index  of orbitals.
!OrbBas() is the number of basis function for IOrb
!OrbSym() is the index of symmetry/irrep  for IOrb

statename = ''
DoTest = .false.
call mma_allocate(CMO1,NCMO,Label='CMO1')
call mma_allocate(CMO2,NCMO,Label='CMO2')
call RDCMO_RASSI(JOB1,CMO1)
call RDCMO_RASSI(JOB2,CMO2)

if (dotest) then
  write(u6,*) 'CMO1'
  do I=0,NCMO,5
    write(u6,'(2X,5F10.6)') (CMO1(I+IPrint),IPrint=1,min(5,NCMO-I))
  end do
  write(u6,*) 'CMO2'
  do I=0,NCMO,5
    write(u6,'(2X,5F10.6)') (CMO2(I+IPrint),IPrint=1,min(5,NCMO-I))
  end do
end if

call mma_allocate(OrbUsedSym,NBST,Label='OrbUsedSym')
call mma_allocate(OrbAct,NISHT+NASHT,Label='OrbAct')
call mma_allocate(OrbBas,NISHT+NASHT,Label='OrbBas')
call mma_allocate(OrbSym,NISHT+NASHT,Label='OrbSym')

! Analyzing the symmetry of the wave function
IUseSym = 0
NSupBas = 0
IOrb = 0
do ISym=1,nIrrep
  if (NASH(ISym) > 0) then
    IUseSym = IUseSym+1
    RealtoUse(ISym) = IUseSym
    UsetoReal(IUseSym) = ISym
    NSupBas = NSupBas+NBASF(ISym)
    if (IUseSym > 1) then
      NUsedBF(IUseSym) = NBASF(UsetoReal(IUseSym-1))+NUsedBF(IUseSym-1)
    else
      NUsedBF(IUseSym) = 0
    end if
    NUseBF(IUseSym) = NBASF(ISym)
    do I=1,NOSH(ISym)
      IOrb = IOrb+1
      if (I > NISH(ISym)) then
        OrbAct(IOrb) = 1
      else
        OrbAct(IOrb) = 0
      end if
      OrbBas(IOrb) = NBASF(ISym)
      OrbSym(IOrb) = ISym
      OrbUsedSym(IOrb) = IUseSym
    end do
  else
    RealtoUse(ISym) = 0
  end if
end do
NSUPCMO = NASHT*NSupBas
NUseSym = IUseSym
NAISHT = NASHT+NISHT
if (DoTest) then
  write(u6,*) 'Reprinting MO information'
  write(u6,*) 'Size of Super-CMO matrix',NSupCMO
  write(u6,'(6X,A20,4X,16I4)') 'MO Index',(IOrb,IOrb=1,NAISHT)
  write(u6,'(6X,A20,4X,16I4)') 'Irrep Belong to',(OrbSym(IOrb),IOrb=1,NAISHT)
  write(u6,'(6X,A20,4X,16I4)') 'Nr. of Basis F',(OrbBas(IOrb),IOrb=1,NAISHT)
  write(u6,'(6X,A20,4X,16I4)') 'Act Orbital?',(OrbAct(IOrb),IOrb=1,NAISHT)
  write(u6,'(6X,A20,4X,16I4)') 'used basis f',(NUsedBF(OrbUsedSym(IOrb)),IOrb=1,NAISHT)
end if
! End of analyzing wave function

! building up a super-CMO matrix (to be C1-like)
call mma_allocate(SUPCMO1,NSUPCMO,Label='SUPCMO1')
call mma_allocate(SUPCMO2,NSUPCMO,Label='SUPCMO1')
SUPCMO1(:) = Zero
SUPCMO2(:) = Zero
call mma_allocate(ONTO,NSUPCMO,Label='ONTO')
call mma_allocate(UNTO,NSUPCMO,Label='UNTO')
icactorb = 0
I = 0
do IOrb=1,NAISHT
  if (OrbAct(IOrb) == 1) then
    icactorb = icactorb+1
    do IPrint=1,(OrbBas(IOrb))
      JPrint = IPrint+NUsedBF(OrbUsedSym(IOrb))
      J = I+IPrint-1
      SUPCMO1(icactorb+(JPRINT-1)*NASHT) = CMO1(1+J)
      SUPCMO2(icactorb+(JPRINT-1)*NASHT) = CMO2(1+J)
    end do
  end if
  I = I+OrbBas(IOrb)
end do
if (DoTest) then
  write(u6,*) 'printing CMO1 in a C1-like format'
  do I=1,NASHT
    do J=1,NSupBas,10
      write(u6,'(4X,10F10.6)') (SUPCMO1(I+(JPrint-1)*NASHT),JPrint=J,min(J+9,NSupBas))
    end do
  end do

  write(u6,*) 'printing CMO2 in a C1-like format'
  do I=1,NASHT
    do J=1,NSupBas,10
      write(u6,'(4X,10F10.6)') (SUPCMO2(I+(JPrint-1)*NASHT),JPrint=J,min(J+9,NSupBas))
    end do
  end do
end if
call mma_deallocate(OrbUsedSym)
call mma_deallocate(OrbAct)
call mma_deallocate(OrbBas)
call mma_deallocate(OrbSym)
! end of building up the super-CMO matrix
!     Start and initialize spaces
statename = str(JSTATE)//'_'//str(ISTATE)
NDge = NASHT**2
call mma_allocate(Umat,NDge,Label='UMat')
call mma_allocate(Vmat,NDge,Label='VMat')
call mma_allocate(Ueig,NDge,Label='Ueig')
call mma_allocate(Veig,NDge,Label='Veig')
call mma_allocate(TDM,NDge,Label='TDM')
call mma_allocate(TDMT,NDge,Label='TDMT')
write(u6,*)
write(u6,'(6X,A)') repeat('*',100)
write(u6,'(6X,A,98X,A)') '*','*'
write(u6,'(6X,A,34X,A31,33X,A)') '*','NATURAL TRANSITION ORBITALS','*'
write(u6,'(6X,A,98X,A)') '*','*'
write(u6,'(6X,A,38X,A14,I2,A12,I2 ,30X,A )') '*','BETWEEN STATE ',JSTATE,' AND STATE ',ISTATE,'*'
write(u6,'(6X,A,98X,A)') '*','*'
write(u6,'(6X,A)') repeat('*',100)
write(u6,*)
write(u6,*)
if (ISpin == 1) then
  N_NTO = 1
  Spin(1) = 'a'
  write(u6,'(10X,a)') 'NTO CALCULATION ONLY DONE FOR ALPHA SPIN BECAUSE THE WAVE FUNCTION IS A SINGLET,'
  write(u6,'(10X,a)') 'SO ALPHA NTOS ARE EQUAL TO BETA ONES'
else
  N_NTO = 2
  Spin(1) = 'a'
  Spin(2) = 'b'
end if
do I_NTO=1,N_NTO
  Ueig(:) = Zero
  Veig(:) = Zero
  ONTO(:) = Zero
  UNTO(:) = Zero
  if (Spin(I_NTO) == 'a') then
    TDM(1:Ndge) = Half*(TRAD(1:Ndge)+TRASD(1:Ndge))
  else
    TDM(1:Ndge) = Half*(TRAD(1:Ndge)-TRASD(1:Ndge))
  end if
  do I=1,NASHT
    do J=1,NASHT
      TDMT(I+NASHT*(J-1)) = TDM(J+NASHT*(I-1))
    end do
  end do
  ! Print out transition density matrix
  if (Spin(I_NTO) == 'a') then
    write(FILENAME,'(a,a)') 'TDM',trim(adjustl(STATENAME))
    LU = ISFREEUNIT(233)
    call Molcas_Open(LU,FILENAME)
    do I=1,NASHT
      write(LU,'(5(1X,ES11.4E2))') (TRAD(NASHT*(I-1)+J),J=1,NASHT)
    end do
    write(LU,*)
    write(LU,*)
    write(LU,*)
    do I=1,NASHT
      write(LU,'(5(1X,ES11.4E2))') (TRASD(NASHT*(I-1)+J),J=1,NASHT)
    end do
    close(LU)
  end if
  ! Generalizing transpose of TDM, TDM_T

  ! Calculating T_trans*T
  call DGEMM_('n','n',NASHT,NASHT,NASHT,One,TDMT,NASHT,TDM,NASHT,Zero,Vmat,NASHT)
  ! Writing Particle Matrix
  write(FILENAME,'(a,a,a)') 'Dhole.',trim(adjustl(STATENAME)),Spin(I_NTO)
  LU = ISFREEUNIT(LU)
  call Molcas_Open(LU,FILENAME)
  do I=1,NASHT
    write(LU,'(10(1X,ES11.4E2))') (Vmat(NASHT*(I-1)+J),J=1,NASHT)
  end do
  close(LU)
  ! Calculating T*T_transpose
  call DGEMM_('n','n',NASHT,NASHT,NASHT,One,TDM,NASHT,TDMT,NASHT,Zero,Umat,NASHT)
  write(FILENAME,'(a,a,a)') 'Dpart.',trim(adjustl(STATENAME)),Spin(I_NTO)
  LU = ISFREEUNIT(LU)
  call Molcas_Open(LU,FILENAME)
  do I=1,NASHT
    write(LU,'(10(1X,ES11.4E2))') (Umat(NASHT*(I-1)+J),J=1,NASHT)
  end do
  close(LU)
  call DSYEV_('V','U',NASHT,Vmat,NASHT,Veig,WGRONK,-1,INFO)
  NScrq = int(WGRONK(1))
  if (Nscrq == 0) then
    Nscrq = max(NDge,100)
    if (DoTest) write(u6,*) 'Size of scratch space is increased to max(NDge,100)'
  end if
  call mma_allocate(Scrq,NScrq,Label='Scrq')
  !   Diagonalizing matrices
  call DSYEV_('V','U',NASHT,Umat,NASHT,Ueig,Scrq,NScrq,INFO)
  call DSYEV_('V','U',NASHT,Vmat,NASHT,Veig,Scrq,NScrq,INFO)
  call mma_deallocate(Scrq)
  ! Printing some matrices
  if (DoTest) then
    write(FILENAME,'(a,a,a)') 'EigVecHole.',trim(adjustl(STATENAME)),Spin(I_NTO)
    LU = ISFREEUNIT(LU)
    call Molcas_Open(LU,FILENAME)
    do I=1,NASHT
      write(LU,'(5(1X,ES11.4E2))') (Vmat(NASHT*(I-1)+J),J=1,NASHT)
    end do
    close(LU)
    write(FILENAME,'(a,a,a)') 'EigVecPart.',trim(adjustl(STATENAME)),Spin(I_NTO)
    LU = ISFREEUNIT(LU)
    call Molcas_Open(LU,FILENAME)
    do I=1,NASHT
      write(LU,'(5(1X,ES11.4E2))') (Umat(NASHT*(I-1)+J),J=1,NASHT)
    end do
    close(LU)
  end if
  ! End of Diagonalizing the matrices

  ! Constructing hole and particle orbitals

  call DGEMM_('t','n',NASHT,NSupBas,NASHT,One,Umat,NASHT,SupCMO1,NASHT,Zero,ONTO,NASHT)
  call DGEMM_('t','n',NASHT,NSupBas,NASHT,One,Vmat,NASHT,SupCMO2,NASHT,Zero,UNTO,NASHT)

  if (DoTest) then
    write(u6,*) 'printing Particle NTO in a C1-like format'
    do I=1,NASHT
      do J=1,NSupBas,10
        write(u6,'(4X,10F10.6)') (ONTO(I+(JPrint-1)*NASHT),JPrint=J,min(J+9,NSupBas))
      end do
    end do

    write(u6,*) 'printing Hole     NTO in a C1-like format'
    do I=1,NASHT
      do J=1,NSupBas,10
        write(u6,'(4X,10F10.6)') (UNTO(I+(JPrint-1)*NASHT),JPrint=J,min(J+9,NSupBas))
      end do
    end do
  end if

  ! Printing NTOs
  call mma_allocate(Symto,NASHT,Label='Symto')
  call mma_allocate(Indto,NASHT,Label='Indto')
  call mma_allocate(Symfr,NASHT,Label='Symfr')
  call mma_allocate(Indfr,NASHT,Label='Indfr')
  NTOType = 'PART'
  call NTOSymAnalysis(NUseSym,NUseBF,NUsedBF,ONTO,NTOType,STATENAME,Ueig,UsetoReal,RealtoUse,Spin(I_NTO),Symto,Indto,SumEigVal)
  NTOType = 'HOLE'
  call NTOSymAnalysis(NUseSym,NUseBF,NUsedBF,UNTO,NTOType,STATENAME,Veig,UsetoReal,RealtoUse,Spin(I_NTO),Symfr,Indfr,SumEigVal)
  ! End of Printing NTOs

  call Get_cArray('Irreps',lIrrep,24)
  lIrrep(1:nIrrep) = adjustr(lIrrep(1:nIrrep))

  ! Putting particle-hole pairs in the output
  write(u6,*)
  if (I_NTO == 1) then
    write(u6,'(10X,a)') 'NATURAL TRANSITION ORBTIAL INFORMATION FOR ALPHA SPIN'
  else
    write(u6,'(10X,a)') 'NATURAL TRANSITION ORBTIAL INFORMATION FOR BETA  SPIN'
  end if
  write(u6,'(6X,A)') repeat('=',100)
  write(u6,'(10X,5A18)') 'EXCITATION','EIGENVALUE','EXCITATION','HOLE NTO','PARTICLE NTO'
  write(u6,'(10X,A18,18X,3A18)') 'AMPLITUDE','CONTRIBUTION(%)','SYMMETRY INDEX','SYMMETRY INDEX'
  write(u6,'(6X,A)') repeat('-',100)
  do IOrb=NASHT,1,-1
    if (Ueig(IOrb) < PrintThres) exit
    write(u6,'(10X,2(10X,F8.5),10X,F8.2,2(A9,I9))') sqrt(Ueig(IOrb)),Ueig(IOrb),Ueig(IOrb)/SumEigVal*100.0_wp,lIrrep(Symfr(IOrb)), &
                                                    Indfr(IOrb),lIrrep(Symto(IOrb)),Indto(IOrb)
  end do

  write(u6,'(6X,A)') repeat('-',100)
  write(u6,'(6X,A,F8.5)') 'SUM OF EIGENVALUES',SumEigVal
  write(u6,'(6X,A)') repeat('=',100)

  call mma_deallocate(Symto)
  call mma_deallocate(Indto)
  call mma_deallocate(Symfr)
  call mma_deallocate(Indfr)
end do
! End of loop over N_NTO (I_NTO=1 for alpha and 2 for beta)

write(u6,*)
write(u6,'(6X,A)') repeat('*',100)
write(u6,'(6X,A,33X,A34,31X,A)') '*','END OF NATURAL TRANSITION ORBITALS','*'
write(u6,'(6X,A,33X,A14,I2,A12,I2 ,35X,A )') '*','BETWEEN STATE ',JSTATE,' AND STATE ',ISTATE,'*'
write(u6,'(6X,A)') repeat('*',100)

call mma_deallocate(VMat)
call mma_deallocate(UMat)
call mma_deallocate(VEig)
call mma_deallocate(UEig)
call mma_deallocate(TDMT)
call mma_deallocate(TDM)
call mma_deallocate(CMO1)
call mma_deallocate(CMO2)

call mma_deallocate(SUPCMO2)
call mma_deallocate(SUPCMO1)
call mma_deallocate(UNTO)
call mma_deallocate(ONTO)

end subroutine NTOCalc
