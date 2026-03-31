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
! Copyright (C) 2020, Bruno Tenorio                                    *
!***********************************************************************

! Print the reduced 2-e TDM in ASCII format.
! Code adapted from trd_print written by P. A. Malmqvist.
! -------------------------------------------------------------
! The spin coupling matrix elements have the following index-code:
!             SPIN=1 means  K2V (AAB+BBB)
!             SPIN=-1 means SDA (AAA+BBA)
! Notice, SPIN here has nothing to do with the spin quantum number. It
! is just a printing code.
! ------------------------------------------------------------
subroutine RTDM2_PRINT(ISTATE,JSTATE,EIJ,NDYSAB,DYSAB,NRT2MAB,RT2M,CMO1,CMO2,AUGSPIN)

use Cntrl, only: LSYM1, LSYM2, OCAA, OCAN
use Symmetry_Info, only: MUL, nIrrep
use rassi_data, only: NASH, NASHT, NBASF, NFRO, NISH, NOSH
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ISTATE, JSTATE, SYM12, NDYSAB, NRT2MAB, AUGSPIN
real(kind=wp) :: EIJ, DYSAB(*), RT2M(*), CMO1(*), CMO2(*)
integer(kind=iwp) :: I, IA, IO, IOFFA(8), IOFFO(8), IOFFTD, ISYI, ISYJ, ISYL, ISYM, J, JA, JO, KPOS, L, LA, LO, LPOS, LU, NAI, &
                     NAJ, NAL, NB, NII, NIJ, NIL, NO, NOI, NOJ, NOL, NORBSYM
character(len=16) :: FNM
character(len=3) :: NUM1, NUM2
integer(kind=iwp), external :: IsFreeUnit

! IOFFA=NR OF ACTIVE ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
IOFFA(1) = 0
do I=1,nIrrep-1
  IOFFA(I+1) = IOFFA(I)+NASH(I)
end do
! IOFFO=NR OF OCC ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
IOFFO(1) = 0
do J=1,nIrrep-1
  IOFFO(J+1) = IOFFO(J)+NOSH(J)
end do
! Subroutine starts
LU = 51
LU = IsFreeUnit(LU)
write(NUM1,'(I3.3)') ISTATE
write(NUM2,'(I3.3)') JSTATE
! AUGSPIN
if (AUGSPIN == 1) then
  FNM = 'r2TM_K2V_'//NUM1//'_'//NUM2
else if (AUGSPIN == -1) then
  FNM = 'r2TM_SDA_'//NUM1//'_'//NUM2
!else if (aab) then ! if AAB (all spin) true
!  select case (AUGSPIN)
!    case (2)
!      FNM = 'r2TM_BBB_'//NUM1//'_'//NUM2
!    case (3)
!      FNM = 'r2TM_AAA_'//NUM1//'_'//NUM2
!    case (4)
!      FNM = 'r2TM_AAB_'//NUM1//'_'//NUM2
!    case (5)
!      FNM = 'r2TM_BBA_'//NUM1//'_'//NUM2
!    case (6)
!      FNM = 'r2TM_ABA_'//NUM1//'_'//NUM2
!    case (7)
!      FNM = 'r2TM_BAB_'//NUM1//'_'//NUM2
!  end select
end if
call Molcas_Open(LU,FNM)
write(LU,*) '# Auger Densities: CMO1, CMO2, 1-e Dyson, 2-e Dyson'
if (AUGSPIN == 1) then
  write(LU,*) '# Spin Matrix (K-2V) for RAES and NAES'
else if (AUGSPIN == -1) then
  write(LU,*) '# Spin matrix SDA for NAES.'
!else if (aab) then ! if AAB (all spin) true
!  select case (AUGSPIN)
!    case (2)
!      write(LU,*) '# Spin BBB 2-el rTDM density matrix.'
!    case (3)
!      write(LU,*) '# Spin AAA 2-el rTDM density matrix.'
!    case (4)
!      write(LU,*) '# Spin AAB 2-el rTDM density matrix.'
!    case (5)
!      write(LU,*) '# Spin BBA 2-el rTDM density matrix.'
!    case (6)
!      write(LU,*) '# Spin ABA 2-el rTDM density matrix.'
!    case (7)
!      write(LU,*) '# Spin BAB 2-el rTDM density matrix.'
!  end select
end if

write(LU,*) '# OCA for scattering atom:'
write(LU,*) OCAN
do I=1,OCAN
  write(LU,*) OCAA(I)
end do
write(LU,*) '# Binding energy (eV)'
write(LU,'(5ES19.12)') EIJ

SYM12 = MUL(LSYM1,LSYM2)
write(LU,*) '# Total Symmetry of the WF product (<N-1,N>):'
write(LU,*) SYM12
write(LU,*) '# States:',ISTATE,JSTATE
write(LU,*) '# Nr of irreps:',nIrrep
write(LU,*) nIrrep
write(LU,'(A17,8I5)') ' # Basis func   :',(NBASF(ISYM),ISYM=1,nIrrep)
write(LU,'(8I5)') (NBASF(ISYM),ISYM=1,nIrrep)
write(LU,'(A17,8I5)') ' # Frozen  orb   :',(NFRO(ISYM),ISYM=1,nIrrep)
write(LU,'(A17,8I5)') ' # Inactive orb  :',(NISH(ISYM),ISYM=1,nIrrep)
write(LU,'(A17,8I5)') ' # Active orb    :',(NASH(ISYM),ISYM=1,nIrrep)
write(LU,'(8I5)') (NASH(ISYM),ISYM=1,nIrrep)
write(LU,'(A17,8I5)') ' # Total num orb :',(NOSH(ISYM),ISYM=1,nIrrep)
write(LU,'(8I5)') (NOSH(ISYM),ISYM=1,nIrrep)
write(LU,*) '# CMO1 Molecular orbitals'
LPOS = 1
do ISYM=1,nIrrep
  NO = NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
  NB = NBASF(ISYM)
  do IO=1,NO
    write(LU,*) '# Symm ',ISYM,'   Orbital ',IO
    do i=0,NB-1
      write(LU,'(5ES19.12)') CMO1(LPOS+NB*(IO-1)+i)
    end do
  end do
  LPOS = LPOS+NB*NO
end do
write(LU,*) '# CMO2 Molecular orbitals'
LPOS = 1
do ISYM=1,nIrrep
  NO = NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
  NB = NBASF(ISYM)
  do IO=1,NO
    write(LU,*) '# Symm ',ISYM,'   Orbital ',IO
    do i=0,NB-1
      write(LU,'(5ES19.12)') CMO2(LPOS+NB*(IO-1)+i)
    end do
  end do
  LPOS = LPOS+NB*NO
end do
! Write Dyson orbitals in CI basis
write(LU,*) '# 1-e Dyson orbital for CI coeff. in MO biorth. basis'
write(LU,*) '# Symmetry Block elements:',NDYSAB
write(LU,*) '# sub-Block info:Sym(I), NumOrb in SymmBlock'
IOFFTD = 0
do ISYI=1,nIrrep
  NOI = NOSH(ISYI)
  NII = NISH(ISYI)
  if (NOI == 0) goto 400
  write(LU,'(A10,8I7)') ' # sub-Block:',ISYI,NOI
  do I=1,NOI
    IA = I+IOFFTD
    ! eliminate small numbers
    if (abs(DYSAB(IA)) < 1.0e-29_wp) DYSAB(IA) = Zero
    write(LU,'(I7,ES22.12)') IA,DYSAB(IA)
  end do
400 continue
  NORBSYM = NOI
  IOFFTD = IOFFTD+NORBSYM
end do
! Write reduced 2-e TDM in CI basis.
IOFFTD = 0
write(LU,*) '# 2-e reduced TDM for CI coeff. in MO biorth. basis'
write(LU,*) '# Symmetry Block elements:',NRT2MAB
write(LU,*) NRT2MAB
write(LU,*) '# sub-Block info:Sym(I,J,L), NumOrb in SymmBlock'
do ISYI=1,nIrrep
  NOI = NOSH(ISYI)
  NAI = NASH(ISYI)
  NII = NISH(ISYI)
  if (NOI == 0) goto 270
  do ISYJ=1,nIrrep
    NOJ = NOSH(ISYJ)
    NAJ = NASH(ISYJ)
    NIJ = NISH(ISYJ)
    if (NOJ == 0) goto 370
    do ISYL=1,nIrrep
      NOL = NOSH(ISYL)
      NAL = NASH(ISYL)
      NIL = NISH(ISYL)
      if (NOL == 0) goto 470
      if (MUL(ISYI,MUL(ISYJ,ISYL)) == SYM12) then
        if (NAI == 0) goto 670
        if (NAJ == 0) goto 670
        if (NAL == 0) goto 670
        write(LU,'(A10,8I7,8I7,8I7,8I7)') ' # sub-Block:',ISYI,ISYJ,ISYL,NOI*NOJ*NOL
        do I=1,NOI
          IA = IOFFA(ISYI)+I-NII
          IO = IOFFO(ISYI)+I
          do J=1,NOJ
            JA = IOFFA(ISYJ)+J-NIJ
            JO = IOFFO(ISYJ)+J
            do L=1,NOL
              LA = IOFFA(ISYL)+L-NIL
              LO = IOFFO(ISYL)+L
              if ((IA <= 0) .or. (JA <= 0) .or. (LA <= 0)) then
                write(LU,'(I7,I7,I7,ES26.12)') IO,JO,LO,Zero
              else
                KPOS = IA+NASHT*((LA+NASHT*(JA-1))-1)
                if (abs(RT2M(KPOS)) < 1.0e-39_wp) RT2M(KPOS) = Zero
                write(LU,'(I7,I7,I7,ES26.12)') IO,JO,LO,RT2M(KPOS)
              end if
            end do
          end do
        end do
670     continue
      end if
470   continue
    end do
370 continue
  end do
270 continue
end do
close(LU)

end subroutine RTDM2_PRINT
