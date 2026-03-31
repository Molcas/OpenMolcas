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
! Copyright (C) 2011, Per Ake Malmqvist                                *
!***********************************************************************

! Print the transition density matrices in ASCII format.
!   Code written by P. A. Malmqvist.
! This code was moved from the main gtdmctl file for clarity.
! - F. Plasser
subroutine TRD_PRINT(ISTATE,JSTATE,DO22,TDMAB,TDM2,CMO1,CMO2,SIJ)

use Cntrl, only: LSYM1, LSYM2
use Symmetry_Info, only: MUL, nIrrep
use rassi_data, only: NASHT, NAES, NASH, NBASF, NFRO, NISH, NOSH
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ISTATE, JSTATE
logical(kind=iwp) :: DO22
real(kind=wp) :: TDMAB(*), TDM2(*), CMO1(*), CMO2(*), SIJ
integer(kind=iwp) :: I, II, IO, ISYM, ISYM1, ISYM2, ISYT, ISYU, ISYV, ISYX, IT, ITABS, ITU, ITUVX, IU, IUABS, IV, IVABS, IVX, &
                     IWBUF, IX, IXABS, JJ, LIMX, LPOS, LSYM12, LU, NA1, NA2, NB, NI1, NI2, NO, NO1, NO2
real(kind=wp) :: WBUF(5)
character(len=12) :: FNM
character(len=3) :: NUM1, NUM2
integer(kind=iwp), external :: IsFreeUnit

LU = IsFreeUnit(50)
write(NUM1,'(I3.3)') ISTATE
write(NUM2,'(I3.3)') JSTATE
FNM = 'TRD2_'//NUM1//'_'//NUM2
call Molcas_Open(LU,FNM)
write(LU,*) '#Transition density file from RASSI.'
write(LU,*) '#  States:'
write(LU,*) ISTATE,JSTATE
write(LU,*) '#  Nr of irreps:'
write(LU,*) nIrrep
write(LU,*) '#  Basis functions:'
write(LU,'(8I5)') (NBASF(ISYM),ISYM=1,nIrrep)
write(LU,*) '#  Frozen orbitals:'
write(LU,'(8I5)') (NFRO(ISYM),ISYM=1,nIrrep)
write(LU,*) '#  Inactive orbitals:'
write(LU,'(8I5)') (NISH(ISYM),ISYM=1,nIrrep)
write(LU,*) '#  Active orbitals:'
write(LU,'(8I5)') (NASH(ISYM),ISYM=1,nIrrep)
write(LU,*) '#  State ',ISTATE,'    CMO coefficients:'
LPOS = 1
do ISYM=1,nIrrep
  NO = NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
  NB = NBASF(ISYM)
  do IO=1,NO
    write(LU,*) '#  Symm ',ISYM,'   Orbital ',IO
    write(LU,'(5ES19.12)') (CMO1(LPOS+NB*(IO-1)+i),i=0,NB-1)
  end do
  LPOS = LPOS+NB*NO
end do
write(LU,*) '#  State ',JSTATE,'    CMO coefficients:'
LPOS = 1
do ISYM=1,nIrrep
  NO = NFRO(ISYM)+NISH(ISYM)+NASH(ISYM)
  NB = NBASF(ISYM)
  do IO=1,NO
    write(LU,*) '#  Symm ',ISYM,'   Orbital ',IO
    write(LU,'(5ES19.12)') (CMO2(LPOS+NB*(IO-1)+i),i=0,NB-1)
  end do
  LPOS = LPOS+NB*NO
end do
write(LU,*) '#  States ',ISTATE,JSTATE,' Overlap:'
write(LU,'(5ES19.12)') SIJ
write(LU,*) '#  States ',ISTATE,JSTATE,' Active TRD1:'
LSYM12 = MUL(LSYM1,LSYM2)
LPOS = 1
do ISYM1=1,nIrrep
  NO1 = NOSH(ISYM1)
  ISYM2 = MUL(ISYM1,LSYM12)
  NO2 = NOSH(ISYM2)
  if (NO1*NO2 > 0) then
    NA1 = NASH(ISYM1)
    NA2 = NASH(ISYM2)
    if (NA1*NA2 > 0) then
      NI1 = NISH(ISYM1)
      NI2 = NISH(ISYM2)
      write(LU,*) '#  Symmetries ',ISYM1,ISYM2
      write(LU,'(5ES19.12)') ((TDMAB(LPOS-1+II+NO1*(JJ-1)),JJ=NI2+1,NO2),II=NI1+1,NO1)
    end if
    LPOS = LPOS+NO1*NO2
  end if
end do

if (DO22) then
  write(LU,*) '#  States ',ISTATE,JSTATE,' Active TRD2:'
  do ISYT=1,nIrrep
    do ISYU=1,nIrrep
      do ISYV=1,ISYT
        LIMX = ISYV
        if (ISYV == ISYT) LIMX = ISYU
        do ISYX=1,LIMX
          !> Write out one symmetry block (4 indices!) of two-electron
          !> transition density matrix elements.
          !> Write a full 'rectangular' array, even if it could be made
          !> smaller by permutation symmetry.
          write(LU,*) '#  Orbital symm:',ISYT,ISYU,ISYV,ISYX
          IWBUF = 0
          do IT=1,NASH(ISYT)
            ITABS = NAES(ISYT)+IT
            do IU=1,NASH(ISYU)
              IUABS = NAES(ISYU)+IU
              ITU = ITABS+NASHT*(IUABS-1)
              do IV=1,NASH(ISYV)
                IVABS = NAES(ISYV)+IV
                do IX=1,NASH(ISYX)
                  IXABS = NAES(ISYX)+IX
                  IVX = IVABS+NASHT*(IXABS-1)
                  if (ITU >= IVX) then
                    ITUVX = (ITU*(ITU-1))/2+IVX
                  else
                    ITUVX = (IVX*(IVX-1))/2+ITU
                  end if
                  IWBUF = IWBUF+1
                  WBUF(IWBUF) = TDM2(ITUVX)
                  if (IWBUF == 5) then
                    write(LU,'(5ES19.12)') (WBUF(I),I=1,IWBUF)
                    IWBUF = 0
                  end if
                end do
              end do
            end do
          end do
          if (IWBUF > 0) then
            write(LU,'(5ES19.12)') (WBUF(I),I=1,IWBUF)
            IWBUF = 0
          end if
          ! End of writing a symmetry block.
        end do
      end do
    end do
  end do
end if
close(LU)

end subroutine TRD_PRINT
