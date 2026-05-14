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

subroutine PRCNFTAB(CnfTab,MXPRT)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: CnfTab(*), MxPrt
integer(kind=iwp) :: ICNF, IERR, IFORM, IGAS, ISUM, ISYM, ITYPE, KCNFSTA, KSTA, LENCNF, LENGTH, LGASLIM, LGASORB, LIM1, LIM2, &
                     LINFO, LSYM, MAXOP, MINOP, NCLS, NCNF, NEL, NGAS, NOPN, NORB, NPRT, NSYM, NTAB
character(len=144) :: TEXT

! Sanity test:
NPRT = min(MXPRT,10000)
! header:
NTAB = CnfTab(1)
ITYPE = CnfTab(2)
NEL = CnfTab(3)
NORB = CnfTab(4)
MINOP = CnfTab(5)
MAXOP = CnfTab(6)
NSYM = CnfTab(7)
LSYM = CnfTab(8)
NGAS = CnfTab(9)
IFORM = CnfTab(10)
write(u6,*) '---------------------------------------------------'
write(u6,*) '       Configuration Table Printout'
write(u6,*) ' Table header:'
write(u6,'(1x,a,i16)') ' Table size             NTAB=',NTAB
write(u6,'(1x,a,i16)') ' Nr of electrons         NEL=',NEL
write(u6,'(1x,a,i16)') ' Nr of orbitals         NORB=',NORB
write(u6,'(1x,a,i16)') ' Min nr of open shells MINOP=',MINOP
write(u6,'(1x,a,i16)') ' Max nr of open shells MAXOP=',MAXOP
write(u6,'(1x,a,i16)') ' Point group order      NSYM=',NSYM
write(u6,'(1x,a,i16)') ' Selected symmetry      LSYM=',LSYM
write(u6,'(1x,a,i16)') ' Nr of GAS restrictions NGAS=',NGAS
write(u6,'(1x,a,i16)') ' Configuration format  IFORM=',IFORM
if (ITYPE /= 37) then
  write(u6,*) ' PRCNFTAB error: This is not a configuration table!'
  call ABEND()
end if
IERR = 0
if (NEL < 0) IERR = IERR+1
if (NEL > 99) IERR = IERR+1
if (NORB < 0) IERR = IERR+1
if (NORB > 199) IERR = IERR+1
if (MINOP < 0) IERR = IERR+1
if (MAXOP < MINOP) IERR = IERR+1
if (NSYM < 1) IERR = IERR+1
if (NSYM > 8) IERR = IERR+1
if (LSYM < 0) IERR = IERR+1
if (LSYM > NSYM) IERR = IERR+1
if (NGAS < 0) IERR = IERR+1
if (NGAS > 19) IERR = IERR+1
if (IFORM < 1) IERR = IERR+1
if (IFORM > 4) IERR = IERR+1
if (IERR > 0) then
  write(u6,*) ' PRCNFTAB error: Those values are unacceptable!'
  call ABEND()
end if
LGASORB = 11
LGASLIM = LGASORB+(NSYM+1)*(NGAS+1)
if (NGAS > 0) then
  write(u6,*)
  write(u6,*) ' Orbitals by symmetry in GAS partitions:'
  write(u6,'(A,8I5)') ' Space       Tot     ',(ISYM,ISYM=1,NSYM)
  write(u6,'(1X,A,5X,I5,5X,8I5)') 'Active',CnfTab(LGASORB),(CnfTab(LGASORB+ISYM),ISYM=1,NSYM)
  do IGAS=1,NGAS
    write(u6,'(1X,A,I2,5X,I5,5X,8I5)') 'GAS',IGAS,CnfTab(LGASORB+(NSYM+1)*IGAS),(CnfTab(LGASORB+ISYM+(NSYM+1)*IGAS),ISYM=1,NSYM)
  end do
  write(u6,*) ' GAS orbital partitions, min and max population:'
  do IGAS=1,NGAS
    ISUM = CnfTab(LGASORB+(NSYM+1)*IGAS)
    LIM1 = CnfTab(LGASLIM+2*(IGAS-1))
    LIM2 = CnfTab(LGASLIM+1+2*(IGAS-1))
    write(u6,'(1X,A,I2,5X,3I5)') 'GAS',IGAS,ISUM,LIM1,LIM2
  end do
end if
LINFO = LGASLIM+2*NGAS
write(u6,*)
write(u6,*) ' INFO table starts at LINFO=',LINFO
do NOPN=MINOP,MAXOP
  NCLS = (NEL-NOPN)/2
  do ISYM=1,NSYM
    NCNF = CnfTab(LINFO+0+3*(ISYM-1+NSYM*(NOPN-MINOP)))
    KCNFSTA = CnfTab(LINFO+1+3*(ISYM-1+NSYM*(NOPN-MINOP)))
    LENCNF = CnfTab(LINFO+2+3*(ISYM-1+NSYM*(NOPN-MINOP)))
    if (NCNF /= 0) then
      write(u6,*)
      write(u6,*) '  NOPN ISYM       Nr of conf Start point  Words/config'
      write(u6,'(1X,2I4,5X,3I12)') NOPN,ISYM,NCNF,KCNFSTA,LENCNF
      if (NCNF > NPRT) then
        write(u6,'(1X,A,I5)') ' The first NPRT configurations. NPRT=',NPRT
      else
        write(u6,*) ' Configurations:'
      end if
      KSTA = KCNFSTA
      do ICNF=1,min(NCNF,NPRT)
        call CNF2TXT(IFORM,NORB,NCLS,NOPN,CnfTab(KSTA),LENGTH,TEXT)
        KSTA = KSTA+LENCNF
        if (LENGTH <= 72) then
          write(u6,'(8X,A)') TEXT(1:LENGTH)
        else
          write(u6,'(8X,A)') TEXT(1:72)
          write(u6,'(8X,A)') TEXT(73:LENGTH)
        end if
      end do
      if (NCNF > NPRT) write(u6,*) ' ( ...and more. This list was truncated.)'
    end if
  end do
end do

end subroutine PRCNFTAB
