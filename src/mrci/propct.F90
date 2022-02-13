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

subroutine PROPCT()

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: I, IDDMO, IDISK, IDUMMY(7,8), IEND, IOPT, IPC, IPROP, IRTC, ISTA, ISTATE, ISYMLB, J, JSTATE, LAFOLD, LCMO, &
                     LCNO, LDAO, LOCC, LPINT, LSCR, LSFOLD, LTDAO
real(kind=wp) :: DUMMY(1)
character(len=100) :: REALNAME
character(len=30) :: REMARK
character(len=8) :: FNAME, LABEL
real(kind=wp), allocatable :: PROP(:,:,:)
#include "mrci.fh"
#include "WrkSpc.fh"

!LCMO = LPRP
!LCNO = LCMO+NCMO
!LOCC = LCNO+NCMO
!LTDAO = LOCC+NBAST
!LDAO = LTDAO
!LAFOLD = LDAO+NBAST**2
!LSFOLD = LAFOLD+NBTRI
!LPINT = LSFOLD+NBTRI
!LSCR = LPINT+NBTRI+4
!NSCR = max(NBTRI,NBMAX**2)
!LTOP = LSCR+NSCR-1
call GETMEM('CMO','ALLO','REAL',LCMO,NCMO)
call GETMEM('CNO','ALLO','REAL',LCNO,NCMO)
call GETMEM('OCC','ALLO','REAL',LOCC,NBAST)
call GETMEM('DAO','ALLO','REAL',LDAO,NBAST**2)
LTDAO = LDAO
call GETMEM('AFOLD','ALLO','REAL',LAFOLD,NBTRI)
call GETMEM('SFOLD','ALLO','REAL',LSFOLD,NBTRI)
call GETMEM('PINT','ALLO','REAL',LPINT,NBTRI+4)
NSCR = max(NBTRI,NBMAX**2)
call GETMEM('SCR','ALLO','REAL',LSCR,NSCR)
! LOOP OVER OPERATORS:
IOPT = 8
NPROP = 0
do I=1,100
  ! PICK UP OPERATOR LABELS FROM ONE-ELECTRON FILE:
  LABEL = 'UNDEF'
  call iRDONE(IRTC,1+IOPT,LABEL,IPC,IDUMMY,ISYMLB)
  if (IRTC /= 0) exit
  IOPT = 16
  !PAM96 if (iand(1,ISYMLB) == 0) cycle
  if (mod(ISYMLB,2) == 0) cycle
  NPROP = NPROP+1
  PNAME(NPROP) = LABEL
  IPCOMP(NPROP) = IPC
  PTYPE(NPROP) = 'HERM'
  if (LABEL == 'VELOCITY') PTYPE(NPROP) = 'ANTI'
  if (LABEL == 'ANGMOM  ') PTYPE(NPROP) = 'ANTI'
end do
call MMA_ALLOCATE(PROP,NRROOT,NRROOT,NPROP,'PROP')
call DCOPY_(NRROOT*NRROOT*NPROP,[Zero],0,PROP,1)
IDDMO = 0
do ISTATE=1,NRROOT
  ! PICK UP DMO
  !PAM04 call dDAFILE(LUEIG,2,HWork(LDMO),NBTRI,IDDMO)
  call dDAFILE(LUEIG,2,Work(LDMO),NBTRI,IDDMO)
  ! PICK UP CMO
  IDISK = ITOC17(1)
  call dDAFILE(LUONE,2,Work(LCMO),NCMO,IDISK)
  ! COMPUTE & WRITE NATURAL ORBITALS
  !PAM04 call NATORB_MRCI(HWork(LCMO),HWork(LDMO),HWork(LCNO),
  call NATORB_MRCI(Work(LCMO),Work(LDMO),Work(LCNO),Work(LOCC),Work(LSCR))
  write(FNAME,'(A5,I2.2)') 'CIORB',ISTATE
  REALNAME = FNAME
  !PAM04 if (ISTATE == 1) call Add_Info('CI_DENS1',HWork(LDMO),1,5)
  if (ISTATE == 1) call Add_Info('CI_DENS1',Work(LDMO),1,5)
  REMARK = '* MRCI  '
  !** Gusarov , include 1st root acpf energy to CiOrb file:
  !if (ICPF == 1) REMARK = '* ACPF  '
  if (ICPF == 1) write(REMARK,'("* ACPF  ",f22.16)') ESMALL(1)+ESHIFT

  call WRVEC(REALNAME,LUVEC,'CO',NSYM,NBAS,NBAS,Work(LCNO),Work(LOCC),Dummy,iDummy,REMARK)
  write(u6,*)
  write(u6,'(A,I2)') ' NATURAL ORBITALS OF STATE NR. ',ISTATE
  write(u6,*) ' FULL SET OF ORBITALS ARE SAVED ON FILE ',REALNAME
  call PRORB(Work(LCNO),Work(LOCC))
  write(u6,*) ' ',('*',I=1,70)
  ! CREATE DAO
  call MKDAO(Work(LCNO),Work(LOCC),Work(LDAO))
  ! CALL PMATEL TO CALCULATE CHARGES AND PROPERTIES.
  ! PUT PROPERTIES INTO APPROPRIATE MATRICES.
  call PMATEL(ISTATE,ISTATE,PROP,Work(LPINT),Work(LSCR),Work(LCNO),Work(LOCC),Work(LSFOLD),Work(LAFOLD),Work(LDAO))
end do
! ENERGIES SAVED FROM PREVIOUS OUTPUT REPEATED HERE FOR CONVENIENCE:
write(u6,*)
write(u6,*) ' SUMMARY OF ENERGIES:'
do ISTA=1,NRROOT,4
  IEND = min(ISTA+3,NRROOT)
  write(u6,'(1X,A,I8,3I16)') '               ROOT:',(I,I=ISTA,IEND)
  write(u6,'(1X,A,4F16.8)') '       TOTAL ENERGY:',(ENGY(I,1),I=ISTA,IEND)
  if (ICPF == 0) then
    write(u6,'(1X,A,4F16.8)') 'DAVIDSON CORRECTION:',(ENGY(I,2),I=ISTA,IEND)
    write(u6,'(1X,A,4F16.8)') '    ACPF CORRECTION:',(ENGY(I,3),I=ISTA,IEND)
  end if
  write(u6,*)
end do
! ---------------------------------------------------
!PAM Grep-able energy output for convenience:
write(u6,*)
write(u6,*) ' Energies, machine-readable format:'
if (ICPF == 0) then
  do I=1,NRROOT
    write(u6,'(1X,A,I3,3(5X,A,F16.8))') ' CI State ',I,'Total energy:',ENGY(I,1),'QDav:',ENGY(I,2),'QACPF:',ENGY(I,3)
  end do
else
  do I=1,NRROOT
    write(u6,'(1X,A,I3,5X,A,F16.8)') ' ACPF State ',I,'Total energy:',ENGY(I,1)
  end do
end if
write(u6,*)
! ---------------------------------------------------
if (NPROP > 0) then
  ! WRITE EXPECTATION VALUES:
  write(u6,*)
  write(u6,*) ' EXPECTATION VALUES OF VARIOUS OPERATORS:'
  write(u6,*) '(Note: Electronic multipoles include a negative sign.)'
  do IPROP=1,NPROP
    if (PTYPE(IPROP) == 'ANTI') cycle
    do ISTA=1,NRROOT,4
      IEND = min(ISTA+3,NRROOT)
      write(u6,*)
      write(u6,'(1X,A,A8,A,I4)') '   PROPERTY :',PNAME(IPROP),'   COMPONENT:',IPCOMP(IPROP)
      write(u6,'(1X,A,3F16.8)') '    GAUGE ORIGIN:',(PORIG(I,IPROP),I=1,3)
      write(u6,'(1X,A,I8,3I16)') '            ROOT:',(I,I=ISTA,IEND)
      write(u6,'(1X,A,4F16.8)') '      ELECTRONIC:',(PROP(I,I,IPROP),I=ISTA,IEND)
      write(u6,'(1X,A,4F16.8)') '         NUCLEAR:',(PNUC(IPROP),I=ISTA,IEND)
      write(u6,'(1X,A,4F16.8)') '           TOTAL:',(PNUC(IPROP)+PROP(I,I,IPROP),I=ISTA,IEND)
    end do
  end do
  write(u6,*)
end if
if (ITRANS /= 0) then
  do ISTATE=2,NRROOT
    do JSTATE=1,ISTATE-1
      ! PICK UP TDMA
      !PAM04 call dDAFILE(LUEIG,2,HWork(LTDMO),NBAST**2,IDDMO)
      call dDAFILE(LUEIG,2,Work(LTDMO),NBAST**2,IDDMO)
      ! CREATE TDAO
      !PAM04 call MKTDAO(HWork(LCMO),HWork(LTDMO),HWork(LTDAO),HWork(LSCR))
      call MKTDAO(Work(LCMO),Work(LTDMO),Work(LTDAO),Work(LSCR))
      ! CALL PMATEL TO CALCULATE TRANSITION PROPERTIES
      ! PUT PROPERTIES INTO APPROPRIATE MATRICES.
      if (NPROP /= 0) call PMATEL(ISTATE,JSTATE,PROP,Work(LPINT),Work(LSCR),Work(LCNO),Work(LOCC),Work(LSFOLD),Work(LAFOLD), &
                                  Work(LTDAO))
    end do
  end do
  if (NPROP /= 0) then
    ! WRITE PROPERTY MATRICES.
    write(u6,*)
    write(u6,*) ' MATRIX ELEMENTS OF VARIOUS OPERATORS:'
    write(u6,*) ' (INCLUDING ANY NUCLEAR CONTRIBUTIONS)'
    do IPROP=1,NPROP
      do I=1,NRROOT
        PROP(I,I,IPROP) = PROP(I,I,IPROP)+PNUC(IPROP)
      end do
    end do
    do IPROP=1,NPROP
      do ISTA=1,NRROOT,4
        IEND = min(ISTA+3,NRROOT)
        write(u6,*)
        write(u6,'(1X,A,A8,A,I4)') '   PROPERTY :',PNAME(IPROP),'   COMPONENT:',IPCOMP(IPROP)
        write(u6,'(1X,A,3F16.8)') '    GAUGE ORIGIN:',(PORIG(I,IPROP),I=1,3)
        write(u6,'(1X,A,I8,3I16)') '            ROOT:',(I,I=ISTA,IEND)
        do J=1,NRROOT
          write(u6,'(15X,I2,4F16.8)') J,(PROP(J,I,IPROP),I=ISTA,IEND)
        end do
      end do
    end do
    write(u6,*)
  end if
end if
call GETMEM('CMO','FREE','REAL',LCMO,NCMO)
call GETMEM('CNO','FREE','REAL',LCNO,NCMO)
call GETMEM('OCC','FREE','REAL',LOCC,NBAST)
call GETMEM('DAO','FREE','REAL',LDAO,NBAST**2)
call GETMEM('AFOLD','FREE','REAL',LAFOLD,NBTRI)
call GETMEM('SFOLD','FREE','REAL',LSFOLD,NBTRI)
call GETMEM('PINT','FREE','REAL',LPINT,NBTRI+4)
call GETMEM('SCR','FREE','REAL',LSCR,NSCR)
call MMA_DEALLOCATE(PROP)

return

end subroutine PROPCT
