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

use mrci_global, only: DMO, ENGY, ESHIFT, ESMALL, ICPF, IPCOMP, ITOC17, ITRANS, LUEIG, LUONE, LUVEC, NBAS, NBAST, NBMAX, NBTRI, &
                       NCMO, NRROOT, NPROP, NSYM, PNAME, PNUC, PORIG, PTYPE, TDMO
use OneDat, only: sOpSiz, sRdFst
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: I, IDDMO, IDISK, IDUMMY(7,8), IEND, IOPT, IPC, IPROP, IRTC, ISTA, ISTATE, ISYMLB, J, JSTATE, NSCR
real(kind=wp) :: DUMMY(1)
character(len=100) :: REALNAME
character(len=30) :: REMARK
character(len=8) :: FNAME, LABEL
real(kind=wp), allocatable :: AFOLD(:), CMO(:), CNO(:), DAO(:,:), OCC(:), PINT(:), PROP(:,:,:), SCR(:), SFOLD(:)

call mma_allocate(CMO,NCMO,label='CMO')
call mma_allocate(CNO,NCMO,label='CNO')
call mma_allocate(OCC,NBAST,label='OCC')
call mma_allocate(DAO,NBAST,NBAST,label='DAO')
call mma_allocate(AFOLD,NBTRI,label='AFOLD')
call mma_allocate(SFOLD,NBTRI,label='SFOLD')
call mma_allocate(PINT,NBTRI+4,label='PINT')
NSCR = max(NBTRI,NBMAX**2)
call mma_allocate(SCR,NSCR,label='SCR')
! LOOP OVER OPERATORS:
IOPT = ibset(0,sRdFst)
NPROP = 0
do I=1,100
  ! PICK UP OPERATOR LABELS FROM ONE-ELECTRON FILE:
  LABEL = 'UNDEF'
  call iRDONE(IRTC,ibset(IOPT,sOpSiz),LABEL,IPC,IDUMMY,ISYMLB)
  if (IRTC /= 0) exit
  IOPT = 16
  if (mod(ISYMLB,2) == 0) cycle
  NPROP = NPROP+1
  PNAME(NPROP) = LABEL
  IPCOMP(NPROP) = IPC
  PTYPE(NPROP) = 'HERM'
  if (LABEL == 'VELOCITY') PTYPE(NPROP) = 'ANTI'
  if (LABEL == 'ANGMOM  ') PTYPE(NPROP) = 'ANTI'
end do
call mma_allocate(PROP,NRROOT,NRROOT,NPROP,'PROP')
PROP(:,:,:) = Zero
IDDMO = 0
do ISTATE=1,NRROOT
  ! PICK UP DMO
  call dDAFILE(LUEIG,2,DMO,NBTRI,IDDMO)
  ! PICK UP CMO
  IDISK = ITOC17(1)
  call dDAFILE(LUONE,2,CMO,NCMO,IDISK)
  ! COMPUTE & WRITE NATURAL ORBITALS
  call NATORB_MRCI(CMO,DMO,CNO,OCC,SCR)
  write(FNAME,'(A5,I2.2)') 'CIORB',ISTATE
  REALNAME = FNAME
  if (ISTATE == 1) call Add_Info('CI_DENS1',DMO,1,5)
  REMARK = '* MRCI  '
  !** Gusarov , include 1st root acpf energy to CiOrb file:
  if (ICPF == 1) write(REMARK,'("* ACPF  ",f22.16)') ESMALL(1)+ESHIFT

  call WRVEC(REALNAME,LUVEC,'CO',NSYM,NBAS,NBAS,CNO,OCC,Dummy,iDummy,REMARK)
  write(u6,*)
  write(u6,'(A,I2)') ' NATURAL ORBITALS OF STATE NR. ',ISTATE
  write(u6,*) ' FULL SET OF ORBITALS ARE SAVED ON FILE ',REALNAME
  call PRORB(CNO,OCC)
  write(u6,*) ' ',repeat('*',70)
  ! CREATE DAO
  call MKDAO(CNO,OCC,DAO)
  ! CALL PMATEL TO CALCULATE CHARGES AND PROPERTIES.
  ! PUT PROPERTIES INTO APPROPRIATE MATRICES.
  call PMATEL(ISTATE,ISTATE,PROP,PINT,SCR,CNO,OCC,SFOLD,AFOLD,DAO)
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
      call dDAFILE(LUEIG,2,TDMO,NBAST**2,IDDMO)
      ! CREATE TDAO
      call MKTDAO(CMO,TDMO,DAO,SCR)
      ! CALL PMATEL TO CALCULATE TRANSITION PROPERTIES
      ! PUT PROPERTIES INTO APPROPRIATE MATRICES.
      if (NPROP /= 0) call PMATEL(ISTATE,JSTATE,PROP,PINT,SCR,CNO,OCC,SFOLD,AFOLD,DAO)
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
call mma_deallocate(CMO)
call mma_deallocate(CNO)
call mma_deallocate(OCC)
call mma_deallocate(DAO)
call mma_deallocate(AFOLD)
call mma_deallocate(SFOLD)
call mma_deallocate(PINT)
call mma_deallocate(SCR)
call MMA_DEALLOCATE(PROP)

return

end subroutine PROPCT
