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

subroutine MECTL(PROP,OVLP,HAM,ESHFT)

use rassi_aux, only: ipglob
use rassi_global_arrays, only: HDIAG
use Cntrl, only: FnEig, iComp, IfDCPL, IFHAM, IfHDia, IfShft, IPUSED, LuEig, NPROP, NSTATE, PNAME, PNUC, PORIG, PrMER, PRXVR, ToFile
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: PROP(NSTATE,NSTATE,NPROP), OVLP(NSTATE,NSTATE), HAM(NSTATE,NSTATE), ESHFT(NSTATE)
integer(kind=iwp) :: i, iDisk, IEND, IfHD, IFON, iProp, ISTA, iState, j, jState, nAtom, nCol, NST
real(kind=wp) :: PLIMIT, PMAX, X
real(kind=wp), allocatable :: DerCpl(:), NucChg(:)

! Print results:
NCOL = 4
if (PRXVR .and. (IPGLOB >= 0)) then
  write(u6,*)
  call CollapseOutput(1,'Expectation values for input states')
  write(u6,'(3X,A)') '-----------------------------------'
  write(u6,*)
  write(u6,*) ' EXPECTATION VALUES OF 1-ELECTRON OPERATORS'
  write(u6,*) ' FOR THE RASSCF INPUT WAVE FUNCTIONS:'
  write(u6,*)
  do IPROP=1,NPROP
    if (IPUSED(IPROP) /= 0) then

      ! Skip printing if all the diagonal values are very small
      !  (presumed zero for reasons of selection rules)
      PLIMIT = 1.0e-10_wp
      PMAX = Zero
      do I=1,NSTATE
        PMAX = max(PMAX,abs(PROP(I,I,IPROP)+PNUC(IPROP)*OVLP(I,I)))
      end do
      if (PMAX >= PLIMIT) then

        do ISTA=1,NSTATE,NCOL
          IEND = min(NSTATE,ISTA+NCOL-1)
          write(u6,*)
          write(u6,'(1X,A,A8,A,I4)') 'PROPERTY: ',PNAME(IPROP),'   COMPONENT:',ICOMP(IPROP)
          write(u6,'(1X,A,3ES17.8)') 'ORIGIN    :',(PORIG(I,IPROP),I=1,3)
          write(u6,'(1X,A,I8,3I17)') 'STATE     :',(I,I=ISTA,IEND)
          write(u6,*)
          write(u6,'(1X,A,4(1x,G16.9))') 'ELECTRONIC:',(PROP(I,I,IPROP),I=ISTA,IEND)
          write(u6,'(1X,A,4(1x,G16.9))') 'NUCLEAR   :',(PNUC(IPROP)*OVLP(I,I),I=ISTA,IEND)
          write(u6,'(1X,A,4(1x,G16.9))') 'TOTAL     :',(PROP(I,I,IPROP)+PNUC(IPROP)*OVLP(I,I),I=ISTA,IEND)
          write(u6,*)
        end do

      end if
    end if
  end do
  call CollapseOutput(0,'Expectation values for input states')
end if

IFON = 1
X = Zero
do I=2,NSTATE
  do J=1,I-1
    X = max(X,abs(OVLP(I,J)))
  end do
end do
if (X >= 1.0e-6_wp) IFON = 0
IFHD = 1
X = Zero
do I=2,NSTATE
  do J=1,I-1
    X = max(X,abs(HAM(I,J)))
  end do
end do
if (X >= 1.0e-6_wp) IFHD = 0
if (IFON == 0) IFHD = 0

if (IPGLOB >= 2) then
  if (IFHAM) then
    write(u6,*)
    write(u6,*) ' HAMILTONIAN MATRIX FOR THE ORIGINAL STATES:'
    write(u6,*)
    if (IFHD == 1) then
      write(u6,*) ' Diagonal, with energies'
      write(u6,'(5(1X,F15.8))') (HAM(J,J),J=1,NSTATE)
      !do J=1,NSTATE
      !  write(u6,*) HAM(J,J)
      !end  do
    else
      do ISTA=1,NSTATE,5
        IEND = min(ISTA+4,NSTATE)
        write(u6,'(10X,5(8X,A3,I4,A3))') (' | ',I,' > ',I=ISTA,IEND)
        do J=1,NSTATE
          write(u6,'(A3,I4,A3,5(2X,F16.8))') ' < ',J,' | ',(HAM(I,J),I=ISTA,IEND)
        end do
      end do
    end if
  end if
end if

if (IPGLOB >= 2) then
  write(u6,*)
  write(u6,*) '     OVERLAP MATRIX FOR THE ORIGINAL STATES:'
  write(u6,*)
  if (IFON == 1) then
    write(u6,*) ' Diagonal, with elements'
    write(u6,'(5(1X,F15.8))') (OVLP(J,J),J=1,NSTATE)
  else
    do ISTATE=1,NSTATE
      write(u6,'(5(1X,F15.8))') (OVLP(ISTATE,J),J=1,ISTATE)
    end do
  end if
end if

! Addition by A.Ohrn. If ToFile keyword has been specified, we put
! numbers on the auxiliary rassi-to-qmstat file.
if (ToFile) then
  if (.not. IfHam) then
    write(u6,*)
    write(u6,*) 'You ask me to print hamiltonian, but there is none to print!'
    call Abend()
  end if
  call DaName(LuEig,FnEig)
  iDisk = 0
  do iState=1,nState
    do jState=1,iState
      call dDaFile(LuEig,1,Ham(iState,jState),1,iDisk)
    end do
  end do
  do iState=1,nState
    do jState=1,iState
      call dDaFile(LuEig,1,OvLp(iState,jState),1,iDisk)
    end do
  end do
  ! File is closed in eigctl.
  call DaClos(LuEig)
end if
! End of addition by A.Ohrn.

if (IFHAM .and. (IFHDIA .or. IFSHFT)) then
  do ISTATE=1,NSTATE
    if (.not. IFSHFT) ESHFT(ISTATE) = Zero
    if (IFHDIA) ESHFT(ISTATE) = ESHFT(ISTATE)+(HDIAG(ISTATE)-HAM(ISTATE,ISTATE))
  end do
  do ISTATE=1,NSTATE
    do JSTATE=1,NSTATE
      HAM(ISTATE,JSTATE) = HAM(ISTATE,JSTATE)+Half*(ESHFT(ISTATE)+ESHFT(JSTATE))*OVLP(ISTATE,JSTATE)
    end do
  end do
  IFHD = 1
  do I=2,NSTATE
    do J=1,J-1
      if (abs(HAM(I,J)) >= 1.0e-10_wp) IFHD = 0
    end do
  end do
  if (IFON == 0) IFHD = 0
  if (IPGLOB >= 2) then
    write(u6,*)
    write(u6,*) ' USER-MODIFIED HAMILTONIAN FOR THE ORIGINAL STATES:'
    write(u6,*) '(With user shifts, and/or replaced diagonal'
    write(u6,*) ' elements, including overlap corrections.)'
    if (IFHD == 1) then
      write(u6,*) ' Diagonal, with energies'
      write(u6,'(5(1X,F15.8))') (HAM(J,J),J=1,NSTATE)
      !do J=1,NSTATE
      !  write(u6,*) HAM(J,J)
      !end do
    else
      do ISTATE=1,NSTATE
        write(u6,'(5(1X,F15.8))') (HAM(ISTATE,J),J=1,ISTATE)
      end do
    end if
  end if
end if
!PAM00 End of updated HDIA/SHIFT section.

if ((IPGLOB > 0) .and. PRMER) then
  write(u6,*)
  call CollapseOutput(1,'Matrix elements for input states')
  write(u6,'(3X,A)') '--------------------------------'
  write(u6,*)
  write(u6,*) ' MATRIX ELEMENTS OF 1-ELECTRON OPERATORS'
  write(u6,*) ' FOR THE RASSCF INPUT WAVE FUNCTIONS:'
  write(u6,*)
  do IPROP=1,NPROP
    if (IPUSED(IPROP) /= 0) then
      do ISTA=1,NSTATE,NCOL
        IEND = min(NSTATE,ISTA+NCOL-1)
        write(u6,*)
        write(u6,'(1X,A,A8,A,I4)') 'PROPERTY: ',PNAME(IPROP),'   COMPONENT:',ICOMP(IPROP)
        write(u6,'(1X,A,3ES17.8)') 'ORIGIN    :',(PORIG(I,IPROP),I=1,3)
        write(u6,'(1X,A,I8,3I17)') 'STATE     :',(I,I=ISTA,IEND)
        write(u6,*)
        do J=1,NSTATE
          write(u6,'(1X,I5,6X,4(1x,G16.9))') J,(PROP(J,I,IPROP)+PNUC(IPROP)*OVLP(J,I),I=ISTA,IEND)
        end do
        write(u6,*)
      end do
    end if
  end do
  call CollapseOutput(0,'Matrix elements for input states')
end if
!nf
if (IfDCpl) then
  call Get_iScalar('Unique atoms',natom)
  call mma_allocate(NucChg,natom,Label='NucChg')
  call Get_dArray('Effective nuclear Charge',NucChg,nAtom)
  nST = nState*(nState+1)/2
  call mma_allocate(DerCpl,3*natom*nST,Label='DerCpl')
  call AppDerCpl(natom,nST,NucChg,Prop,DerCpl,HAM)
  call mma_deallocate(DerCpl)
  call mma_deallocate(NucChg)
end if
!nf
call Put_dArray('SFS_HAM',HAM,NSTATE**2)
call Put_dArray('SFS_OVLP',OVLP,NSTATE**2)

write(u6,*)

end subroutine MECTL
