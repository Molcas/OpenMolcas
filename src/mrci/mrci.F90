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

subroutine MRCI(IRETURN)
!***********************************************************************
!  MULTI REFERENCE SDCI AND AVERAGE CPF PROGRAM.                       *
!***********************************************************************
! UNITS USED IN THE PROGRAM
! UNIT u5, INPUT
! UNIT u6, OUTPUT
! UNIT  2=LUPROP, (DA,ONEINT) FOR PROPERTY CALCULATIONS
! UNIT 10=LUSYMB, (DA,CIGUGA) SYMBOLIC FORMULAS
! UNIT 50=LUTRA, (DA,TRAINT) TRANSFORMED MO 2-EL INTEGRALS
! UNIT 60, (DA) SORTED AIBJ, ABIJ AND AIJK INTEGRALS
! UNIT 70, (DA) SORTED IJKL AND ABCI INTEGRALS
! UNIT 80, (DA) SORTED ABCD INTEGRALS
! UNIT 17=LUONE, (DA,TRAONE) ONE ELECTRON INTEGRALS
! UNIT 18=LUVEC, (Formatted, sequential!) MRCI ORBITALS OUT
! UNIT 23=LUEIG, (DA) WORKSPACE FOR MALMQUIST DIAGONALIZATION.
! UNIT 25, (DA) FOCK MATRIX AND DIAGONAL CSF MATRIX ELEMENTS
! UNIT 27, (DA) SCRATCH IN IIJJ. ALSO, REFERENCE CI VECTOR.
! UNIT 28=LUREST, (DA,MRCIVECT) CI VECTOR
!***********************************************************************

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: IRETURN
#include "mrci.fh"

! Prologue, print program header
!
!call SETTIM()
!call XUFLOW()
!
! ( Workspace allocated in Start() )
!
!call IniMem()
!PAM04 call GetMem('WrkSpc','Max ','Real',MemOff,MaxMem)
!PAM04 ! PAM July 2004:
!PAM04 ! Actually allocate only half the memory, minus some spare overhead.
!PAM04 ! This is a temporary measure. I will systematically change the
!PAM04 ! present static allocation in order to use GETMEM instead.
!PAM04 ! Changed statement:
!PAM04 !MaxMem = MaxMem-3*1000
!PAM04 MaxMem = (MaxMem-3*1000)/2
!PAM04 call GetMem('WrkSpc','Allo','Real',MemOff,MaxMem)
!PAM04 write(u6,*) " Allocated 'WrkSpc'. memoff, maxmem=",memoff,maxmem

!PAM04 Now try completely without "WrkSpc" array:
call mma_maxdble(MemTot)

! Open files

LUVEC = 18

LUSYMB = 10
call DANAME(LUSYMB,'CIGUGA')
LUTRA = 50
call DANAME_MF(LUTRA,'TRAINT')
LUONE = 17
call DANAME(LUONE,'TRAONE')
LUREST = 28
call DANAME(LUREST,'MRCIVECT')
! Temporaries:
Lu_60 = 60
call DANAME_MF(Lu_60,'TIABIJ')
Lu_70 = 70
call DANAME_MF(Lu_70,'TIABCI')
Lu_80 = 80
call DANAME_MF(Lu_80,'TIABCD')
LUEIG = 23
call DANAME(LUEIG,'FT23F001')
Lu_25 = 25
call DANAME(Lu_25,'FT25F001')
Lu_27 = 27
call DANAME(Lu_27,'FT27F001')

! main body

!PAM04 iMemOff = ip_of_iWork(Work(MemOff))
!PAM04 call SDCI(Work(MemOff),iWork(iMemOff))
call SDCI_MRCI()

! Epilogue, end
!                                                                      *
!***********************************************************************
!                                                                      *
! Close open dafiles.

call DACLOS(LUSYMB)
call DACLOS(LUTRA)
call DACLOS(LUONE)
call DACLOS(LUREST)
call DACLOS(Lu_60)
call DACLOS(Lu_70)
call DACLOS(Lu_80)
call DACLOS(LUEIG)
call DACLOS(Lu_25)
call DACLOS(Lu_27)
!                                                                      *
!***********************************************************************
!                                                                      *
call FastIO('STATUS')
IRETURN = 0

return

end subroutine MRCI
