************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE MRCI(IRETURN)
************************************************************************
*  MULTI REFERENCE SDCI AND AVERAGE CPF PROGRAM.                       *
************************************************************************
C UNITS USED IN THE PROGRAM
C UNIT  5, INPUT
C UNIT  6, OUTPUT
C UNIT  2=LUPROP, (DA,ONEINT) FOR PROPERTY CALCULATIONS
C UNIT 10=LUSYMB, (DA,CIGUGA) SYMBOLIC FORMULAS
C UNIT 50=LUTRA, (DA,TRAINT) TRANSFORMED MO 2-EL INTEGRALS
C UNIT 60, (DA) SORTED AIBJ, ABIJ AND AIJK INTEGRALS
C UNIT 70, (DA) SORTED IJKL AND ABCI INTEGRALS
C UNIT 80, (DA) SORTED ABCD INTEGRALS
C UNIT 17=LUONE, (DA,TRAONE) ONE ELECTRON INTEGRALS
C UNIT 18=LUVEC, (Formatted, sequential!) MRCI ORBITALS OUT
C UNIT 23=LUEIG, (DA) WORKSPACE FOR MALMQUIST DIAGONALIZATION.
C UNIT 25, (DA) FOCK MATRIX AND DIAGONAL CSF MATRIX ELEMENTS
C UNIT 27, (DA) SCRATCH IN IIJJ. ALSO, REFERENCE CI VECTOR.
C UNIT 28=LUREST, (DA,MRCIVECT) CI VECTOR
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "WrkSpc.fh"

#include "SysDef.fh"

#include "mrci.fh"
*
*     Prologue, print program header
*
*      CALL SETTIM
      CALL XUFLOW
      Call qEnter('MRCI')
*
*     ( Workspace allocated in Start() )
*
*      Call IniMem
*PAM04      Call GetMem('WrkSpc','Max ','Real',MemOff,MaxMem)
*PAM04* PAM July 2004:
*PAM04* Actually allocate only half the memory, minus some spare overhead.
*PAM04* This is a temporary measure. I will systematically change the
*PAM04* present static allocation in order to use GETMEM instead.
*PAM04* Changed statement:
*PAM04*     MaxMem=MaxMem-3*1000
*PAM04      MaxMem=(MaxMem-3*1000)/2
*PAM04      Call GetMem('WrkSpc','Allo','Real',MemOff,MaxMem)
*PAM04      write(6,*)' Allocated ''WrkSpc''. memoff, maxmem=',memoff,maxmem

*PAM04 Now try completely without ''WrkSpc'' array:
      Call GetMem('HowMuch','Max ','Real',LDummy,MemTot)
*
*     Open files
*
      LUVEC=18
*
      LUSYMB=10
      CALL DANAME(LUSYMB,'CIGUGA')
      LUTRA=50
      CALL DANAME_MF(LUTRA,'TRAINT')
      LUONE=17
      CALL DANAME(LUONE,'TRAONE')
      LUREST=28
      CALL DANAME(LUREST,'MRCIVECT')
C Temporaries:
      Lu_60=60
      CALL DANAME_MF(Lu_60 ,'TIABIJ')
      Lu_70=70
      CALL DANAME_MF(Lu_70 ,'TIABCI')
      Lu_80=80
      CALL DANAME_MF(Lu_80 ,'TIABCD')
      LUEIG=23
      CALL DANAME(LUEIG,'FT23F001')
      Lu_25=25
      CALL DANAME(Lu_25 ,'FT25F001')
      Lu_27=27
      CALL DANAME(Lu_27 ,'FT27F001')
*
*     main body
*
*PAM04      iMemOff=ip_of_iWork(Work(MemOff))
*PAM04      CALL SDCI(Work(MemOff),iWork(iMemOff))
      CALL SDCI_MRCI()
*
*     Epilogue, end
*
*                                                                      *
************************************************************************
*                                                                      *
*     Close open dafiles.
*
      CALL DACLOS(LUSYMB)
      CALL DACLOS(LUTRA )
      CALL DACLOS(LUONE )
      CALL DACLOS(LUREST)
      CALL DACLOS(Lu_60 )
      CALL DACLOS(Lu_70 )
      CALL DACLOS(Lu_80 )
      CALL DACLOS(LUEIG )
      CALL DACLOS(Lu_25 )
      CALL DACLOS(Lu_27 )
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('MRCI')
      Call qStat(' ')
      Call FastIO('STATUS')
      IRETURN=0
      RETURN
      END
