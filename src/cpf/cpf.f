************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1986, Per E. M. Siegbahn                               *
*               1986, Margareta R. A. Blomberg                         *
************************************************************************
************************************************************************
*                                                                      *
* PER SIEGBAHN                                                         *
* MARGARETA BLOMBERG                                                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY                                  *
* UNIVERSITY OF LUND                                                   *
* SWEDEN                                                               *
*                                                                      *
************************************************************************
      SUBROUTINE CPF(IRETURN)
************************************************************************
*                                                                      *
*                                C P F                                 *
* MODIFIED TO IBM BY ROLAND LINDH 02/17/88                             *
* MODIFIED TO MOLCAS-2 BY ROLAND LINDH 03/26/91                        *
* MODIFIED TO MOLCAS-3 BY M.P. FUELSCHER 08/31/93                      *
* MODIFIED TO MOLCAS-4 BY P.A. MALMQVIST AND N.W.MORTIARTY 10/25/96    *
* MODIFIED TO MOLCAS 4.1 BY R. LINDH 02/24/98 (Multi fileing)          *
************************************************************************
C
C     UNITS USED IN THE PROGRAM
C     UNIT  5 , INPUT
C     UNIT  6 , OUTPUT
C     UNIT 10 , SYMBOLIC FORMULAS
C     UNIT 50 , TRANSFORMED MO 2-EL INTEGRALS
C     UNIT 60 , SORTED AIBJ, ABIJ AND AIJK INTEGRALS
C     UNIT 70 , SORTED IJKL AND ABCI INTEGRALS
C     UNIT 80 , SORTED ABCD INTEGRALS
C     UNIT 17 , ONE ELECTRON INTEGRALS
C     UNIT 19 , (Formatted sequential!) CPF-ORBITALS OUT
C     UNIT 25 , FOCK MATRIX AND DIAGONAL CSF MATRIX ELEMENTS
C     UNIT 26 , CI VECTOR
C     UNIT 27 , SCRATCH IN IIJJ
C     UNIT 30 ,
C
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "files_cpf.fh"
#include "WrkSpc.fh"
*
*     Prologue
*
      CALL QENTER('CPFMCPF')
*     CALL SETTIM
C     CALL HELLO
*
*     (Workspace allocated in Start() )
*
      Call GetMem('WrkSpc','Max ','Real',MemOff,MEMORY)
      MEMORY=INT(MEMORY*0.80D0)
      Call GetMem('WrkSpc','Allo','Real',MemOff,MEMORY)
*
*     Open files
*
      Lu_CIGuga=10
      CALL DANAME(Lu_CIGuga,'CIGUGA')
      Lu_TraInt=50
      CALL DANAME_MF(Lu_TraInt,'TRAINT')
      Lu_TraOne=17
      CALL DANAME(Lu_TraOne,'TRAONE')
      Lu_CI=26
      CALL DANAME(Lu_CI,'CPFVECT')
      Lu_CPFORB=19
C Temporaries:
      Lu_TiABIJ=60
      CALL DANAME_MF(Lu_TiABIJ,'TIABIJ')
      Lu_TiABCI=70
      CALL DANAME_MF(Lu_TiABCI,'TIABCI')
      Lu_TiABCD=80
      CALL DANAME_MF(Lu_TiABCD,'TIABCD')
      Lu_25=25
      CALL DANAME(Lu_25,'FT25F001')
      Lu_27=27
      CALL DANAME(Lu_27,'FT27F001')
      Lu_30=30
      CALL DANAME(Lu_30    ,'FT30F001')
*
*     Body
*
      iMemOff=ip_of_iWork_d(Work(MemOff))
      CALL SDCI_CPF(Work(MemOff),iWork(iMemOff),MEMORY)
*
*     Deallocate the workspace
*
      Call GetMem('WrkSpc','Free','Real',MemOff,MEMORY)
*
*     Epilogue, end
*
*                                                                      *
************************************************************************
*                                                                      *
*     Close open dafiles
*
      CALL DACLOS(Lu_CIGuga)
      CALL DACLOS(Lu_TraInt)
      CALL DACLOS(Lu_TraOne)
      CALL DACLOS(Lu_CI)
      CALL DACLOS(Lu_TiABIJ)
      CALL DACLOS(Lu_TiABCI)
      CALL DACLOS(Lu_TiABCD)
      CALL DACLOS(Lu_25)
      CALL DACLOS(Lu_27)
      CALL DACLOS(Lu_30)
*                                                                      *
************************************************************************
*                                                                      *
      CALL QEXIT('CpfMcpf')
      CALL QSTAT(' ')
      CALL FASTIO('STATUS')
      IRETURN=0
      RETURN
      END
