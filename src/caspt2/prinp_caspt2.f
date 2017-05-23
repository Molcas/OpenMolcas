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
* Copyright (C) 1994, Markus P. Fuelscher                              *
*               1994, Per Ake Malmqvist                                *
************************************************************************
      Subroutine PrInp_CASPT2
************************************************************************
*                                                                      *
*     purpose:                                                         *
*     - echo the input parameters                                      *
*                                                                      *
*     calling parameters: none                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and P.-AA. Malmqvist                              *
*     University of Lund, Sweden, 1994                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
      Character(8)   Fmt1,Fmt2
      Character(120)  Line,BlLine,StLine
      Character(3) lIrrep(8)
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
      Call qEnter('PrInp_CASPT2')
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*
      Line=' '
      lLine=Len(Line)
      Do i=1,lLine
        BlLine(i:i)=' '
        StLine(i:i)='*'
      End Do
      lPaper=132
      left=(lPaper-lLine)/2
      WRITE(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
      WRITE(Fmt2,'(A,I3.3,A)') '(',left,'X,'
*----------------------------------------------------------------------*
*     Print the ONEINT file identifier                                 *
*----------------------------------------------------------------------*
      IF(IPRGLB.GE.VERBOSE) THEN
      WRITE(6,*)
      WRITE(6,Fmt1) 'Header of the ONEINT file:'
      WRITE(6,Fmt1) '--------------------------'
      WRITE(Line,'(36A2)') (Header(i),i=1,36)
      Call LeftAd(Line)
      WRITE(6,Fmt1)  Line(:mylen(Line))
      WRITE(Line,'(36A2)') (Header(i),i=37,72)
      Call LeftAd(Line)
      WRITE(6,Fmt1)  Line(:mylen(Line))
      WRITE(6,*)
      END IF
*----------------------------------------------------------------------*
*     Print cartesian coordinates of the system                        *
*----------------------------------------------------------------------*
      IF(IPRGLB.GE.VERBOSE) THEN
      Call PrCoor
      END IF
*----------------------------------------------------------------------*
*     Print orbital and wavefunction specifications                    *
*----------------------------------------------------------------------*
      IF(IPRGLB.GE.USUAL) THEN
      WRITE(6,*)
      Line=' '
      WRITE(Line(left-2:),'(A)') 'Wave function specifications:'
      CALL CollapseOutput(1,Line)
      WRITE(6,Fmt1)'-----------------------------'
      WRITE(6,*)
      WRITE(6,Fmt2//'A,T45,I6)')'Number of closed shell electrons',
     &                           2*NISHT
      WRITE(6,Fmt2//'A,T45,I6)')'Number of electrons in active shells',
     &                           NACTEL
      WRITE(6,Fmt2//'A,T45,I6)')'Max number of holes in RAS1 space',
     &                           NHOLE1
      WRITE(6,Fmt2//'A,T45,I6)')'Max number of electrons in RAS3 space',
     &                           NELE3
      WRITE(6,Fmt2//'A,T45,I6)')'Number of inactive orbitals',
     &                           NISHT
      WRITE(6,Fmt2//'A,T45,I6)')'Number of active orbitals',
     &                           NASHT
      WRITE(6,Fmt2//'A,T45,I6)')'Number of secondary orbitals',
     &                           NSSHT
      WRITE(6,Fmt2//'A,T45,F6.1)')'Spin quantum number',
     &                           0.5D0*DBLE(ISPIN-1)
      WRITE(6,Fmt2//'A,T45,I6)')'State symmetry',
     &                           LSYM
      WRITE(6,Fmt2//'A,T40,I11)')'Number of CSFs',
     &                           NCONF
      WRITE(6,Fmt2//'A,T45,I6)')'Number of root(s) available',
     &                           NROOTS
      WRITE(6,Fmt2//'A,T45,I6)')'Root passed to geometry opt.',
     &                           iRlxRoot
      IF(IFMIX) THEN
        WRITE(6,Fmt2//'A,T45,10I3)')'A file JOBMIX will be created.'
      END IF
      IF(NSTATE.GT.1) THEN
        WRITE(6,Fmt1) 'This is a MULTI-STATE CASSCF reference'
        WRITE(6,Fmt2//'A,T45,I6)')'Number of CI roots used',
     &                           NSTATE
        WRITE(6,Fmt2//'A,(T45,10I3))')'These are:',
     &                           (MSTATE(I),I=1,NSTATE)
      ELSE
        If ( ISCF.eq.0 ) then
           WRITE(6,Fmt1) 'This is a CASSCF or RASSCF reference function'
        Else If ( ISCF.eq.1 ) then
           WRITE(6,Fmt1) 'This is a closed shell RHF reference function'
        Else
           WRITE(6,Fmt1)
     &     'This is a high spin open shell RHF reference function'
        End If
      END IF
      CALL CollapseOutput(0,'Wave function specifications:')
      END IF
*
      Call Get_cArray('Irreps',lIrrep,24)
      Do iSym = 1, nSym
         Call RightAd(lIrrep(iSym))
      End Do
*
      IF(IPRGLB.GE.USUAL  ) THEN
      WRITE(6,*)
      Line=' '
      WRITE(Line(left-2:),'(A)') 'Orbital specifications:'
      CALL CollapseOutput(1,Line)
      WRITe(6,Fmt1)'-----------------------'
      WRITE(6,*)
      WRITE(6,Fmt2//'A,T47,8I4)') 'Symmetry species',
     &                            (iSym,iSym=1,nSym)
      WRITE(6,Fmt2//'A,T47,8(1X,A))') '                ',
     &                            (lIrrep(iSym),iSym=1,nSym)
      WRITE(6,Fmt2//'A,T47,8I4)') 'Frozen orbitals',
     &                            (nFro(iSym),iSym=1,nSym)
      WRITE(6,Fmt2//'A,T47,8I4)') 'Inactive orbitals',
     &                            (nIsh(iSym),iSym=1,nSym)
      WRITE(6,Fmt2//'A,T47,8I4)') 'Active orbitals',
     &                            (nAsh(iSym),iSym=1,nSym)
      WRITE(6,Fmt2//'A,T47,8I4)') 'Secondary orbitals',
     &                            (nSsh(iSym),iSym=1,nSym)
      WRITE(6,Fmt2//'A,T47,8I4)') 'Deleted orbitals',
     &                            (nDel(iSym),iSym=1,nSym)
      WRITE(6,Fmt2//'A,T47,8I4)') 'Number of basis functions',
     &                            (nBas(iSym),iSym=1,nSym)
      CALL CollapseOutput(0,'Orbital specifications:')
      END IF
*----------------------------------------------------------------------*
*     Print routing information                                        *
*----------------------------------------------------------------------*
      IF(IPRGLB.GE.USUAL  ) THEN
      If ( RFpert ) then
         WRITE(6,*)
         WRITE(6,Fmt1)'Reaction field specifications:'
         WRITE(6,Fmt1)'------------------------------'
         WRITE(6,*)
         WRITE(6,'(6X,A)')'An external reaction field was determined'//
     &     ' previously and added to the one-electron hamiltonian.'
         WRITE(6,'(6X,A)')'It will not be reevaluated even though'//
     &     ' dynamic correlation may change the density.'
         WRITE(6,*)
      End If
      WRITE(6,*)
      WRITE(6,Fmt2//'A,A)')'Type of Fock operator: ',TRIM(FOCKTYPE)
      WRITE(6,Fmt2//'A,A)')'Type of 0-order Hamiltonian: ',TRIM(HZERO)
      WRITE(6,Fmt2//'A,F9.2)')'  "IP-EA" denominator shift = ',BSHIFT
      If ( SHIFT.ne.0.0d0 ) then
         WRITE(6,Fmt2//'A,F9.2)')
     &     'Extra denominator shift = ',SHIFT
      End if
      If ( SHIFTI.ne.0.0d0 ) then
         WRITE(6,Fmt2//'A,F9.2)')
     &     'Extra imaginary denominator shift = ',SHIFTI
      End if
      If ( ORBIN.eq.'TRANSFOR' ) then
         WRITE(6,Fmt1)
     &     'The input orbitals and the CI vector will be transformed.'
      Else
         WRITE(6,Fmt1)
     &'The input orbitals and the CI vector will be used as they are.'
      End If
      WRITE(6,*)
C Compute necessary quantities for subsequent gradient calc:
      IF(IFDENS) THEN
        WRITE(6,*)' The wave functions (P_CAS) H W |Psi_0> and '//
     &       '(P_CAS) W H |Psi_0>'
        WRITE(6,*)' will be computed and written out for subsequent'
        WRITE(6,*)' use in a gradient calculation.'
      END IF
      END IF

      Call qExit('PrInp_CASPT2')
      Return
      End
