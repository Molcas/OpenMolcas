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
* Copyright (C) 2021, Rulin Feng                                       *
************************************************************************
*       ****************************************************
*                 Desymmetrize a given array(matrix)
*       ****************************************************
*        This routine is made to expand a given symmetry-adapted array
*        into a C1 symmetry. SYMLAB is a bit flag,
*        e.g., for symmetry 3, symlab=4=2^(3-1). A is input array with
*        size SIZA, B is output array with size of NBST**2.
*        Not tested for general cases, so for use of SO-NTOs only.
*
*                                                      -RF 8/24,2021
      SUBROUTINE DESYM_SONTO(A,SIZA,B,SYMLAB)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
#include "rassi.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='DESYM_SONTO')
      INTEGER SIZA,SYMLAB
      REAL*8 ME
      REAL*8 A(SIZA)
      REAL*8 B(NBST**2)

c Initialize
      Call DCOPY_(NBST**2,[0.0D0],0,B,1)

      Call GETMEM('SCR','ALLO','REAL',LSCR,SIZA)
      Call DCOPY_(SIZA,[0.0D00],0,WORK(LSCR),1)

c Diagonal symmetry blocks.
c Dont need to do anything, just leave it be
      IF(SYMLAB.EQ.1) THEN
        Call DCOPY_(SIZA,A(:),1,WORK(LSCR),1)
c Non-diagonal symmetry blocks
c note that only half of the total matrix has been stored
      ELSE
        ITD=0
        IOF=0
        Do ISY1=1,NSYM
          NB1=NBASF(ISY1)
          Do ISY2=1,NSYM
            NB2=NBASF(ISY2)
            ISY12_ma=MUL(ISY1,ISY2)
            ISY12_ma_bi=2**(ISY12_ma-1)
            If(ISY12_ma_bi.EQ.SYMLAB) then
              IF(ISY1.GT.ISY2) THEN
                Do J=1,NB2
                  Do I=1,NB1
                    ITD=ITD+1
                    TDM=A(ITD)
                    IJ=IOF+J+NB2*(I-1)
                    WORK(LSCR-1+IJ)=TDM
                  Enddo
                Enddo
                IOF=IOF+NB1*NB2
              Endif
            Endif
          Enddo
        Enddo
      Endif
c Expand into C1

      ITD=0
      NB1_i=0
      NB1_f=0
      Do ISY1=1,NSYM
        NB1=NBASF(ISY1)
        NB1_f=NB1_i+NB1
        NB2_i=0
        NB2_f=0
        Do ISY2=1,NSYM
          ISY12_ma=MUL(ISY1,ISY2)
          ISY12_ma_bi=2**(ISY12_ma-1)
          NB2=NBASF(ISY2)
          NB2_f=NB2_i+NB2
          If(ISY12_ma_bi.EQ.SYMLAB) then
            Do J=NB2_i+1,NB2_f
              Do I=NB1_i+1,NB1_f
                If(I.LE.J) then
                  ITD=ITD+1
                  ME=WORK(LSCR-1+ITD)
                  IJ=I+NBST*(J-1)
                  JI=J+NBST*(I-1)
                  B(IJ)=ME
                  B(JI)=ME
                Endif
              Enddo
            Enddo
          Endif
          NB2_i=NB2_i+NB2
        Enddo
        NB1_i=NB1_i+NB1
      Enddo
      Call GETMEM('SCR','FREE','REAL',LSCR,SIZA)

      RETURN
      END
