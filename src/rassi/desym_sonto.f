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
      use stdalloc, only: mma_allocate, mma_deallocate
      use Symmetry_Info, only: nSym=>nIrrep, MUL
      IMPLICIT None
#include "rassi.fh"
      INTEGER SIZA,SYMLAB
      REAL*8 A(SIZA)
      REAL*8 B(NBST**2)

      REAL*8 ME
      REAL*8, Allocatable:: SCR(:)
      INTEGER ITD, IOF, ISY1, NB1, ISY2, NB2, ISY12_MA, I, J, IJ,
     &        ISY12_MA_BI, NB1_I, NB1_F, NB2_I, NB2_F, JI
      REAL*8 TDM

c Initialize
      B(:)=0.0D0

      Call mma_allocate(SCR,SIZA,Label='SCR')
      SCR(:)=0.0D0

c Diagonal symmetry blocks.
c Dont need to do anything, just leave it be
      IF(SYMLAB.EQ.1) THEN
        Call DCOPY_(SIZA,A(:),1,SCR,1)
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
                    SCR(IJ)=TDM
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
                  ME=SCR(ITD)
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
      Call mma_deallocate(SCR)

      END SUBROUTINE DESYM_SONTO
