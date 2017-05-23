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
* Copyright (C) 2014, Naoki Nakatani                                   *
************************************************************************
      Subroutine BLOCK_DENSI_RASSCF(jRoot,D,DS,PS,PA,PT)
      Implicit Real*8 (A-H,O-Z)
* Load 2RDM from Block DMRG and compute D, DS, PS, and PA
* Written by N. Nakatani, Oct. 2014
*
* D   : 1-El density matrix
*
* DS  : spin density matrix
*
* PS  : symmetrized 2-El density matrix
*
* PA  : anti-symmetrized 2-El density matrix
*
* PT  : working space for 2-El density matrix (NAC**4)
*
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
*
      Dimension D(NACPAR),DS(NACPAR),PS(NACPR2),PA(NACPR2)
      Dimension PT(NAC,NAC,NAC,NAC)
*
      Call DCOPY_(NACPAR,0.0D0,0,D, 1)
      Call DCOPY_(NACPAR,0.0D0,0,DS,1)
      Call DCOPY_(NACPR2,0.0D0,0,PS,1)
      Call DCOPY_(NACPR2,0.0D0,0,PA,1)
*
      If(NACTEL.GT.1) Then
        Call block_load2pdm(NAC,PT,jRoot,jRoot)
* No spin density is obtained from Spin-adapted DMRG (Block's default)
* TODO: 1-El density matrix should be computed with spin-orbital index
*       to store spin density with Block code...
        IJ_pack=1
        Do J=1,NAC
          Do I=1,J
            D1sum=0.0D0
            Do K=1,NAC
              D1sum=D1sum+PT(K,K,I,J)
            End Do
            D(IJ_pack)=D1sum/(NACTEL-1)
*           DS(IJ_pack)=0.0D0
            IJ_pack=IJ_pack+1
          End Do
        End Do

        IJKL_pack=0
        Do I=1,NAC
          Do J=1,I
            Do K=1,I
              LLIM=K
              If(K.EQ.I)LLIM=J
              DO L=1,LLIM
                IJKL_pack=IJKL_pack+1
                If(K.EQ.L) Then
                  PS(IJKL_pack)=0.5D0*PT(L,K,J,I)
*                 PA(IJKL_pack)=0.0D0
                Else
                  PS(IJKL_pack)=0.5D0*(PT(L,K,J,I)+PT(K,L,J,I))
                  PA(IJKL_pack)=0.5D0*(PT(L,K,J,I)-PT(K,L,J,I))
                End If
              End Do
            End Do
          End Do
        End Do
      Else
* special case for NACTEL = 1
        Call block_load1pdm(NAC,PT,jRoot,jRoot)
        IJ_pack=1
        Do J=1,NAC
          Do I=1,J
            D(IJ_pack)=PT(I,J,1,1)
            IJ_pack=IJ_pack+1
          End Do
        End Do
      End If

      RETURN
      END
