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
* Copyright (C) 2016, Sebastian Wouters                                *
*               2016, Quan Phung                                       *
************************************************************************
! CheMPS2-Molcas main interface
! Based on Block interface, writen by N. Nakatani
! Written by Quan Phung and Sebastian Wouters, Leuven, Aug 2016
      Subroutine CHEMPS2_DENSI_RASSCF(jRoot,D,DS,PS,PA,PT)
      Implicit Real*8 (A-H,O-Z)
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
        Call chemps2_load2pdm(NAC,PT,jRoot)
        IJ_pack=1
        Do J=1,NAC
          Do I=1,J
            D1sum=0.0D0
            Do K=1,NAC
              D1sum=D1sum+PT(K,K,I,J)
            End Do
            D(IJ_pack)=D1sum/(NACTEL-1)
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
         write(6,*) 'CheMPS2 does not allow 1 electron.'
      End If

      RETURN
      END
