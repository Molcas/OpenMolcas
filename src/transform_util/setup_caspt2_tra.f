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
      Subroutine SetUp_CASPT2_Tra(nSym_,nBas_,nOrb_,nIsh_,nAsh_,
     &                            nFro_,nDel_,ipCMO,lthCMO,
     &                            LuIntM_,LuHlf1_,LuHlf2_,LuHlf3_)
      Implicit Real*8 (a-h,o-z)

#include "rasdim.fh"
#include "caspt2.fh"
*
      Integer nBas_(8),nOrb_(8),nAsh_(8),nIsh_(8), nFro_(8),nDel_(8)
*                                                                      *
************************************************************************
*                                                                      *
      nSym=nSym_
      Do i = 1, nSym
         nBas(i) = nBas_(i)
         nOrb(i) = nOrb_(i)
         nFro(i) = nFro_(i)
         nDel(i) = nDel_(i)
         nAsh(i) = nAsh_(i)
         nIsh(i) = nIsh_(i)
         nOsh(i) = nAsh_(i) + nIsh_(i)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Do i = 1, 8
         Do j = 1, 8
            Mul(i,j)=iEor(i-1,j-1)+1
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      LCMO=ipCMO
      nCMO=lthCMO
*                                                                      *
************************************************************************
*                                                                      *
*---- Open the temporary files here! This is where they get their final
*     unit numbers.
*
*
      CALL DANAME_MF_wa(LuHlf1_,'LUHLF1')
      CALL DANAME_MF_wa(LuHlf2_,'LUHLF2')
      CALL DANAME_MF_wa(LuHlf3_,'LUHLF3')
      LuHlf1=LuHlf1_
      LuHlf2=LuHlf2_
      LuHlf3=LuHlf3_
*
*     Observe that LuIntM should be opened prior to this call!
*
      LuIntM=LuIntM_
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
