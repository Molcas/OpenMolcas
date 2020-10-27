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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      Subroutine CSF2SD(CSF,SD,is)
      use ipPage, only: Diskbased
*
*  Transforms a CSF vector to slater determinants
*
      Use Arrays, only: DTOC, CNSM
      implicit Real*8(a-h,o-z)
#include "detdim.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "cicisp_mclr.fh"
#include "Input.fh"
#include "spinfo_mclr.fh"
*
      Real*8 CSF(*),SD(*)
      Real*8, Allocatable:: CTM(:)
*

      iiCOPY=0
      iprdia=0
      nConf=Max(ncsf(is),ndtasm(iS))
      isym=iEor(is-1,State_Sym-1)+1
      i=2
      If (isym.eq.1) i=1

      If (diskbased) Then
         CALL CSDTVC_MCLR(CSF,SD,1,DTOC,CNSM(i)%ICTS,IS,iiCOPY,IPRDIA)
      Else
         Call mma_allocate(CTM,nConf,Label='CTM')
         CTM(:)=Zero
         CTM(1:ncsf(is))=CSF(1:ncsf(is))

         CALL CSDTVC_MCLR(CTM,SD,1,DTOC,CNSM(i)%ICTS,IS,iiCOPY,IPRDIA)

         Call mma_deallocate(CTM)
      End If
*
      Return
      End
