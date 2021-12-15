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
* Copyright (C) 1993, Markus P. Fuelscher                              *
*               1995, Martin Schuetz                                   *
************************************************************************
      Subroutine RdMBPT(ipCMO,lthCMO,ipEOrb,lthEOr)
************************************************************************
*                                                                      *
*     Read the MBPTOUT file genereated by the SCF program              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*     modified by:                                                     *
*     M.G. Schuetz                                                     *
*     University of Lund, Sweden, 1995                                 *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)

*     declaration of calling arguments
      Integer ipCMO,ipEOrb,lthCMO,lthEOr

#include "real.fh"
#include "mxdim.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "mbpt2aux.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "files_mbpt2.fh"
      Real*8, Allocatable:: CMO_t(:)

*     declaration of local variables...
      Logical Debug
      Data Debug/.False./


#include "SysDef.fh"

*
*...  Read nSym, nBas, nOrb, nOcc, nFro, CMO and orbital energies from COMFILE
*
         Call Get_iScalar('nSym',nSym)
         Call Get_iArray('nBas',nBas,nSym)
         Call Get_iArray('nOrb',nOrb,nSym)
         Call Get_iArray('nIsh',nOcc,nSym)
         Call Get_iArray('nFro',nFro,nSym)
         If (Debug) Then
            Write (6,'(A,8I5)') 'nSym:',nSym
            Write (6,'(A,8I5)') 'nBas:',(nBas(i),i=1,nSym)
            Write (6,'(A,8I5)') 'nOrb:',(nOrb(i),i=1,nSym)
            Write (6,'(A,8I5)') 'nOcc:',(nOcc(i),i=1,nSym)
            Write (6,'(A,8I5)') 'nFro:',(nFro(i),i=1,nSym)
         End If
         lthCMO=0
         Do iSym=1,nSym
            If (nFro(iSym).ne.0) Then
               Write (6,*) 'Some orbitals where frozen in the SCF!'
               Call Abend()
            End If
            nDel(iSym)=nBas(iSym)-nOrb(iSym)
            nExt(iSym)=nOrb(iSym)-nOcc(iSym)
            nDsto(iSym)=nDel(iSym)
            lthCMO = lthCMO + nBas(iSym)*nOrb(iSym)
         End Do
*
         Call mma_allocate(CMO_t,lthCMO,Label='CMO_t')
         Call Get_CMO_(CMO_t,lthCMO)
         Call GetMem('CMO   ','Allo','Real',ipCMO,lthCMO)
*
*...  set MO coefficients of the deleted orbitals to zero
*     Observe that these are not included at all in the basis
         iStart   = ipCMO
         iStart_t = 1
         Do iSym=1,nSym
            call dcopy_(nOrb(iSym)*nBas(iSym),CMO_t(iStart_t),1,
     &                                       Work(iStart),1)
            iStart   = iStart   + nOrb(iSym)*nBas(iSym)
            iStart_t = iStart_t + nOrb(iSym)*nBas(iSym)
            call dcopy_((nBas(iSym)-nOrb(iSym))*nBas(iSym),
     &                  [Zero],0,Work(iStart),1)
            iStart   = iStart   +
     &                 (nBas(iSym)-nOrb(iSym))*nBas(iSym)
         End Do
         Call mma_deallocate(CMO_t)
*
         Call Get_OrbE_(ipEOrb_t,lthEOr)
         nnB=lthEOr
         Call GetMem('EOrb  ','Allo','Real',ipEOrb,lthEOr)
*
*...  set energies of the deleted orbitals to zero
*
         iStart   = ipEOrb
         iStart_t = ipEOrb_t
         Do iSym=1,nSym
            call dcopy_(nOrb(iSym),Work(iStart_t),1,Work(iStart),1)
            iStart  = iStart   + nOrb(iSym)
            iStart_t= iStart_t + nOrb(iSym)
*
            call dcopy_(nBas(iSym)-nOrb(iSym),[Zero],0,Work(iStart),1)
            iStart  = iStart   + nBas(iSym)-nOrb(iSym)
         End Do
         Call GetMem('EOrb_t','Free','Real',ipEOrb_t,lthEOr)
*
      Return
      End
