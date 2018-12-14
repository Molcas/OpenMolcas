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
* Copyright (C) 2003, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This routine put scalar double data to the runfile.                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
*  Put_dScalar
*
*> @brief
*>   To add/update scalar data in runfile
*> @author Per-Olof Widmark
*>
*> @details
*> This routine is used to put scalar data of type
*> ``Real*8`` into the runfile. The data items are
*> identified by the \p label. Below is a list of the
*> data items that are recognized. The labels are
*> case insensitive and significant to 16 characters.
*>
*> For development purposes you can use an unsupported
*> label. Whenever such a field is accessed a warning
*> message is printed in the output, to remind the
*> developer to update this routine.
*>
*> List of known labels:
*>
*> - '``CASDFT energy``'  Energy for the last CASDFT calculation.
*> - '``CASPT2 energy``'  Energy for the last CASPT2 calculation.
*> - '``CASSCF energy``'  Energy for the last CASSCF calculation.
*> - '``Ener_ab``'
*> - '``KSDFT energy``'   Energy for the last KS-DFT calculation.
*> - '``Last energy``'    Last energy computed.
*> - '``PC Self Energy``' Self energy for point charges.
*> - '``PotNuc``'         Nuclear repusion energy.
*> - '``RF Self Energy``' Self energy in the Kirkwood model.
*> - '``SCF energy``'     Energy for the last SCF calculation.
*> - '``EThr``'           Energy convergence threshold.
*> - '``Thrs``'
*> - '``UHF energy``'
*> - '``DFT exch coeff``' Scaling factor for exchange terms of a density
*>   functional
*> - '``DFT corr coeff``' Scaling factor for correlation terms of a
*>   density functional
*>
*> @param[in] Label Name of field
*> @param[in] Data  Data to put on runfile
************************************************************************
      Subroutine Put_dScalar(Label,Data)
      Implicit None
#include "pg_ds_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Label
      Real*8        Data
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Real*8       RecVal(nTocDS)
      Character*16 RecLab(nTocDS)
      Integer      RecIdx(nTocDS)
      Save         RecVal
      Save         RecLab
      Save         RecIdx
*
      Character*16 CmpLab1
      Character*16 CmpLab2
      Integer      nData
      Integer      item
      Integer      iTmp
      Integer      i
*----------------------------------------------------------------------*
* Initialize local variables                                           *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* Do setup if this is the first call.                                  *
*----------------------------------------------------------------------*
*--- start pow mod ---
*     Write(6,'(3a)') 'Runfile: put_dscalar field "',Label,'"'
*--- end pow mod ---
      Call ffRun('dScalar labels',nData,iTmp)
      If(nData.eq.0) Then
         Do i=1,nTocDS
            RecLab(i)=' '
            RecVal(i)=0.0d0
            RecIdx(i)=sNotUsed
         End Do
*
*        Observe that label is at most 16 characters!
*
*                     1234567890123456
         RecLab(  1)='CASDFT energy   '
         RecLab(  2)='CASPT2 energy   '
         RecLab(  3)='CASSCF energy   '
         RecLab(  4)='Ener_ab         '
         RecLab(  5)='KSDFT energy    '
         RecLab(  6)='Last energy     '
         RecLab(  7)='PC Self Energy  '
         RecLab(  8)='PotNuc          '
         RecLab(  9)='RF Self Energy  '
         RecLab( 10)='SCF energy      '
         RecLab( 11)='Thrs            '
         RecLab( 12)='UHF energy      '
         RecLab( 13)='E_0_NN          '
         RecLab( 14)='W_or_el         '
         RecLab( 15)='W_or_Inf        '
         RecLab( 16)='EThr            '
         RecLab( 17)='Cholesky Thresho' !ld
         RecLab( 18)='Total Nuclear Ch' !arge
         RecLab( 19)='Numerical Gradie' !nt rDelta
         RecLab( 20)='MpProp Energy   '
         RecLab( 21)='UHFSPIN         '
         RecLab( 22)='S delete thr    '
         RecLab( 23)='T delete thr    '
         RecLab( 24)='MD_Etot0        '
         RecLab( 25)='MD_Time         '
         RecLab( 26)='LDF Accuracy    '
         RecLab( 27)='NAD dft energy  '
         RecLab( 28)='GradLim         '
         RecLab( 29)='StepFactor      '
         RecLab( 30)='Average energy  '
         RecLab( 31)='Timestep        '
         RecLab( 32)='MD_Etot         '
         RecLab( 33)='Max error       '
         RecLab( 34)='Total Charge    ' ! stores the total number of electrons in the computed system.
         RecLab( 35)='DFT exch coeff  '
         RecLab( 36)='DFT corr coeff  '
*                     1234567890123456
*
*        If u go beyond 64: update pg_ds_info.fh and this line!
         Call cWrRun('dScalar labels',RecLab,16*nTocDS)
         Call dWrRun('dScalar values',RecVal,nTocDS)
         Call iWrRun('dScalar indices',RecIdx,nTocDS)
      Else
         Call cRdRun('dScalar labels',RecLab,16*nTocDS)
         Call dRdRun('dScalar values',RecVal,nTocDS)
         Call iRdRun('dScalar indices',RecIdx,nTocDS)
      End If
*----------------------------------------------------------------------*
* Locate item                                                          *
*----------------------------------------------------------------------*
      item=-1
      CmpLab1=Label
      Call UpCase(CmpLab1)
      Do i=1,nTocDS
         CmpLab2=RecLab(i)
         Call UpCase(CmpLab2)
         If(CmpLab1.eq.CmpLab2) item=i
      End Do
*
* Do we create a new temporary field?
*
      If(item.eq.-1) Then
         Do i=1,nTocDS
            If(RecLab(i).eq.' ') item=i
         End Do
         If(item.ne.-1) Then
            RecLab(item)=Label
            RecIdx(item)=sSpecialField
            Call cWrRun('dScalar labels',RecLab,16*nTocDS)
            Call iWrRun('dScalar indices',RecIdx,nTocDS)
         End If
      End If
*
* Is this a temporary field?
*
      If(item.ne.-1) Then
         If(Recidx(item).eq.sSpecialField) Then
            Write(6,*) '***'
            Write(6,*) '*** Warning, writing temporary dScalar field'
            Write(6,*) '***   Field: ',Label
            Write(6,*) '***'
#ifdef _BIGOT_
            Call AbEnd()
#endif
         End If
      End If
*----------------------------------------------------------------------*
* Write data to disk.                                                  *
*----------------------------------------------------------------------*
      If(item.eq.-1) Then
         Call SysAbendMsg('put_dScalar','Could not locate',Label)
      End If
      RecVal(item)=Data
      Call dWrRun('dScalar values',RecVal,nTocDS)
      If(RecIdx(item).eq.0) Then
         RecIdx(item)=sRegularField
         Call iWrRun('dScalar indices',RecIdx,nTocDS)
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Do i=1,num_DS_init
         If(iLbl_DS_inmem(i).eq.CmpLab1) then
             i_DS_inmem(i)=Data
             DS_init(i)=1
             return
         End If
      End Do

      Return
      End
