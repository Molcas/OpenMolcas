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
* This routine get scalar double data from the runfile.                *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
*  Get_dScalar
*
*> @brief
*>   Get scalar data form runfile
*> @author Per-Olof Widmark
*>
*> @details
*> This routine gets scalar double data from the runfile.
*>
*> @param[in]  Label Name of field
*> @param[out] Data  Data to get from runfile
*>
*> @see ::Put_dScalar
************************************************************************
      Subroutine Get_dScalar(Label,Data)
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
      Character*16 CmpLab

      Integer      dfirst,i

      DATA dfirst /0/
      SAVE dfirst

      If(dfirst.eq.0) Then
         dfirst=1
         num_DS_init=0
         Do i=1,nTocDS
            iLbl_DS_inmem(i)=' '
            DS_init(i)=0
         End Do
      End If

      CmpLab=Label
      Call UpCase(CmpLab)

      Do i=1,num_DS_init
         If(iLbl_DS_inmem(i).eq.CmpLab) Then
           If(DS_init(i).ne.0) Then
             Data=i_DS_inmem(i)
             return
           End If
         End If
      End Do

      Call Get_dScalar_(Label,Data)
      num_DS_init=num_DS_init+1

      If(num_DS_init.gt.nTocDS) Then
#ifdef _DEBUG_
        Do i=1,num_DS_init
          Write(6,*) iLbl_DS_inmem(i), DS_init(i), i_DS_inmem(i), CmpLab
        End Do
#endif
        Call Abend()
      End If

      iLbl_DS_inmem(num_DS_init)=CmpLab
      DS_init(num_DS_init)=1
      i_DS_inmem(num_DS_init)=Data
      return
      End

      Subroutine Get_dScalar_(Label,Data)
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
*
      Character*16 CmpLab1
      Character*16 CmpLab2
      Integer      item
      Integer      i
*----------------------------------------------------------------------*
* Initialize local variables                                           *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* Read info from runfile.                                              *
*----------------------------------------------------------------------*
      Call cRdRun('dScalar labels',RecLab,16*nTocDS)
      Call dRdRun('dScalar values',RecVal,nTocDS)
      Call iRdRun('dScalar indices',RecIdx,nTocDS)
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

      If(item.ne.-1) Then
         If(Recidx(item).eq.sSpecialField) Then
            Write(6,*) '***'
            Write(6,*) '*** Warning, reading temporary dScalar field'
            Write(6,*) '***   Field: ',Label
            Write(6,*) '***'
#ifdef _BIGOT_
            Call AbEnd()
#endif
         End If
      End If
*----------------------------------------------------------------------*
* Transfer data to arguments                                           *
*----------------------------------------------------------------------*
      i_run_DS_used(item)=i_run_DS_used(item)+1
      If(item.eq.-1) Then
         Call SysAbendMsg('get_dScalar','Could not locate',Label)
      End If
      If(RecIdx(item).eq.0) Then
         Call SysAbendMsg('get_dScalar','Data not defined',Label)
      End If
      Data=RecVal(item)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
