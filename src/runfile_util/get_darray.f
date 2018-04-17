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
* This routine get array double data from the runfile.                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
*  Get_dArray
*
*> @brief
*>   Read array data from runfile
*> @author Per-Olof Widmark
*>
*> @details
*> This routine gets array double data from the runfile.
*>
*> @param[in]  Label Name of field
*> @param[out] Data  Data to read from runfile
*> @param[in]  nData Length of array
*>
*> @see ::Put_dArray
************************************************************************
      Subroutine Get_dArray(Label,Data,nData)
      Implicit None
#include "pg_da_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Label
      Integer       nData
      Real*8        Data(nData)
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Character*16 RecLab(nTocDA)
      Integer      RecIdx(nTocDA)
      Integer      RecLen(nTocDA)
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
      Call cRdRun('dArray labels',RecLab,16*nTocDA)
      Call iRdRun('dArray indices',RecIdx,nTocDA)
      Call iRdRun('dArray lengths',RecLen,nTocDA)
*----------------------------------------------------------------------*
* Locate item                                                          *
*----------------------------------------------------------------------*
      item=-1
      CmpLab1=Label
      Call UpCase(CmpLab1)
      Do i=1,nTocDA
         CmpLab2=RecLab(i)
         Call UpCase(CmpLab2)
         If(CmpLab1.eq.CmpLab2) item=i
      End Do

      If(item.ne.-1) Then
         If(Recidx(item).eq.sSpecialField) Then
            Write(6,*) '***'
            Write(6,*) '*** Warning, reading temporary dArray field'
            Write(6,*) '***   Field: ',Label
            Write(6,*) '***'
         End If
      End If
*----------------------------------------------------------------------*
* Transfer data to arguments                                           *
*----------------------------------------------------------------------*
      i_run_DA_used(item)=i_run_DA_used(item)+1
      If(item.eq.-1) Then
         Call SysAbendMsg('get_dArray','Could not locate:',Label)
      End If
      If(RecIdx(item).eq.0) Then
         Call SysAbendMsg('get_dArray','Data not defined:',Label)
      End If
      If(Reclen(item).ne.nData) Then
         Call SysAbendMsg('get_dArray','Data of wrong length:',Label)
      End If
      Call dRdRun(RecLab(item),Data,nData)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
