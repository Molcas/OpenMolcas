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
* This routine query the existence of array data on runfile.           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
*  Qpg_cArray
*
*> @brief
*>   Check if a field exist on the runfile
*> @author Per-Olof Widmark
*>
*> @details
*> This routine query the existence of array data on runfile.
*>
*> @param[in]  Label Name of field
*> @param[out] Found Was the field found
*> @param[out] nData Length of field
*>
*> @see ::Put_cArray
************************************************************************
      Subroutine Qpg_cArray(Label,Found,nData)
      Implicit None
#include "pg_ca_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Label
      Logical       Found
      Integer       nData
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Character*16 RecLab(nTocCA)
      Integer      RecIdx(nTocCA)
      Integer      RecLen(nTocCA)
*
      Character*16 CmpLab1
      Character*16 CmpLab2
      Integer      item
      Integer      nTmp,iTmp
      Integer      i
*----------------------------------------------------------------------*
* Initialize local variables                                           *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* Read info from runfile.                                              *
*----------------------------------------------------------------------*
      Call ffRun('cArray labels',nTmp,iTmp)
      If(nTmp.eq.0) Then
         Found=.False.
         nData=0
         Return
      End If
      Call cRdRun('cArray labels',RecLab,16*nTocCA)
      Call iRdRun('cArray indices',RecIdx,nTocCA)
      Call iRdRun('cArray lengths',RecLen,nTocCA)
*----------------------------------------------------------------------*
* Locate item                                                          *
*----------------------------------------------------------------------*
      item=-1
      CmpLab1=Label
      Call UpCase(CmpLab1)
      Do i=1,nTocCA
         CmpLab2=RecLab(i)
         Call UpCase(CmpLab2)
         If(CmpLab1.eq.CmpLab2) item=i
      End Do
*
* Do we read an old temporary field?
*
      If(item.ne.-1) Then
         If(RecIdx(item).eq.sSpecialField) Then
            Write(6,*) '***'
            Write(6,*) '*** Warning, querying temporary cArray field'
            Write(6,*) '***   Field: ',Label
            Write(6,*) '***'
#ifdef _BIGOT_
            Call AbEnd()
#endif
         End If
      End If
*----------------------------------------------------------------------*
* Did we manage to find it?                                            *
*----------------------------------------------------------------------*
      Found=.true.
      If(item.eq.-1) Found=.false.
      If(item.ne.-1.and.RecIdx(item).eq.0) Found=.false.
      If(Found) Then
         nData=RecLen(item)
      Else
         nData=0
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
