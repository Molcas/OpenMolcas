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
* This routine get scalar integer data from the runfile.               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
*  Get_iScalar
*
*> @brief
*>   Get scalar data from runfile
*> @author Per-Olof Widmark
*>
*> @details
*> This routine gets scalar integer data from the runfile.
*>
*> @param[in]  Label Name of field
*> @param[out] Data  Data to get from runfile
*>
*> @see ::Put_iScalar
************************************************************************
      Subroutine Get_iScalar(Label,Data)
#include "pg_is_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Label
      Integer       Data
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Integer      ifirst,i
*
      Character*16 CmpLab
      DATA ifirst /0/
      SAVE ifirst

      If(ifirst.eq.0) then
         ifirst=1
         num_IS_init=0
         Do i=1,nTocIS
            iLbl_IS_inmem(i)=' '
            IS_init(i)=0
         End Do
      End If

      CmpLab=Label
      Call UpCase(CmpLab)

      Do i=1,num_IS_init
         If(iLbl_IS_inmem(i).eq.CmpLab) Then
           If(IS_init(i).ne.0) Then
             Data=i_IS_inmem(i)
             return
           End If
         End If
      End Do

      Call Get_iScalar_(Label,Data)
      num_IS_init=num_IS_init+1

      If(num_IS_init.gt.nTocIS) Then
#ifdef _DEBUG__
        Do i=1,num_IS_init
           write(6,*) iLbl_IS_inmem(i), IS_init(i), i_IS_inmem(i)
        End Do
#endif
        Call Abend()
      End If


      iLbl_IS_inmem(num_IS_init)=CmpLab
      IS_init(num_IS_init)=1
      i_IS_inmem(num_IS_init)=Data
      return
      End

      Subroutine Get_iScalar_(Label,Data)
      Implicit None
#include "pg_is_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Label
      Integer       Data
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Integer      RecVal(nTocIS)
      Character*16 RecLab(nTocIS)
      Integer      RecIdx(nTocIS)
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
      Call cRdRun('iScalar labels',RecLab,16*nTocIS)
      Call iRdRun('iScalar values',RecVal,nTocIS)
      Call iRdRun('iScalar indices',RecIdx,nTocIS)
*----------------------------------------------------------------------*
* Locate item                                                          *
*----------------------------------------------------------------------*
      item=-1
      CmpLab1=Label
      Call UpCase(CmpLab1)
      Do i=1,nTocIS
         CmpLab2=RecLab(i)
         Call UpCase(CmpLab2)
         If(CmpLab1.eq.CmpLab2) item=i
      End Do

      If(item.ne.-1) Then
         If(Recidx(item).eq.sSpecialField) Then
            Write(6,*) '***'
            Write(6,*) '*** Warning, reading temporary iScalar field'
            Write(6,*) '***   Field: ',Label
            Write(6,*) '***'
         End If
      End If
*----------------------------------------------------------------------*
* Transfer data to arguments                                           *
*----------------------------------------------------------------------*
      i_run_IS_used(item)=i_run_IS_used(item)+1
      If(item.eq.-1) Then
         Call SysAbendMsg('get_iScalar','Could not locate',Label)
      End If
      If(RecIdx(item).eq.0) Then
         Call SysAbendMsg('get_iScalar','Data not defined',Label)
      End If
      Data=RecVal(item)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
