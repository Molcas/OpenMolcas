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
* This routine put array integer data to the runfile.                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
*  Put_iArray
*
*> @brief
*>   Add/update array data in runfile
*> @author Per-Olof Widmark
*>
*> @details
*> This routine is used to put array data of type
*> ``Integer`` into the runfile. The data items are
*> identified by the label. Below is a list of the
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
*> - '``Center Index``'
*> - '``Ctr Index Prim``'       Idem with primitive basis set.
*> - '``nAsh``'                 The number of active orbitals per irreducible representation.
*> - '``nBas``'                 The number of basis functions per irreducible representation.
*> - '``nDel``'                 The number of deleted orbitals per irreducible representation.
*> - '``nFro``'                 The number of frozen orbitals per irreducible representation, i.e. orbitals that are not optimized.
*> - '``nIsh``'                 The number of inactive orbitals per irreducible representation.
*> - '``nIsh beta``'
*> - '``nOrb``'                 The total number of orbitals per irreducible representation.
*> - '``Orbital Type``'
*> - '``Slapaf Info 1``'        Misc. information for module SLAPAF.
*> - '``Symmetry operations``'  The symmetry operations of the point group.
*> - '``Non valence orbitals``' The total number of non valence orbitals per irreducible representation.
*> - '``MkNemo.hDisp``'         The hash matrix for displacements as specified in the mknemo module.
*>
*> @param[in] Label Name of field
*> @param[in] Data  Data to put on runfile
*> @param[in] nData Length of array
************************************************************************
      Subroutine Put_iArray(Label,Data,nData)
      Implicit None
#include "pg_ia_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Label
      Integer       nData
      Integer       Data(nData)
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Character*16 RecLab(nTocIA)
      Integer      RecIdx(nTocIA)
      Integer      RecLen(nTocIA)
      Save         RecLab
      Save         RecIdx
      Save         RecLen
*
      Character*16 CmpLab1
      Character*16 CmpLab2
      Integer      nTmp
      Integer      item
      Integer      iTmp
      Integer      i
*----------------------------------------------------------------------*
* Initialize local variables                                           *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* Do setup if this is the first call.                                  *
*----------------------------------------------------------------------*
      Call ffRun('iArray labels',nTmp,iTmp)
      If(nTmp.eq.0) Then
         Do i=1,nTocIA
            RecLab(i)=' '
            RecIdx(i)=sNotUsed
            RecLen(i)=0
         End Do
*
*        Observe that label is at most 16 characters!
*
*                     1234567890123456
         RecLab(  1)='Center Index    '
         RecLab(  2)='nAsh            '
         RecLab(  3)='nBas            '
         RecLab(  4)='nDel            '
         RecLab(  5)='nFro            '
         RecLab(  6)='nIsh            '
         RecLab(  7)='nIsh beta       '
         RecLab(  8)='nOrb            '
         RecLab(  9)='Orbital Type    '
         RecLab( 10)='Slapaf Info 1   '
         RecLab( 11)='Symmetry operati' !ons
         RecLab( 12)='nIsh_ab         '
         RecLab( 13)='nStab           '
         RecLab( 14)='Quad_c          '
         RecLab( 15)='Quad_i          '
         RecLab( 16)='RFcInfo         '
         RecLab( 17)='RFiInfo         '
         RecLab( 18)='RFlInfo         '
         RecLab( 19)='SCFInfoI        '
         RecLab( 20)='SewCInfo        '
         RecLab( 21)='SewIInfo        '
         RecLab( 22)='SewLInfo        '
         RecLab( 23)='SCFInfoI_ab     '
         RecLab( 24)='nExp            '
         RecLab( 25)='nBasis          '
         RecLab( 26)='ipCff           '
         RecLab( 27)='ipExp           '
         RecLab( 28)='IndS            '
         RecLab( 29)='ip_Occ          '
         RecLab( 30)='ipAkl           '
         RecLab( 31)='ipBk            '
         RecLab( 32)='nOpt            '
         RecLab( 33)='Prjct           '
         RecLab( 34)='Transf          '
         RecLab( 35)='iCoSet          '
         RecLab( 36)='LP_A            '
         RecLab( 37)='NumCho          ' ! Number of Cholesky vectors.
         RecLab( 38)='nFroPT          ' ! Number of Frozen for PT
         RecLab( 39)='nDelPT          ' ! Number of Deleted for PT
         RecLab( 40)='BasType         '
         RecLab( 41)='Spread of Coord.'
         RecLab( 42)='Unit Cell Atoms '
         RecLab( 43)='iSOShl          '
         RecLab( 44)='Non valence orbi' !tals
         RecLab( 45)='LoProp nInts    '
         RecLab( 46)='LoProp iSyLbl   '
         RecLab( 47)='nDel_go         '
         RecLab( 48)='nBas_Prim       '
         RecLab( 49)='IsMM            '
         RecLab( 50)='Atom -> Basis   '
         RecLab( 51)='nBasis_Cntrct   '
         RecLab( 52)='ipCff_Cntrct    '
         RecLab( 53)='ipCff_Prim      '
         RecLab( 54)='SCF nOcc        '
         RecLab( 55)='SCF nOcc_ab     '
         RecLab( 56)='ipFockOp        '
         RecLab( 57)='IrrCmp          '
         RecLab( 58)='iAOtSO          '
         RecLab( 59)='iSOInf          '
         RecLab( 60)='FragShell       '
         RecLab( 61)='AuxShell        '
         RecLab( 62)='nVec_RI         '
         RecLab( 63)='MkNemo.hDisp    '
         RecLab( 64)='Index ZMAT      '
         RecLab( 65)='NAT ZMAT        '
         RecLab( 66)='                ' !Free slot
         RecLab( 67)='nDisp           '
         RecLab( 68)='DegDisp         '
         RecLab( 69)='LBList          '
         RecLab( 71)='Ctr Index Prim  '
         RecLab( 72)='MLTP_SINGLE     '
         RecLab( 73)='JBNUM_SINGLE    '
         RecLab( 74)='LROOT_SINGLE    '
         RecLab( 75)='GeoInfo         '
         RecLab( 76)='Cholesky BkmDim '
         RecLab( 77)='Cholesky BkmVec '
         RecLab( 78)='Atom Types      '
         RecLab( 79)='LA Def          '
         RecLab( 80)='Basis IDs       '
         RecLab( 81)='Desym Basis IDs '
         RecLab( 82)='primitive ids   '
         RecLab( 83)='Root Mapping    '
         RecLab( 84)='Fermion IDs     '
         RecLab( 85)='IsMM Atoms      '
         RecLab( 86)='Un_cen Charge   '
*                     1234567890123456

* Do not go beyond 128 without changing the length of RecLab in include
* file too!
         Call cWrRun('iArray labels',RecLab,16*nTocIA)
         Call iWrRun('iArray indices',RecIdx,nTocIA)
         Call iWrRun('iArray lengths',RecLen,nTocIA)
      Else
         Call cRdRun('iArray labels',RecLab,16*nTocIA)
         Call iRdRun('iArray indices',RecIdx,nTocIA)
         Call iRdRun('iArray lengths',RecLen,nTocIA)
      End If
*----------------------------------------------------------------------*
* Locate item                                                          *
*----------------------------------------------------------------------*
      item=-1
      CmpLab1=Label
      Call UpCase(CmpLab1)
      Do i=1,nTocIA
         CmpLab2=RecLab(i)
         Call UpCase(CmpLab2)
         If(CmpLab1.eq.CmpLab2) item=i
      End Do
*
* Do we create a new temporary field?
*
      If(item.eq.-1) Then
         Do i=1,nTocIA
            If(RecLab(i).eq.' ') item=i
         End Do
         If(item.ne.-1) Then
            RecLab(item)=Label
            RecIdx(item)=sSpecialField
            Call cWrRun('iArray labels',RecLab,16*nTocIA)
            Call iWrRun('iArray indices',RecIdx,nTocIA)
         End If
      End If
*
* Is this a temporary field?
*
      If(item.ne.-1) Then
         If(Recidx(item).eq.sSpecialField) Then
            Write(6,*) '***'
            Write(6,*) '*** Warning, writing temporary iArray field'
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
         Call SysAbendMsg('put_iArray','Could not locate',Label)
      End If
      Call iWrRun(RecLab(item),Data,nData)
      If(RecIdx(item).eq.0) Then
         RecIdx(item)=sRegularField
         Call iWrRun('iArray indices',RecIdx,nTocIA)
      End If
      If(RecLen(item).ne.nData) Then
         RecLen(item)=nData
         Call iWrRun('iArray lengths',RecLen,nTocIA)
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
