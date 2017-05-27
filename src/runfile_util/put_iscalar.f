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
* Copyright (C) Per-Olof Widmark                                       *
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
*
* <DOC>
*   <Name>Put\_iScalar</Name>
*   <Syntax>Call Put\_iScalar(Label,Data)</Syntax>
*   <Arguments>
*     \Argument{Label}{Name of field}{Character*(*)}{in}
*     \Argument{Data}{Data to put on runfile}{Integer}{in}
*   </Arguments>
*   <Purpose>To add/update scalar data in runfile.</Purpose>
*   <Dependencies></Dependencies>
*   <Author>Per-Olof Widmark</Author>
*   <Modified_by></Modified_by>
*   <Side_Effects></Side_Effects>
*   <Description>
*     This routine is used to put scalar data of type
*     Integer into the runfile. The data items are
*     identified by the label. Below is a list of the
*     data items that are recognized. The labels are
*     case insensitive and significant to 16 characters.
*
*     For development purposes you can use an unsupported
*     label. Whenever such a field is accessed a warning
*     message is printed in the output, to remind the
*     developer to update this routine.
*
*     List of known labels:
*     \begin{itemize}
*     \item `Multiplicity' is the spin multiplicity of
*           the last SCf or RASSCF calculation.
*     \item `nMEP' Number of points on the minimum energy path.
*     \item `No of Internal coordinates' The number of internal
*                coordinates for the molecule that is allowed
*                within the given point group.
*     \item `nSym' is the number of irreducible representations
*           of the molecule.
*     \item `PCM info length' Length of the block containing
*                 misc. info for the PCM model.
*     \item `Relax CASSCF root' signals which root to perform
*           geometry optimization for in a state average CASSCF
*           geometry optimization.
*     \item `SA ready' signals that SA wavefunction is ready for
*           gradient calculations.
*     \item `System BitSwitch' is a bit switch controlling
*           various functions. Will be replaced!
*     \item `Unique atoms'
*     \item `nActel` is the number of active electrons in CASSCF
*           calculation
*     \item `MkNemo.nMole' is the number of molecules as specified
*           in the mknemo module.
*     \item `nLambda' is the number of constraints in the PCO
*     \item `DNG' force numerical gradients
*     \item `HessIter' Last iteration where the analytical Hessian
*           was computed
*     \item `CHCCLarge` segmentation of VOs in CHCC
*     \item `Seed` is the seed number for random number generator
*           used in surface hoping.
*     \end{itemize}
*   </Description>
* </DOC>
*
************************************************************************
      Subroutine Put_iScalar(Label,Data)
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
      Call ffRun('iScalar labels',nData,iTmp)
      If(nData.eq.0) Then
         Do i=1,nTocIS
            RecLab(i)=' '
            RecVal(i)=0
            RecIdx(i)=sNotUsed
         End Do
*
*        Observe that label is at most 16 characters!
*
*                     1234567890123456
         RecLab(  1)='Multiplicity    '
         RecLab(  2)='nMEP            '
         RecLab(  3)='No of Internal c' !oordinates
         RecLab(  4)='nSym            '
         RecLab(  5)='PCM info length '
         RecLab(  6)='Relax CASSCF roo' !t
         RecLab(  7)='System BitSwitch'
         RecLab(  8)='Unique atoms    '
         RecLab(  9)='LP_nCenter      '
         RecLab( 10)='ChoIni          '
         RecLab( 11)='Unit Cell NAtoms'
         RecLab( 12)='Cholesky Reorder'
         RecLab( 13)='ChoVec Address  '
         RecLab( 14)='SA ready        '
         RecLab( 15)='NumGradRoot     '
         RecLab( 16)='Number of roots '
         RecLab( 17)='LoProp Restart  '
         RecLab( 18)='MpProp nOcOb    '
         RecLab( 19)='Highest Mltpl   '
         RecLab( 20)='nActel          '
         RecLab( 21)='Run_Mode        '
         RecLab( 22)='Grad ready      '
         RecLab( 23)='ISPIN           '
         RecLab( 24)='SCF mode        '
         RecLab( 25)='MkNemo.nMole    '
         RecLab( 26)='N ZMAT          '
         RecLab( 27)='Bfn Atoms       '
         RecLab( 28)='FMM             '
         RecLab( 29)='Pseudo atoms    '
         RecLab( 30)='nChDisp         '
         RecLab( 31)='iOff_Iter       '
         RecLab( 32)='Columbus        '
         RecLab( 33)='ColGradMode     '
         RecLab( 34)='IRC             '
         RecLab( 35)='MaxHops         '
         RecLab( 36)='nRasHole        '
         RecLab( 37)='nRasElec        '
         RecLab( 38)='Rotational Symme' !try Number
         RecLab( 39)='Saddle Iter     '
         RecLab( 40)='iMass           '
         RecLab( 41)='mp2prpt         ' ! True(=1) if mbpt2 was run with prpt
         RecLab( 42)='NJOB_SINGLE     '
         RecLab( 43)='MXJOB_SINGLE    '
         RecLab( 44)='NSS_SINGLE      '
         RecLab( 45)='NSTATE_SINGLE   '
         RecLab( 46)='LDF Status      ' ! Initialized or not
         RecLab( 47)='DF Mode         ' ! Local (1) or non-local (0) DF
         RecLab( 48)='agrad           ' ! Forces analytical gradients
         RecLab( 49)='LDF Constraint  ' ! Constraint type for LDF
         RecLab( 50)='OptimType       ' ! Optimization type in hyper
         RecLab( 51)='LSYM            ' ! symmetry of the CAS root(s)
         RecLab( 52)='RF CASSCF root  '
         RecLab( 53)='RF0CASSCF root  '
         RecLab( 54)='nCoordFiles     ' ! number of xyz-files in gateway
         RecLab( 55)='nLambda         '
         RecLab( 56)='DNG             '
         RecLab( 57)='HessIter        '
c         RecLab( 58)='GEO_nConnect    '
         RecLab( 58)='CHCCLarge       ' ! Segmentation of VOs in CHCC
         RecLab( 59)='TS Search       '
         RecLab( 60)='Number of Hops  '
         RecLab( 61)='hopped          '
         RecLab( 62)='Invert constrain' !ts
         RecLab( 63)='Keep old gradien' !t
         RecLab( 64)='embpot          ' ! Flag whether an embedding potential is present
         RecLab( 65)='nPrim           '
         RecLab( 66)='Seed            '
         RecLab( 67)='Track Done      '
         RecLab( 68)='MaxHopsTully    '
         RecLab( 69)='EFP             ' ! Flag Effective fragment potentials
         RecLab( 70)='nEFP_fragments  '
         RecLab( 71)='Coor_Type       ' ! EFP fragment coordinate format
         RecLab( 72)='nEFP_Coor       ' ! Associated number of coordinates per fragment
*                     1234567890123456
*
*        Note, when the counter here exceeds 128 update this line
*        and the nTocIS parameter in pg_is_info.fh!
*
         Call cWrRun('iScalar labels',RecLab,16*nTocIS)
         Call iWrRun('iScalar values',RecVal,nTocIS)
         Call iWrRun('iScalar indices',RecIdx,nTocIS)
      Else
         Call cRdRun('iScalar labels',RecLab,16*nTocIS)
         Call iRdRun('iScalar values',RecVal,nTocIS)
         Call iRdRun('iScalar indices',RecIdx,nTocIS)
      End If
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
*
* Do we create a new temporary field?
*
      If(item.eq.-1) Then
         Do i=1,nTocIS
            If(RecLab(i).eq.' ') item=i
         End Do
         If(item.ne.-1) Then
            RecLab(item)=Label
            RecIdx(item)=sSpecialField
            Call cWrRun('iScalar labels',RecLab,16*nTocIS)
            Call iWrRun('iScalar indices',RecIdx,nTocIS)
         End If
      End If
*
* Is this a temporary field?
*
      If(item.ne.-1) Then
         If(Recidx(item).eq.sSpecialField) Then
            Write(6,*) '***'
            Write(6,*) '*** Warning, writing temporary iScalar field'
            Write(6,*) '***   Field: ',Label
            Write(6,*) '***'
         End If
      End If
*----------------------------------------------------------------------*
* Write data to disk.                                                  *
*----------------------------------------------------------------------*
      If(item.eq.-1) Then
         Call SysAbendMsg('put_iScalar','Could not locate',Label)
      End If
      RecVal(item)=Data
      Call iWrRun('iScalar values',RecVal,nTocIS)
      If(RecIdx(item).eq.0) Then
         RecIdx(item)=sRegularField
         Call iWrRun('iScalar indices',RecIdx,nTocIS)
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Do i=1,num_IS_init
         If(iLbl_IS_inmem(i).eq.CmpLab1) Then
             i_IS_inmem(i)=Data
             IS_init(i)=1
             return
         End If
      End Do

      Return
      End
