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

      SUBROUTINE RasScf_Mcontrol(id_call)

      Use Para_Info, Only: MyRank
      Implicit Real*8 (a-h,o-z)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
      Logical  timings,DoLock,Deco
      Logical  DoCholesky
      Integer  ALGO,Nscreen
      Real*8   dmpk
*
      Common /CHLCAS / DoCholesky,ALGO
      COMMON /CHOTIME / timings
      Common /CHOLK / DoLocK,Deco,dmpk,Nscreen
*
      Integer id_call
      Character*512 List
      Character*32 Value
*******************************************************

      icount=0

      If (id_call .eq. 1) Then

* --- Label definitions

         write(List,100) 'RASSCF_started_OK:(-:-):',ALGO,timings,dmpK,
     &                        nScreen,MaxIt,ThrE,ThrSX,ThrTE


* --- Initialize the control system

         Call MolcasControlInit(List)
         Return


      Else

* --- Read the molcas control file
*1
         Call MolcasControl('Cho_ALGO',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) ALGO
            write(6,*)'--- Warning: Cho_ALGO changed by user to the ',
     &'value ',ALGO
            icount = icount + 1
         EndIf
*2
         Call MolcasControl('Chotime',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) timings
            write(6,*)'--- Warning: Cholesky timings visualization ',
     &'changed by user to the value ',timings
            icount = icount + 1
         EndIf
*3
         Call MolcasControl('nScreen',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) nScreen
            write(6,*)'--- Warning: Cholesky LK option nSCREEN changed',
     &' by user to the value ',nScreen
            icount = icount + 1
         EndIf
*4
         Call MolcasControl('dmpK',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) dmpK
            write(6,*)'--- Warning: Cholesky LK option DMPK changed by',
     &' user to the value ',dmpK
            icount = icount + 1
         EndIf
*5
         Call MolcasControl('MaxIter',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) MaxIt
            write(6,*)'--- Warning: MaxIt changed by user to the value '
     &,MaxIt
            icount = icount + 1
         EndIf
*6
         Call MolcasControl('ThrE',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) ThrE
            write(6,*)'--- Warning: ThrE changed by user to the value '
     &,ThrE
            icount = icount + 1
         EndIf
*7
         Call MolcasControl('ThrSX',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) ThrSX
            write(6,*)'--- Warning: ThrSX changed by user to the value '
     &,ThrSX
            icount = icount + 1
         EndIf
*8
         Call MolcasControl('ThrTE',Value)
         If (Value(1:4).ne.'    ') Then
            read(Value,*,err=101,end=102) ThrTE
            write(6,*)'--- Warning: ThrTE changed by user to the value '
     &,ThrTE
            icount = icount + 1
         EndIf

      EndIf

      icount0=icount

* --- Get the true updated counter in parallel runs
*
      Call gaIgOP_SCAL(icount,'max')

      If (MyRank.eq.0) Then
         If (icount.gt.icount0) Then
           write(6,*)' Steering will NOT be activated this time because'
           write(6,*)' molcas.control file must be changed on node_0 !!'
           Call gaIgOP_SCAL(icount,'min')
         EndIf
      EndIf

      If (icount.gt.0) Then

* --- Trick to broadcast the values across nodes
*
         If (MyRank.ne.0) Then
            ALGO=0
            MaxIt=0
            nScreen=0
            dmpK=0.0d0
            ThrE=0.0d0
            ThrSX=0.0d0
            ThrTE=0.0d0
         EndIf
         Call gaIgOP_SCAL(ALGO,'+')
         Call gaIgOP_SCAL(nScreen,'+')
         Call gaIgOP_SCAL(MaxIt,'+')
         Call gadgOP_SCAL(dmpK,'+')
         Call gadgOP_SCAL(ThrE,'+')
         Call gadgOP_SCAL(ThrSX,'+')
         Call gadgOP_SCAL(ThrTE,'+')

* --- Update label values (note that "timings" is updated locally)
*
         write(List,100) 'RASSCF_modified_by_user:',ALGO,timings,dmpK,
     &                       nScreen,MaxIt,ThrE,ThrSX,ThrTE

* --- Initialize the control system with the new values
*
         Call MolcasControlInit(List)
         Return

      EndIf

      Return

100   FORMAT(A24,
     &',Cho_ALGO=',I2,
     &',Chotime=',L2,
     &',dmpK=',E11.4,
     &',nScreen=',I4,
     &',MaxIter=',I4,
     &',ThrE=',E11.4,
     &',ThrSX=',E11.4,
     &',ThrTE=',E11.4)

101   write(6,*) 'RasScf_Mcontrol: error in data Input. ( icount= ',
     &           icount,' )'
102   write(6,*) 'RasScf_Mcontrol: reached end of file. ( icount= ',
     &           icount,' )'


      End
