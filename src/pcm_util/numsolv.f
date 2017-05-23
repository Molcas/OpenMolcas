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
      Integer Function NumSolv(Solvent)
      Implicit Real*8 (a-h,o-z)
      Character*30 Solvent, Solvent_
      IdSolv = 0
      Solvent_=Solvent
      Call Upcase(Solvent_)
      If(Solvent_.eq.'WATER')               IdSolv = 1
      If(Solvent_.eq.'ACETONITRILE')        IdSolv = 2
      If(Solvent_.eq.'METHANOL')            IdSolv = 3
      If(Solvent_.eq.'ETHANOL')             IdSolv = 4
      If(Solvent_.eq.'ISOQUINOLINE')        IdSolv = 5
      If(Solvent_.eq.'QUINOLINE')           IdSolv = 6
      If(Solvent_.eq.'CHLOROFORM')          IdSolv = 7
      If(Solvent_.eq.'ETHYLETHER')          IdSolv = 8
      If(Solvent_.eq.'METHYLENECHLORIDE')   IdSolv = 9
      If(Solvent_.eq.'DICHLOROETHANE')      IdSolv =10
      If(Solvent_.eq.'CARBONTETRACHLORIDE') IdSolv =11
      If(Solvent_.eq.'BENZENE')             IdSolv =12
      If(Solvent_.eq.'TOLUENE')             IdSolv =13
      If(Solvent_.eq.'CHLOROBENZENE')       IdSolv =14
      If(Solvent_.eq.'NITROMETHANE')        IdSolv =15
      If(Solvent_.eq.'HEPTANE')             IdSolv =16
      If(Solvent_.eq.'CYCLOHEXANE')         IdSolv =17
      If(Solvent_.eq.'ANILINE')             IdSolv =18
      If(Solvent_.eq.'ACETONE')             IdSolv =19
      If(Solvent_.eq.'TETRAHYDROFURAN')     IdSolv =20
      If(Solvent_.eq.'DIMETHYLSULFOXIDE')   IdSolv =21
      If(Solvent_.eq.'ARGON')               IdSolv =22
      If(Solvent_.eq.'KRYPTON')             IdSolv =23
      If(Solvent_.eq.'XENON')               IdSolv =24

      If(IdSolv.eq.0) then
*       Call ErrTra
        Write (6,*) ' Unrecognized solvent: ',Solvent
        Write (6,10)'WATER', 'ACETONITRILE', 'METHANOL',
     &  'ETHANOL', 'ISOQUINOLINE', 'QUINOLINE', 'CHLOROFORM',
     &  'ETHYLETHER', 'METHYLENECHLORIDE', 'DICHLOROETHANE',
     &  'CARBONTETRACHLORIDE', 'BENZENE', 'TOLUENE',
     &  'CHLOROBENZENE', 'NITROMETHANE' , 'HEPTANE', 'CYCLOHEXANE',
     &  'ANILINE', 'ACETONE', 'TETRAHYDROFURAN', 'DIMETHYLSULFOXIDE',
     &  'ARGON', 'KRYPTON', 'XENON'
  10    format(' Allowed solvents are:',/,24(A30,/))
        Call Abend()
      EndIf
      NumSolv = IdSolv
      Return
      end
