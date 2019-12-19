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
      Subroutine DrvEFP(First)
#ifdef _EFP_
      use EFP_Module
      use EFP
      use iso_c_binding, only: c_int, c_char, c_ptr, c_funptr,
     &                         c_size_t
      Implicit Real*8 (a-h,o-z)
      External Molcas_ELECTRON_DENSITY_FIELD_FN
      Character(len=180) :: CurrDir, MolDir
      Logical First
      Type(c_ptr) :: cptr1
      integer(c_int) :: irc
      Character(kind=c_char):: Name*180, PATH*180
      Integer(c_int) :: Molcas_ELECTRON_DENSITY_FIELD_FN
      Integer(c_size_t) :: frag_idx
      Integer(c_size_t), Target :: n_atoms
*#define _DEBUG_
#ifdef _DEBUG_
      type(efp_energy), Target :: Energy
      Integer(c_int) :: do_gradient
#endif
*
*     Initate the EFP object,
*
      If (First) Then
*
*        Create a new EFP object
*
         EFP_Instance=EFP_Create()
*
*        Now add the potentials
*
#ifdef _DEBUG_
         Write (6,*) 'Initiation of EFP'
         Write (6,*) 'nEFP_fragments=',nEFP_fragments
#endif
         iFrag=0
         Do i = 1, nEFP_fragments
*
*           Find a unique fragment potential.
*           Procastinate.
*
            Do j = i+1, nEFP_fragments
*
*              Branch out if there is an entry later with the same
*              potential
               If (FRAG_TYPE(i).eq.FRAG_TYPE(j)) Go To 999
            End Do
*
*           At this point we have a unique potential label in FRAG_TYPE(i)
*
            If (Index(FRAG_TYPE(i),'_l').ne.0) Then
*
*              .efg file found in the library directory
*
               iLast= Index(FRAG_TYPE(i),'_l')
               Name=FRAG_TYPE(i)(1:iLast-1)
               Call GetEnvf('MOLCAS',MolDir)
               iLast= Index(MolDir,' ')
               If (MolDir(iLast-1:iLast-1).ne.'/')Then
                  MolDir(iLast:iLast)='/'
                  iLast=iLast+1
               End If
               Path=MolDir(1:iLast-1)//'External/efp/fraglib/'
            Else
*
*              .efg file found in the $CurrDir directory
*
               iLast= Index(FRAG_TYPE(i),' ')
               Name=FRAG_TYPE(i)(1:iLast-1)
               Call GetEnvf('CurrDir',CurrDir)
               iLast= Index(CurrDir,' ')
               Path=CurrDir(1:iLast-1)//'/'
            End If
*
            iLast= Index(Path,' ')
            Path=Path(1:iLast-1)//Name
            iLast= Index(Path,' ')
            Path=Path(1:iLast-1)//'.efp'//CHAR(0)
*
            irc = EFP_ADD_POTENTIAL(EFP_Instance,Path)
            If (irc.ne.0) Then
               Write (6,*) 'EFP potential file error.'
               Write (6,*) Path
               Write (6,*) 'Return code:',irc
               Call Abend()
            End If
*
*           Loop over all fragment again an add those that are of the
*           current type
*
            iLast= Index(FRAG_TYPE(i),' ')
            Name=FRAG_TYPE(i)(1:iLast-1)
            iLast= Index(NAME,' ')
            NAME=NAME(1:iLast-1)//CHAR(0)
            Do j = 1, i
               If (FRAG_TYPE(j).ne.FRAG_TYPE(i)) Cycle
               irc  = efp_add_fragment(EFP_Instance,Name)
               If (irc.ne.0) Then
                  Write (6,*) 'EFP_ADD_FRAGMET error.'
                  Write (6,*) 'Return code:',irc
                  Call Abend()
               End If
               cptr1=c_Loc(EFP_COORS(1,j))
               irc  = efp_set_frag_Coordinates(EFP_Instance,iFrag,
     &                                         Coor_type,cptr1)
               If (irc.ne.0) Then
                  Write (6,*) 'EFP_SET_FRAG_COORDINATES error.'
                  Write (6,*) 'Return code:',irc
                  Call Abend()
               End If
               iFrag=iFrag+1
            End Do
*
*
 999        Continue
         End Do
*
         irc=EFP_PREPARE(EFP_Instance)
         If (irc.ne.0) Then
            Write (6,*) 'EFP_PREPARE error.'
            Write (6,*) 'Return code:',irc
            Call Abend()
         End If
*
#ifdef _DEBUG_
         do_gradient=0
         irc=EFP_COMPUTE(EFP_Instance,do_gradient)
         If (irc.ne.0) Then
            Write (6,*) 'EFP_COMPUTE error.'
            Write (6,*) 'Return code:',irc
            Call Abend()
         End If
*
         irc=EFP_GET_ENERGY(EFP_Instance,c_LOC(Energy))
         Write (6,*) Energy%Total
#endif
*
         irc=EFP_SET_ELECTRON_DENSITY_FIELD_FN(EFP_Instance,
     &        c_funloc(Molcas_ELECTRON_DENSITY_FIELD_FN))
*
      End If
*
*     Add EFP charges to the nuclear repulsion term
*
      Do frag_idx = 1, nEFP_fragments
*
*        Pick up the number of atoms in the fragment
*
         irc=EFP_GET_FRAG_ATOM_COUNT(EFP_Instance,frag_idx,
     &                               c_loc(n_atoms))
*
*        Pick up the fragment coordinates and charges
*
*         irc=EFP_GET_FRAG_ATOMS(EFP_Instance,...
*
      End Do
#else
*     Dummy routine
      Logical First
      If (First .or. .Not.First) Return
#endif
      Return
      End
