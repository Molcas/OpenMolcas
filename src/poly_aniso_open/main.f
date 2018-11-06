* $ this file belongs to the Molcas repository $
      program main
#ifdef _FPE_TRAP_
      Use, Intrinsic :: IEEE_Exceptions
#endif
      Integer iReturn
      Character(20) Module_Name
      Parameter (Module_Name = 'poly_aniso')
#ifdef _FPE_TRAP_
      Call IEEE_Set_Halting_Mode(IEEE_Usual,.True._4)
#endif

      Call Start(Module_Name)
      Call poly_aniso(iReturn)
      Call QStat(' ')
      Call Finish(iReturn)
      End
