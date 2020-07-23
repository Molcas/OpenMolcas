foreach( _config neci dneci mneci)

    neci_add_test(
        TARGET test_${_config}_loop
        SOURCES test_loop.F90 test_loop_testcases.F90 test_loop_helpers.F90
        LIBS lib${_config} fruit
        LINKER_LANGUAGE Fortran
        META_TARGET test_${_config}
        DEFINITIONS
            "${${PROJECT_NAME}_${_config}_DEFINITIONS}"
    )

endforeach()


