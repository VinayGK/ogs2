# cmake-lint: disable=C0103,R0912,R0915
function(OgsTest)
    if(NOT OGS_BUILD_CLI OR NOT OGS_BUILD_TESTING)
        return()
    endif()

    set(options DISABLED)
    set(oneValueArgs PROJECTFILE RUNTIME)
    set(multiValueArgs WRAPPER PROPERTIES LABELS)
    cmake_parse_arguments(
        OgsTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
    )

    get_filename_component(OgsTest_DIR "${OgsTest_PROJECTFILE}" DIRECTORY)
    get_filename_component(OgsTest_NAME "${OgsTest_PROJECTFILE}" NAME)
    get_filename_component(OgsTest_NAME_WE "${OgsTest_PROJECTFILE}" NAME_WE)

    if(OgsTest_UNPARSED_ARGUMENTS)
        message(
            FATAL_ERROR
                "Unparsed argument(s) '${OgsTest_UNPARSED_ARGUMENTS}' to OgsTest call."
        )
    endif()

    if(NOT DEFINED OgsTest_RUNTIME)
        set(OgsTest_RUNTIME 1)
    elseif(OgsTest_RUNTIME GREATER 750)
        # Set a timeout on jobs larger than the default ctest timeout of 1500
        # (s). The allowed runtime is twice as long as the given RUNTIME
        # parameter.
        math(EXPR timeout "${OgsTest_RUNTIME} * 2")
        set(timeout TIMEOUT ${timeout})
    endif()

    if(DEFINED OGS_CTEST_MAX_RUNTIME)
        if(${OgsTest_RUNTIME} GREATER ${OGS_CTEST_MAX_RUNTIME})
            return()
        endif()
    endif()
    if(${OgsTest_RUNTIME} GREATER ${ogs.ctest.large_runtime})
        string(APPEND OgsTest_NAME_WE "-LARGE")
    endif()

    if(DEFINED OGS_EXCLUDE_CTESTS)
        foreach(regexp ${OGS_EXCLUDE_CTESTS})
            if("${OgsTest_NAME}" MATCHES "${regexp}")
                message(
                    STATUS "Disabled by OGS_EXCLUDE_CTESTS: ${OgsTest_NAME}"
                )
                return()
            endif()
        endforeach()
    endif()

    set(OgsTest_SOURCE_DIR "${Data_SOURCE_DIR}/${OgsTest_DIR}")
    set(TEST_NAME "ogs-${OgsTest_DIR}/${OgsTest_NAME_WE}")
    set(_is_petsc_np1_source FALSE)

    # Add wrapper postfix (-mpi for mpirun).
    if(OgsTest_WRAPPER)
        string(REGEX MATCH "^[^ ]+" WRAPPER ${OgsTest_WRAPPER})
        if(WRAPPER STREQUAL "mpirun")
            set(TEST_NAME "${TEST_NAME}-mpi")
            list(FIND OgsTest_WRAPPER "-np" _np_arg_pos)
            if(_np_arg_pos GREATER -1)
                math(EXPR _np_value_pos "${_np_arg_pos} + 1")
                list(LENGTH OgsTest_WRAPPER _wrapper_len)
                if(_np_value_pos LESS _wrapper_len)
                    list(GET OgsTest_WRAPPER ${_np_value_pos} _np_value)
                    if(OGS_ENABLE_PETSC_NP4_VARIANTS AND _np_value STREQUAL "1")
                        set(_is_petsc_np1_source TRUE)
                    endif()
                endif()
            endif()
            list(APPEND OgsTest_WRAPPER --bind-to none)
        endif()
    endif()

    set(_exe_args -r ${OgsTest_SOURCE_DIR}
                  ${OgsTest_SOURCE_DIR}/${OgsTest_NAME}
    )
    set(_enable_petsc_np4_pair FALSE)
    if(_is_petsc_np1_source)
        set(_enable_petsc_np4_pair TRUE)
        set(_petsc_np4_project_file "${OgsTest_SOURCE_DIR}/${OgsTest_NAME}")
        if(EXISTS "${_petsc_np4_project_file}")
            _ogs_project_supports_partitions("${_petsc_np4_project_file}" 4
                                             _supports_np4_partitions
                                             _supports_np4_reason
            )
            if(NOT _supports_np4_partitions)
                set(_enable_petsc_np4_pair FALSE)
                message(
                    STATUS
                        "Skipping PETSc np4 pairing for ${TEST_NAME}: ${_supports_np4_reason}."
                )
            endif()
        endif()
    endif()

    current_dir_as_list(ProcessLib labels)
    if(${AddTest_LABELS})
        list(APPEND labels ${AddTest_LABELS})
    else()
        list(APPEND labels default)
    endif()

    if(${OgsTest_RUNTIME} LESS_EQUAL ${ogs.ctest.large_runtime})
        list(APPEND labels small)
    else()
        list(APPEND labels large)
    endif()
    if(_enable_petsc_np4_pair)
        list(APPEND labels petsc_np1_source)
    endif()

    set(_has_omp_variant FALSE)
    list(JOIN OGS_OPENMP_PARALLEL_ASM_PROCESSES ";|;" match_parallel_asm_processes)
    # OpenMP tests for specific processes only. TODO (CL) Once all processes can
    # be assembled OpenMP parallel, the condition should be removed.
    if(";${labels};" MATCHES ";${match_parallel_asm_processes};")
        set(_has_omp_variant TRUE)
    endif()

    set(_add_non_omp_variant TRUE)
    if(NOT OGS_ENABLE_NON_OMP_TEST_VARIANTS AND _has_omp_variant)
        set(_add_non_omp_variant FALSE)
    endif()

    if(_add_non_omp_variant)
        _ogs_add_test(${TEST_NAME})
        if(_enable_petsc_np4_pair)
            _ogs_add_np4_test_variant(${TEST_NAME})
        endif()
    endif()

    if(_has_omp_variant)
        _ogs_add_test(${TEST_NAME}-omp)
        _set_omp_test_properties_for(${TEST_NAME}-omp)
        if(_enable_petsc_np4_pair)
            _ogs_add_np4_test_variant(${TEST_NAME}-omp)
        endif()
    endif()
endfunction()

# Adds a ctest and sets properties
macro(_ogs_add_test TEST_NAME)
    # TEST_NAME is unique, shortened hash added to the working directory of
    # the test to prevent race conditions.
    set(_unique_string "${TEST_NAME}")
    string(SHA1 _unique_hash "${_unique_string}")
    string(SUBSTRING "${_unique_hash}" 0 8 _short_hash)

    # Create unique directory name
    set(OgsTest_BINARY_DIR "${Data_BINARY_DIR}/${OgsTest_DIR}_${_short_hash}")
    file(MAKE_DIRECTORY ${OgsTest_BINARY_DIR})
    file(TO_NATIVE_PATH "${OgsTest_BINARY_DIR}" OgsTest_BINARY_DIR_NATIVE)
    string(REPLACE "/" "_" TEST_NAME_UNDERSCORE ${TEST_NAME})

    isTestCommandExpectedToSucceed(${TEST_NAME} ${OgsTest_PROPERTIES})
    message(
        DEBUG
        "Is test '${TEST_NAME}' expected to succeed? → ${TEST_COMMAND_IS_EXPECTED_TO_SUCCEED}"
    )

    set(_ogs_exe $<TARGET_FILE:ogs>)
    if(OGS_BUILD_WHEEL)
        # When testing the installed wheel assume executable is in PATH
        # from venv.
        set(_ogs_exe ogs)
    endif()

    add_test(
        NAME ${TEST_NAME}
        COMMAND
            ${CMAKE_COMMAND} -DEXECUTABLE=${_ogs_exe}
            "-DEXECUTABLE_ARGS=${_exe_args}"
            "-DWRAPPER_COMMAND=${OgsTest_WRAPPER}"
            -DWORKING_DIRECTORY=${OgsTest_BINARY_DIR}
            "-DLOG_FILE_BASENAME=${TEST_NAME_UNDERSCORE}.txt"
            "-DLOG_ROOT=${PROJECT_BINARY_DIR}/logs"
            "-DTEST_COMMAND_IS_EXPECTED_TO_SUCCEED=${TEST_COMMAND_IS_EXPECTED_TO_SUCCEED}"
            -P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
    )

    set_tests_properties(
        ${TEST_NAME}
        PROPERTIES ${OgsTest_PROPERTIES}
                   ENVIRONMENT
                   VTKDIFF_EXE=$<TARGET_FILE:vtkdiff>
                   COST
                   ${OgsTest_RUNTIME}
                   DISABLED
                   ${OgsTest_DISABLED}
                   LABELS
                   "${labels}"
                   ${timeout}
    )
endmacro()

macro(_ogs_add_np4_test_variant BASE_TEST_NAME)
    set(_wrapper_backup ${OgsTest_WRAPPER})
    set(_labels_backup ${labels})
    set(_exe_args_backup ${_exe_args})

    set(_np4_wrapper "")
    set(_expect_np_value FALSE)
    foreach(_arg ${OgsTest_WRAPPER})
        if(_expect_np_value)
            list(APPEND _np4_wrapper 4)
            set(_expect_np_value FALSE)
            continue()
        endif()
        list(APPEND _np4_wrapper ${_arg})
        if(_arg STREQUAL "-np")
            set(_expect_np_value TRUE)
        endif()
    endforeach()

    set(_register_np4_variant TRUE)
    set(_np4_setup_test "")
    set(_np4_project_file "${OgsTest_SOURCE_DIR}/${OgsTest_NAME}")
    set(_np4_cfg_ready FALSE)

    if(EXISTS "${_np4_project_file}")
        _ogs_project_supports_partitions("${_np4_project_file}" 4
                                         _supports_np4_partitions
                                         _supports_np4_reason
        )
        if(NOT _supports_np4_partitions)
            set(_register_np4_variant FALSE)
        endif()
        _ogs_project_has_cfg4_for_all_meshes("${_np4_project_file}" 4
                                             _np4_cfg_ready
        )
    endif()

    if(OGS_PETSC_NP4_VARIANTS_REQUIRE_CFG4_READY)
        if(NOT EXISTS "${_np4_project_file}")
            set(_register_np4_variant FALSE)
        elseif(NOT _np4_cfg_ready AND (NOT OGS_PETSC_NP4_AUTO_PREPARE_CFG4
                                       OR NOT TARGET partmesh)
        )
            set(_register_np4_variant FALSE)
        endif()
    endif()

    if(NOT _register_np4_variant)
        if(DEFINED _supports_np4_reason AND NOT _supports_np4_reason STREQUAL "")
            message(
                STATUS
                    "Skipping PETSc np4 variant for ${BASE_TEST_NAME}: ${_supports_np4_reason}."
            )
        else()
            message(
                STATUS
                    "Skipping PETSc np4 variant for ${BASE_TEST_NAME}: cfg4 meshes are not ready and auto-preparation is unavailable."
            )
        endif()
    endif()

    set(OgsTest_WRAPPER ${_np4_wrapper})
    set(_np4_labels ${labels})
    list(REMOVE_ITEM _np4_labels petsc_np1_source)
    list(APPEND _np4_labels petsc_np4_variant)
    set(labels ${_np4_labels})

    if(_register_np4_variant)
        if(EXISTS "${_np4_project_file}")
            set(_np4_mesh_token "${BASE_TEST_NAME}-np4-meshes")
            string(SHA1 _np4_mesh_hash "${_np4_mesh_token}")
            string(SUBSTRING "${_np4_mesh_hash}" 0 8 _np4_mesh_short_hash)
            set(_np4_mesh_dir
                "${Data_BINARY_DIR}/${OgsTest_DIR}_${_np4_mesh_short_hash}_np4_meshes"
            )
            list(APPEND _exe_args -m ${_np4_mesh_dir})
            _ogs_add_np4_mesh_setup_test(
                _np4_setup_test "${BASE_TEST_NAME}-np4" "${_np4_project_file}"
                "${_np4_mesh_dir}" "${OgsTest_DISABLED}"
            )
        endif()

        _ogs_add_test(${BASE_TEST_NAME}-np4)
        if(NOT _np4_setup_test STREQUAL "")
            _ogs_append_test_dependency("${BASE_TEST_NAME}-np4"
                                        "${_np4_setup_test}"
            )
        endif()

        if("${BASE_TEST_NAME}" MATCHES "-omp$")
            _set_omp_test_properties_for(${BASE_TEST_NAME}-np4)
        endif()
    endif()

    set(OgsTest_WRAPPER ${_wrapper_backup})
    set(_exe_args ${_exe_args_backup})
    set(labels ${_labels_backup})
endmacro()

# macro(_set_omp_test_properties) defined in AddTest.cmake
