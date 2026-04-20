#
# AddTest
# -------
#
# Creates application test runs. Order of arguments can be arbitrary.
#
# ~~~
# AddTest(
#   NAME <name of the the test>
#   PATH <working directory> # relative to SourceDir/Tests/Data
#   EXECUTABLE <executable target> # optional, defaults to ogs
#   EXECUTABLE_ARGS <arguments>
#   WRAPPER <time|mpirun> # optional
#   WRAPPER_ARGS <arguments> # optional
#   TESTER <diff|vtkdiff|vtkdiff-mesh|gmldiff|memcheck|numdiff> # optional
#   TESTER_ARGS <argument> # optional
#   REQUIREMENTS # optional simple boolean expression which has to be true to
#                  enable the test, e.g.
#                  OGS_USE_PETSC AND (FOO OR BAR)
#   PYTHON_PACKAGES package_x=1.2.3 package_y=0.1.x # optional
#   RUNTIME <in seconds> # optional for optimizing ctest duration
#                          values should be taken from envinf job
#   LABELS <labelA;labelB;...> # optional, defaults to "default"
#   WORKING_DIRECTORY # optional, specify the working directory of the test
#   DISABLED # optional, disables the test
#   PROPERTIES <test properties> # optional
# )
# ~~~
#
# Conditional arguments:
#
# ~~~
#   diff-tester
#     - DIFF_DATA <list of files to diff>
#         the given file is compared to a file with the same name from Tests/Data
#
#   vtkdiff-tester
#     - DIFF_DATA
#         <vtk file a> <vtk file b> <data array a name> <data array b name> <absolute tolerance> <relative tolerance>
#         Can be given multiple times; the given data arrays in the vtk files are
#         compared using the given absolute and relative tolerances.
#       OR
#     - DIFF_DATA
#         GLOB <globbing expression, e.g. xyz*.vtu> <data array a name> <data array b name> <absolute tolerance> <relative tolerance>
#         Searches for all matching files in the working directory (PATH).
#         Matched files are then compared against files with the same name in
#         the benchmark output directory.
#
#   gmldiff-tester
#     - DIFF_DATA
#         <gml file> <absolute tolerance> <relative tolerance>
#         Can be given multiple times; the point coordinates in the gml files are
#         compared using the given absolute and relative tolerances.
# ~~~
# cmake-lint: disable=C0103,R0911,R0912,R0915
function(AddTest)

    # parse arguments
    set(options DISABLED)
    set(oneValueArgs
        EXECUTABLE
        PATH
        NAME
        WRAPPER
        TESTER
        ABSTOL
        RELTOL
        RUNTIME
        DEPENDS
        WORKING_DIRECTORY
    )
    set(multiValueArgs
        EXECUTABLE_ARGS
        DATA
        DIFF_DATA
        WRAPPER_ARGS
        TESTER_ARGS
        REQUIREMENTS
        PYTHON_PACKAGES
        PROPERTIES
        LABELS
    )
    cmake_parse_arguments(
        AddTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
    )

    set(AddTest_SOURCE_PATH "${Data_SOURCE_DIR}/${AddTest_PATH}")
    set(AddTest_BINARY_PATH "${Data_BINARY_DIR}/${AddTest_PATH}")
    set(AddTest_STDOUT_FILE_PATH
        "${AddTest_BINARY_PATH}/${AddTest_NAME}_stdout.txt"
    )

    # set defaults
    if(NOT DEFINED AddTest_EXECUTABLE)
        message(FATAL_ERROR "Test ${AddTest_NAME}: No EXECUTABLE set!")
    endif()
    if(NOT DEFINED AddTest_REQUIREMENTS)
        set(AddTest_REQUIREMENTS TRUE)
    endif()
    if(NOT DEFINED AddTest_RUNTIME)
        set(AddTest_RUNTIME 1)
    elseif(AddTest_RUNTIME GREATER 750)
        # Set a timeout on jobs larger than the default ctest timeout of 1500
        # (s). The allowed runtime is twice as long as the given RUNTIME
        # parameter.
        math(EXPR timeout "${AddTest_RUNTIME} * 2")
        set(timeout TIMEOUT ${timeout})
    endif()
    if(NOT DEFINED AddTest_WORKING_DIRECTORY)
        set(AddTest_WORKING_DIRECTORY ${AddTest_BINARY_PATH})
    endif()

    if("${AddTest_EXECUTABLE}" STREQUAL "ogs")
        set(AddTest_WORKING_DIRECTORY ${AddTest_SOURCE_PATH})
    endif()

    if(DEFINED OGS_CTEST_MAX_RUNTIME)
        if(${AddTest_RUNTIME} GREATER ${OGS_CTEST_MAX_RUNTIME})
            return()
        endif()
    endif()
    if(${AddTest_RUNTIME} GREATER ${ogs.ctest.large_runtime})
        string(PREPEND AddTest_NAME "LARGE_")
    endif()

    if(DEFINED OGS_EXCLUDE_CTESTS)
        foreach(regexp ${OGS_EXCLUDE_CTESTS})
            if("${AddTest_NAME}" MATCHES "${regexp}")
                message(
                    STATUS "Disabled by OGS_EXCLUDE_CTESTS: ${AddTest_NAME}"
                )
                return()
            endif()
        endforeach()
    endif()

    # check requirements, disable if not met

    # When testing the installed wheel assume executable is in PATH
    # from venv.
    if(NOT TARGET ${AddTest_EXECUTABLE} AND NOT OGS_BUILD_WHEEL)
        return()
    endif()
    if(${AddTest_REQUIREMENTS})
        message(DEBUG "Enabling test ${AddTest_NAME}.")
    else()
        return()
    endif()

    set(_is_petsc_np1_source FALSE)
    unset(MPI_PROCESSORS)

    # --- Implement wrappers ---
    if(AddTest_WRAPPER STREQUAL "time")
        if(TIME_TOOL_PATH)
            set(WRAPPER_COMMAND time)
        else()
            set(AddTest_WRAPPER_ARGS "")
        endif()
    elseif(AddTest_WRAPPER STREQUAL "mpirun")
        if(MPIRUN_TOOL_PATH)
            if("${HOSTNAME}" MATCHES "frontend.*")
                list(APPEND AddTest_WRAPPER_ARGS --mca btl_openib_allow_ib 1)
            endif()
            list(APPEND AddTest_WRAPPER_ARGS --bind-to none)
            set(WRAPPER_COMMAND ${MPIRUN_TOOL_PATH})
            if("${AddTest_WRAPPER_ARGS}" MATCHES "-np;([0-9]*)")
                set(MPI_PROCESSORS ${CMAKE_MATCH_1})
                if(OGS_ENABLE_PETSC_NP4_VARIANTS AND MPI_PROCESSORS STREQUAL "1")
                    set(_is_petsc_np1_source TRUE)
                endif()
            endif()
        else()
            message(
                STATUS
                    "ERROR: mpirun was not found but is required for ${AddTest_NAME}!"
            )
            return()
        endif()
    endif()

    # --- Implement testers ---
    # check requirements, disable if not met
    if(AddTest_TESTER STREQUAL "diff" AND NOT DIFF_TOOL_PATH)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "vtkdiff" AND NOT TARGET vtkdiff)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "vtkdiff-mesh" AND NOT TARGET vtkdiff)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "xdmfdiff" AND NOT TARGET xdmfdiff)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "gmldiff" AND NOT ${Python_Interpreter_FOUND})
        return()
    endif()
    if(AddTest_TESTER STREQUAL "memcheck" AND NOT GREP_TOOL_PATH)
        return()
    endif()
    if(AddTest_TESTER STREQUAL "numdiff" AND NOT NUMDIFF_TOOL_PATH)
        return()
    endif()

    if(AddTest_DIFF_DATA)
        string(LENGTH "${AddTest_DIFF_DATA}" DIFF_DATA_LENGTH)
        if(${DIFF_DATA_LENGTH} GREATER 7500)
            message(
                FATAL_ERROR
                    "${AddTest_NAME}: DIFF_DATA to long! Consider using regex-syntax: TODO"
            )
        endif()
    endif()

    if((AddTest_TESTER STREQUAL "diff" OR AddTest_TESTER MATCHES "vtkdiff"
        OR AddTest_TESTER STREQUAL "xdmfdiff") AND NOT AddTest_DIFF_DATA
    )
        message(FATAL_ERROR "AddTest(): ${AddTest_NAME} - no DIFF_DATA given!")
    endif()

    if(AddTest_TESTER STREQUAL "diff")
        set(SELECTED_DIFF_TOOL_PATH ${DIFF_TOOL_PATH})
        set(TESTER_ARGS "-sbB")
    elseif(AddTest_TESTER MATCHES "vtkdiff")
        set(SELECTED_DIFF_TOOL_PATH $<TARGET_FILE:vtkdiff>)
    elseif(AddTest_TESTER STREQUAL "xdmfdiff")
        set(SELECTED_DIFF_TOOL_PATH $<TARGET_FILE:xdmfdiff>)
    elseif(AddTest_TESTER STREQUAL "numdiff")
        set(SELECTED_DIFF_TOOL_PATH ${NUMDIFF_TOOL_PATH})
    endif()

    # -----------
    if(TARGET ${AddTest_EXECUTABLE} AND NOT OGS_BUILD_WHEEL)
        set(AddTest_EXECUTABLE_PARSED $<TARGET_FILE:${AddTest_EXECUTABLE}>)
    else()
        # When testing the installed wheel assume executable is in PATH
        # from venv.
        set(AddTest_EXECUTABLE_PARSED ${AddTest_EXECUTABLE})
    endif()

    # Run the wrapper
    if(DEFINED AddTest_WRAPPER)
        set(AddTest_WRAPPER_STRING "-${AddTest_WRAPPER}")
    endif()
    set(TEST_NAME
        "${AddTest_EXECUTABLE}-${AddTest_NAME}${AddTest_WRAPPER_STRING}"
    )

    # Process placeholders <PATH> <SOURCE_PATH> and <BUILD_PATH>
    string(REPLACE "<PATH>" "${AddTest_PATH}" AddTest_WORKING_DIRECTORY
                   "${AddTest_WORKING_DIRECTORY}"
    )
    string(REPLACE "<PATH>" "${AddTest_PATH}" AddTest_EXECUTABLE_ARGS
                   "${AddTest_EXECUTABLE_ARGS}"
    )

    string(REPLACE "<SOURCE_PATH>" "${Data_SOURCE_DIR}/${AddTest_PATH}"
                   AddTest_WORKING_DIRECTORY "${AddTest_WORKING_DIRECTORY}"
    )
    string(REPLACE "<SOURCE_PATH>" "${Data_SOURCE_DIR}/${AddTest_PATH}"
                   AddTest_EXECUTABLE_ARGS "${AddTest_EXECUTABLE_ARGS}"
    )
    string(REPLACE "<BUILD_PATH>" "${Data_BINARY_DIR}/${AddTest_PATH}"
                   AddTest_WORKING_DIRECTORY "${AddTest_WORKING_DIRECTORY}"
    )
    string(REPLACE "<BUILD_PATH>" "${Data_BINARY_DIR}/${AddTest_PATH}"
                   AddTest_EXECUTABLE_ARGS "${AddTest_EXECUTABLE_ARGS}"
    )

    set(_enable_petsc_np4_pair FALSE)
    if(_is_petsc_np1_source AND "${AddTest_EXECUTABLE}" STREQUAL "ogs")
        set(_enable_petsc_np4_pair TRUE)
        _ogs_find_project_file_from_args(
            _petsc_np4_project_file "${AddTest_SOURCE_PATH}" ${AddTest_EXECUTABLE_ARGS}
        )
        if(NOT _petsc_np4_project_file STREQUAL "")
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
    if(DEFINED AddTest_LABELS)
        list(APPEND labels ${AddTest_LABELS})
    else()
        list(APPEND labels default)
    endif()
    if(${AddTest_RUNTIME} LESS_EQUAL ${ogs.ctest.large_runtime})
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
        _add_test(${TEST_NAME})
        if(_enable_petsc_np4_pair)
            _add_np4_test_variant(${TEST_NAME})
        endif()
    endif()

    if(_has_omp_variant)
        _add_test(${TEST_NAME}-omp)
        _set_omp_test_properties_for(${TEST_NAME}-omp)
        if(_enable_petsc_np4_pair)
            _add_np4_test_variant(${TEST_NAME}-omp)
        endif()
    endif()

    if(NOT AddTest_TESTER OR OGS_COVERAGE)
        return()
    endif()

    # Run the tester
    if(_add_non_omp_variant)
        _add_test_tester(${TEST_NAME})
        if(_enable_petsc_np4_pair AND OGS_PETSC_NP4_INCLUDE_TESTER_VARIANTS)
            _add_test_tester(${TEST_NAME}-np4)
        endif()
    elseif(_has_omp_variant)
        _add_test_tester(${TEST_NAME}-omp)
        if(_enable_petsc_np4_pair AND OGS_PETSC_NP4_INCLUDE_TESTER_VARIANTS)
            _add_test_tester(${TEST_NAME}-omp-np4)
        endif()
    endif()

endfunction()

function(_ogs_resolve_project_file INPUT_PROJECT_FILE OUT_PROJECT_FILE)
    set(_resolved_project_file "${INPUT_PROJECT_FILE}")
    if(EXISTS "${INPUT_PROJECT_FILE}")
        file(READ "${INPUT_PROJECT_FILE}" _project_xml LIMIT 4096)
        if(_project_xml MATCHES "<OpenGeoSysProjectDiff")
            get_filename_component(_project_dir "${INPUT_PROJECT_FILE}" DIRECTORY)
            get_filename_component(_project_stem "${INPUT_PROJECT_FILE}" NAME_WE)
            set(_fallback_project "${_project_dir}/${_project_stem}.prj")
            if(EXISTS "${_fallback_project}")
                set(_resolved_project_file "${_fallback_project}")
            else()
                file(GLOB _project_candidates "${_project_dir}/*.prj")
                list(LENGTH _project_candidates _num_project_candidates)
                if(_num_project_candidates EQUAL 1)
                    list(GET _project_candidates 0 _resolved_project_file)
                endif()
            endif()
        endif()
    endif()
    set(${OUT_PROJECT_FILE} "${_resolved_project_file}" PARENT_SCOPE)
endfunction()

function(_ogs_find_project_file_from_args OUT_PROJECT_FILE BASE_DIR)
    set(_project_file "")
    foreach(_arg ${ARGN})
        if(NOT _arg MATCHES "\\.(prj|xml)$")
            continue()
        endif()

        if(IS_ABSOLUTE "${_arg}")
            set(_candidate "${_arg}")
        else()
            set(_candidate "${BASE_DIR}/${_arg}")
        endif()
        if(NOT EXISTS "${_candidate}")
            continue()
        endif()

        file(READ "${_candidate}" _candidate_xml LIMIT 4096)
        if(NOT _candidate_xml MATCHES "OpenGeoSysProject")
            continue()
        endif()

        set(_project_file "${_candidate}")
        break()
    endforeach()
    set(${OUT_PROJECT_FILE} "${_project_file}" PARENT_SCOPE)
endfunction()

function(_ogs_collect_projectdiff_mesh_replacements PROJECT_FILE OUT_MESHES)
    file(READ "${PROJECT_FILE}" _project_xml)
    string(REGEX MATCHALL
                 "<replace[^>]*sel=\"[^\"]*mesh[^\"\n]*text\\(\\)[^\"]*\"[^>]*>[^<]+</replace>"
                 _replace_matches "${_project_xml}"
    )
    set(_meshes "")
    foreach(_match IN LISTS _replace_matches)
        string(REGEX REPLACE ".*>([^<]+)</replace>.*" "\\1" _mesh "${_match}")
        string(STRIP "${_mesh}" _mesh)
        if(NOT _mesh STREQUAL "")
            list(APPEND _meshes "${_mesh}")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _meshes)
    set(${OUT_MESHES} "${_meshes}" PARENT_SCOPE)
endfunction()

function(_ogs_collect_project_meshes PROJECT_FILE OUT_MESHES)
    set(_meshes "")
    if(EXISTS "${PROJECT_FILE}")
        file(READ "${PROJECT_FILE}" _project_xml LIMIT 4096)
        if(_project_xml MATCHES "<OpenGeoSysProjectDiff")
            _ogs_collect_projectdiff_mesh_replacements("${PROJECT_FILE}" _meshes)
        endif()
    endif()
    if(_meshes)
        set(${OUT_MESHES} "${_meshes}" PARENT_SCOPE)
        return()
    endif()

    _ogs_resolve_project_file("${PROJECT_FILE}" _resolved_project_file)
    file(STRINGS "${_resolved_project_file}" _project_lines)

    set(_inside_meshes_block FALSE)
    set(_before_process_sections TRUE)
    set(_inside_mesh_tag FALSE)
    set(_pending_mesh "")
    set(_resolved_meshes "")
    foreach(_line IN LISTS _project_lines)
        string(STRIP "${_line}" _line)
        if(_line MATCHES "^<meshes[^>]*>")
            set(_inside_meshes_block TRUE)
        elseif(_line MATCHES "^</meshes>")
            set(_inside_meshes_block FALSE)
        elseif(
            _line MATCHES
            "^<(processes|time_loop|media|parameters|process_variables|nonlinear_solvers|linear_solvers|curves|python_script|output|deactivated_subdomains)>"
        )
            set(_before_process_sections FALSE)
        endif()

        if(_inside_mesh_tag)
            if(_line MATCHES "^(.*)</mesh>.*$")
                string(APPEND _pending_mesh "${CMAKE_MATCH_1}")
                set(_inside_mesh_tag FALSE)
                string(STRIP "${_pending_mesh}" _mesh)
                set(_pending_mesh "")
                if(NOT _mesh STREQUAL "")
                    list(APPEND _resolved_meshes "${_mesh}")
                endif()
            else()
                string(APPEND _pending_mesh "${_line}")
            endif()
            continue()
        endif()

        if(NOT _inside_meshes_block AND NOT _before_process_sections)
            continue()
        endif()

        if(_line MATCHES ".*<mesh([ \t][^>]*)?>([^<]+)</mesh>.*")
            set(_mesh "${CMAKE_MATCH_2}")
            string(STRIP "${_mesh}" _mesh)
            if(NOT _mesh STREQUAL "")
                list(APPEND _resolved_meshes "${_mesh}")
            endif()
            continue()
        endif()

        if(_line MATCHES ".*<mesh([ \t][^>]*)?>(.*)$")
            set(_pending_mesh "${CMAKE_MATCH_2}")
            if(_pending_mesh MATCHES "^(.*)</mesh>.*$")
                set(_inside_mesh_tag FALSE)
                set(_mesh "${CMAKE_MATCH_1}")
                string(STRIP "${_mesh}" _mesh)
                set(_pending_mesh "")
                if(NOT _mesh STREQUAL "")
                    list(APPEND _resolved_meshes "${_mesh}")
                endif()
            else()
                set(_inside_mesh_tag TRUE)
            endif()
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _resolved_meshes)
    set(${OUT_MESHES} "${_resolved_meshes}" PARENT_SCOPE)
endfunction()

function(_ogs_resolve_mesh_source PROJECT_DIR MESH_REFERENCE OUT_MESH_SOURCE)
    set(_mesh_source "")
    if(IS_ABSOLUTE "${MESH_REFERENCE}" AND EXISTS "${MESH_REFERENCE}")
        set(_mesh_source "${MESH_REFERENCE}")
    else()
        set(_candidate "${PROJECT_DIR}/${MESH_REFERENCE}")
        if(EXISTS "${_candidate}")
            set(_mesh_source "${_candidate}")
        else()
            get_filename_component(_mesh_ext "${MESH_REFERENCE}" EXT)
            if(_mesh_ext STREQUAL "")
                foreach(_ext .vtu .msh .vtk .vtm .vti)
                    set(_candidate_with_ext "${PROJECT_DIR}/${MESH_REFERENCE}${_ext}")
                    if(EXISTS "${_candidate_with_ext}")
                        set(_mesh_source "${_candidate_with_ext}")
                        break()
                    endif()
                endforeach()
            endif()
        endif()
    endif()

    set(${OUT_MESH_SOURCE} "${_mesh_source}" PARENT_SCOPE)
endfunction()

function(_ogs_mesh_has_min_elements MESH_FILE MIN_ELEMENTS OUT_HAS_MIN_ELEMENTS)
    set(_has_min_elements TRUE)
    if(EXISTS "${MESH_FILE}")
        get_filename_component(_mesh_ext "${MESH_FILE}" EXT)
        string(TOLOWER "${_mesh_ext}" _mesh_ext)

        if(_mesh_ext STREQUAL ".vtu")
            file(READ "${MESH_FILE}" _mesh_header LIMIT 16384)
            if(_mesh_header MATCHES "NumberOfCells=\"([0-9]+)\"")
                if(${CMAKE_MATCH_1} LESS ${MIN_ELEMENTS})
                    set(_has_min_elements FALSE)
                endif()
            endif()
        elseif(_mesh_ext STREQUAL ".msh")
            file(READ "${MESH_FILE}" _mesh_header LIMIT 32768)
            if(_mesh_header MATCHES "\\$Elements[ \t\r\n]+([0-9]+)")
                if(${CMAKE_MATCH_1} LESS ${MIN_ELEMENTS})
                    set(_has_min_elements FALSE)
                endif()
            endif()
        endif()
    endif()

    set(${OUT_HAS_MIN_ELEMENTS} ${_has_min_elements} PARENT_SCOPE)
endfunction()

function(
    _ogs_project_supports_partitions
    PROJECT_FILE
    N_PARTITIONS
    OUT_SUPPORTS_PARTITIONS
    OUT_REASON
)
    set(_supports_partitions TRUE)
    set(_reason "")

    if(NOT EXISTS "${PROJECT_FILE}")
        set(${OUT_SUPPORTS_PARTITIONS} TRUE PARENT_SCOPE)
        set(${OUT_REASON} "" PARENT_SCOPE)
        return()
    endif()

    _ogs_resolve_project_file("${PROJECT_FILE}" _resolved_project_file)
    _ogs_collect_project_meshes("${PROJECT_FILE}" _project_meshes)
    if(NOT _project_meshes)
        set(${OUT_SUPPORTS_PARTITIONS} TRUE PARENT_SCOPE)
        set(${OUT_REASON} "" PARENT_SCOPE)
        return()
    endif()

    list(GET _project_meshes 0 _primary_mesh_reference)
    get_filename_component(_project_dir "${PROJECT_FILE}" DIRECTORY)
    get_filename_component(_resolved_project_dir "${_resolved_project_file}" DIRECTORY)

    _ogs_resolve_mesh_source("${_project_dir}" "${_primary_mesh_reference}"
                             _primary_mesh_source
    )
    if(_primary_mesh_source STREQUAL ""
       AND NOT _resolved_project_dir STREQUAL "${_project_dir}"
    )
        _ogs_resolve_mesh_source("${_resolved_project_dir}" "${_primary_mesh_reference}"
                                 _primary_mesh_source
        )
    endif()

    if(_primary_mesh_source STREQUAL "")
        set(${OUT_SUPPORTS_PARTITIONS} TRUE PARENT_SCOPE)
        set(${OUT_REASON} "" PARENT_SCOPE)
        return()
    endif()

    math(EXPR _minimum_recommended_elements "${N_PARTITIONS} * 2 + 1")
    _ogs_mesh_has_min_elements("${_primary_mesh_source}"
                               "${_minimum_recommended_elements}"
                               _has_min_elements
    )
    if(NOT _has_min_elements)
        set(_supports_partitions FALSE)
        set(
            _reason
            "primary mesh '${_primary_mesh_reference}' has fewer than ${_minimum_recommended_elements} elements"
        )
    endif()

    set(${OUT_SUPPORTS_PARTITIONS} ${_supports_partitions} PARENT_SCOPE)
    set(${OUT_REASON} "${_reason}" PARENT_SCOPE)
endfunction()

function(_ogs_project_has_cfg4_for_all_meshes PROJECT_FILE N_PARTITIONS OUT_HAS_CFG)
    if(NOT EXISTS "${PROJECT_FILE}")
        set(${OUT_HAS_CFG} FALSE PARENT_SCOPE)
        return()
    endif()

    _ogs_resolve_project_file("${PROJECT_FILE}" _resolved_project_file)
    _ogs_collect_project_meshes("${PROJECT_FILE}" _project_meshes)
    if(NOT _project_meshes)
        set(${OUT_HAS_CFG} FALSE PARENT_SCOPE)
        return()
    endif()

    get_filename_component(_project_dir "${PROJECT_FILE}" DIRECTORY)
    get_filename_component(_resolved_project_dir "${_resolved_project_file}" DIRECTORY)
    set(_has_cfg TRUE)
    foreach(_mesh IN LISTS _project_meshes)
        _ogs_resolve_mesh_source("${_project_dir}" "${_mesh}" _mesh_source)
        if(_mesh_source STREQUAL "" AND NOT _resolved_project_dir STREQUAL "${_project_dir}")
            _ogs_resolve_mesh_source("${_resolved_project_dir}" "${_mesh}" _mesh_source)
        endif()
        if(_mesh_source STREQUAL "")
            set(_has_cfg FALSE)
            break()
        endif()
        get_filename_component(_mesh_source_dir "${_mesh_source}" DIRECTORY)
        get_filename_component(_mesh_name "${_mesh_source}" NAME)
        string(REGEX REPLACE "\\.[^.]+$" "" _mesh_base "${_mesh_name}")
        set(_cfg_file
            "${_mesh_source_dir}/${_mesh_base}_partitioned_msh_cfg${N_PARTITIONS}.bin"
        )
        if(NOT EXISTS "${_cfg_file}")
            set(_has_cfg FALSE)
            break()
        endif()
    endforeach()

    set(${OUT_HAS_CFG} ${_has_cfg} PARENT_SCOPE)
endfunction()

function(_ogs_append_test_dependency TEST_NAME DEPENDENCY_TEST)
    get_test_property(${TEST_NAME} DEPENDS _existing_depends)
    if(_existing_depends AND NOT _existing_depends STREQUAL "NOTFOUND")
        set(_depends ${_existing_depends})
    else()
        set(_depends "")
    endif()
    list(APPEND _depends ${DEPENDENCY_TEST})
    list(REMOVE_DUPLICATES _depends)
    set_tests_properties(${TEST_NAME} PROPERTIES DEPENDS "${_depends}")
endfunction()

function(
    _ogs_add_np4_mesh_setup_test
    OUT_SETUP_TEST_NAME
    NP4_TEST_NAME
    PROJECT_FILE
    MESH_OUTPUT_DIR
    DISABLED
)
    set(_partmesh_executable "")
    if(TARGET partmesh)
        set(_partmesh_executable "$<TARGET_FILE:partmesh>")
    endif()

    set(_setup_test_name "${NP4_TEST_NAME}-prepare-meshes")
    get_filename_component(_project_dir "${PROJECT_FILE}" DIRECTORY)

    add_test(
        NAME ${_setup_test_name}
        COMMAND
            ${CMAKE_COMMAND} -DPROJECT_FILE=${PROJECT_FILE}
            -DMESH_OUTPUT_DIR=${MESH_OUTPUT_DIR} -DN_PARTITIONS=4
            "-DPARTMESH_EXECUTABLE=${_partmesh_executable}" -P
            ${PROJECT_SOURCE_DIR}/scripts/cmake/test/PreparePartitionedMesh.cmake
        WORKING_DIRECTORY ${_project_dir}
    )

    set_tests_properties(
        ${_setup_test_name}
        PROPERTIES COST
                   1
                   DISABLED
                   ${DISABLED}
                   LABELS
                   "petsc_np4_setup"
    )

    set(${OUT_SETUP_TEST_NAME} ${_setup_test_name} PARENT_SCOPE)
endfunction()

# Add a ctest and sets properties
macro(_add_test TEST_NAME)

    set(_binary_path ${AddTest_BINARY_PATH})
    if("${TEST_NAME}" MATCHES "-omp")
        set(_binary_path ${_binary_path}-omp)
    endif()
    if("${TEST_NAME}" MATCHES "-np4")
        set(_binary_path ${_binary_path}-np4)
    endif()

    file(TO_NATIVE_PATH "${_binary_path}" AddTest_BINARY_PATH_NATIVE)
    file(MAKE_DIRECTORY ${_binary_path})
    set(_exe_args ${AddTest_EXECUTABLE_ARGS})
    if("${AddTest_EXECUTABLE}" STREQUAL "ogs")
        set(_exe_args -o ${AddTest_BINARY_PATH_NATIVE}
                      ${AddTest_EXECUTABLE_ARGS}
        )
    endif()

    isTestCommandExpectedToSucceed(${TEST_NAME} ${AddTest_PROPERTIES})
    message(
        DEBUG
        "Is test '${TEST_NAME}' expected to succeed? → ${TEST_COMMAND_IS_EXPECTED_TO_SUCCEED}"
    )

    if(OGS_USE_PIP AND DEFINED AddTest_PYTHON_PACKAGES)
        list(APPEND labels additional_python_modules)
        foreach(_package ${AddTest_PYTHON_PACKAGES})
            list(APPEND _uv_run_args --with ${_package})
        endforeach()
    endif()

    add_test(
        NAME ${TEST_NAME}
        COMMAND
            ${CMAKE_COMMAND} -DEXECUTABLE=${AddTest_EXECUTABLE_PARSED}
            "-DEXECUTABLE_ARGS=${_exe_args}" # Quoted because
            # passed as list see https://stackoverflow.com/a/33248574/80480
            -DBINARY_PATH=${_binary_path} -DWRAPPER_COMMAND=${WRAPPER_COMMAND}
            "-DWRAPPER_ARGS=${AddTest_WRAPPER_ARGS}"
            -DWORKING_DIRECTORY=${AddTest_WORKING_DIRECTORY}
            "-DLOG_ROOT=${PROJECT_BINARY_DIR}/logs"
            "-DLOG_FILE_BASENAME=${TEST_NAME}.txt"
            "-DTEST_COMMAND_IS_EXPECTED_TO_SUCCEED=${TEST_COMMAND_IS_EXPECTED_TO_SUCCEED}"
            "-DUV_RUN_ARGS=${_uv_run_args}"
            -P ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestWrapper.cmake
    )

    if(DEFINED AddTest_DEPENDS)
        set(_depends ${AddTest_DEPENDS})
        if(NOT (TEST ${_depends} OR TARGET ${_depends}))
            # If non-OMP variant was skipped, check if -omp variant exists
            if(NOT OGS_ENABLE_NON_OMP_TEST_VARIANTS AND TEST ${_depends}-omp)
                set(_depends ${_depends}-omp)
            else()
                message(
                    FATAL_ERROR
                        "AddTest ${TEST_NAME}: dependency ${AddTest_DEPENDS} does not exist!"
                )
            endif()
        endif()
        set_tests_properties(${TEST_NAME} PROPERTIES DEPENDS ${_depends})
    endif()
    if(DEFINED MPI_PROCESSORS)
        set_tests_properties(
            ${TEST_NAME} PROPERTIES PROCESSORS ${MPI_PROCESSORS}
        )
    endif()

    set_tests_properties(
        ${TEST_NAME}
        PROPERTIES ${AddTest_PROPERTIES}
                   COST
                   ${AddTest_RUNTIME}
                   DISABLED
                   ${AddTest_DISABLED}
                   LABELS
                   "${labels}"
                   ${timeout}
                   ENVIRONMENT
                   "PYDEVD_DISABLE_FILE_VALIDATION=1;UV_PYTHON=$ENV{UV_PYTHON};UV_PROJECT=$ENV{UV_PROJECT};UV_PROJECT_ENVIRONMENT=$ENV{UV_PROJECT_ENVIRONMENT}"
    )
endmacro()

macro(_add_np4_test_variant BASE_TEST_NAME)
    set(_wrapper_args_backup ${AddTest_WRAPPER_ARGS})
    set(_labels_backup ${labels})
    set(_exe_args_backup ${AddTest_EXECUTABLE_ARGS})
    set(_had_mpi_processors FALSE)
    if(DEFINED MPI_PROCESSORS)
        set(_had_mpi_processors TRUE)
        set(_processors_backup ${MPI_PROCESSORS})
    endif()

    set(_np4_wrapper_args "")
    set(_expect_np_value FALSE)
    foreach(_arg ${AddTest_WRAPPER_ARGS})
        if(_expect_np_value)
            list(APPEND _np4_wrapper_args 4)
            set(_expect_np_value FALSE)
            continue()
        endif()
        list(APPEND _np4_wrapper_args ${_arg})
        if(_arg STREQUAL "-np")
            set(_expect_np_value TRUE)
        endif()
    endforeach()

    set(_register_np4_variant TRUE)
    set(_np4_project_file "")
    set(_np4_cfg_ready FALSE)
    set(_np4_setup_test "")

    if("${AddTest_EXECUTABLE}" STREQUAL "ogs")
        _ogs_find_project_file_from_args(
            _np4_project_file "${AddTest_SOURCE_PATH}" ${AddTest_EXECUTABLE_ARGS}
        )
        if(NOT _np4_project_file STREQUAL "")
            _ogs_project_supports_partitions("${_np4_project_file}" 4
                                             _supports_np4_partitions
                                             _supports_np4_reason
            )
            if(NOT _supports_np4_partitions)
                set(_register_np4_variant FALSE)
            endif()
            _ogs_project_has_cfg4_for_all_meshes(
                "${_np4_project_file}" 4 _np4_cfg_ready
            )
        endif()
    endif()

    if(OGS_PETSC_NP4_VARIANTS_REQUIRE_CFG4_READY)
        if(_np4_project_file STREQUAL "")
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

    set(AddTest_WRAPPER_ARGS ${_np4_wrapper_args})
    set(MPI_PROCESSORS 4)
    set(_np4_labels ${labels})
    list(REMOVE_ITEM _np4_labels petsc_np1_source)
    list(APPEND _np4_labels petsc_np4_variant)
    set(labels ${_np4_labels})

    if(_register_np4_variant)
        if(NOT _np4_project_file STREQUAL "")
            set(_np4_mesh_dir "${AddTest_BINARY_PATH}-np4-meshes")
            list(APPEND AddTest_EXECUTABLE_ARGS -m ${_np4_mesh_dir})
            _ogs_add_np4_mesh_setup_test(
                _np4_setup_test "${BASE_TEST_NAME}-np4" "${_np4_project_file}"
                "${_np4_mesh_dir}" "${AddTest_DISABLED}"
            )
        endif()

        _add_test(${BASE_TEST_NAME}-np4)
        if(NOT _np4_setup_test STREQUAL "")
            _ogs_append_test_dependency("${BASE_TEST_NAME}-np4"
                                        "${_np4_setup_test}"
            )
        endif()

        if("${BASE_TEST_NAME}" MATCHES "-omp$")
            _set_omp_test_properties_for(${BASE_TEST_NAME}-np4)
        endif()
    endif()

    set(AddTest_EXECUTABLE_ARGS ${_exe_args_backup})

    set(AddTest_WRAPPER_ARGS ${_wrapper_args_backup})
    if(_had_mpi_processors)
        set(MPI_PROCESSORS ${_processors_backup})
    else()
        unset(MPI_PROCESSORS)
    endif()
    set(labels ${_labels_backup})
endmacro()

# Sets number of threads, adds label 'omp'
macro(_set_omp_test_properties_for TEST_NAME_OMP)
    get_test_property(${TEST_NAME_OMP} ENVIRONMENT _environment)
    if(NOT _environment)
        set(_environment "")
    endif()
    set_tests_properties(
        ${TEST_NAME_OMP}
        PROPERTIES ENVIRONMENT "OGS_ASM_THREADS=4;${_environment}" PROCESSORS 4
                   LABELS "${labels};omp"
    )
endmacro()

macro(_set_omp_test_properties)
    _set_omp_test_properties_for(${TEST_NAME}-omp)
endmacro()

# Adds subsequent tester ctest.
macro(_add_test_tester TEST_NAME)

    unset(TESTER_COMMAND)
    set(_binary_path ${AddTest_BINARY_PATH})
    if("${TEST_NAME}" MATCHES "-omp")
        set(_binary_path ${_binary_path}-omp)
    endif()
    if("${TEST_NAME}" MATCHES "-np4")
        set(_binary_path ${_binary_path}-np4)
    endif()

    set(TESTER_NAME "${TEST_NAME}-${AddTest_TESTER}")

    if(AddTest_TESTER STREQUAL "diff")
        foreach(FILE ${AddTest_DIFF_DATA})
            get_filename_component(FILE_EXPECTED ${FILE} NAME)
            list(
                APPEND
                TESTER_COMMAND
                "${SELECTED_DIFF_TOOL_PATH} \
                ${TESTER_ARGS} ${AddTest_TESTER_ARGS} ${AddTest_SOURCE_PATH}/${FILE_EXPECTED} \
                ${_binary_path}/${FILE}"
            )
        endforeach()
    elseif(AddTest_TESTER STREQUAL "vtkdiff")
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLengthMod4 "${DiffDataLength} % 4")
        math(EXPR DiffDataLengthMod6 "${DiffDataLength} % 6")
        if(${DiffDataLengthMod4} EQUAL 0 AND NOT ${DiffDataLengthMod6} EQUAL 0)
            message(WARNING "DEPRECATED AddTest call with four arguments.\
Use six arguments version of AddTest with absolute and relative tolerances"
            )
            if(NOT AddTest_ABSTOL)
                set(AddTest_ABSTOL 1e-16)
            endif()
            if(NOT AddTest_RELTOL)
                set(AddTest_RELTOL 1e-16)
            endif()
            set(TESTER_ARGS "--abs ${AddTest_ABSTOL} --rel ${AddTest_RELTOL}")
            math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
            foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 4)
                list(GET AddTest_DIFF_DATA "${DiffDataIndex}"
                     REFERENCE_VTK_FILE
                )
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" VTK_FILE)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_A)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+3")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_B)

                list(
                    APPEND
                    TESTER_COMMAND
                    "${SELECTED_DIFF_TOOL_PATH} \
                ${AddTest_SOURCE_PATH}/${REFERENCE_VTK_FILE} \
                ${_binary_path}/${VTK_FILE} \
                -a ${NAME_A} -b ${NAME_B} \
                ${TESTER_ARGS} ${AddTest_TESTER_ARGS}"
                )
            endforeach()
        elseif(${DiffDataLengthMod6} EQUAL 0)
            if(${AddTest_ABSTOL} OR ${AddTest_RELTOL})
                message(
                    FATAL_ERROR
                        "ABSTOL or RELTOL arguments must not be present."
                )
            endif()
            math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
            foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 6)
                list(GET AddTest_DIFF_DATA "${DiffDataIndex}"
                     REFERENCE_VTK_FILE
                )
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" VTK_FILE)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_A)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+3")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_B)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+4")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" ABS_TOL)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+5")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" REL_TOL)

                if("${REFERENCE_VTK_FILE}" STREQUAL "GLOB")
                    list(APPEND TESTER_COMMAND
                         "${VTK_FILE} ${NAME_A} ${NAME_B} ${ABS_TOL} ${REL_TOL}"
                    )
                    set(GLOB_MODE TRUE)
                else()
                    list(
                        APPEND
                        TESTER_COMMAND
                        "${SELECTED_DIFF_TOOL_PATH} \
                    ${AddTest_SOURCE_PATH}/${REFERENCE_VTK_FILE} \
                    ${_binary_path}/${VTK_FILE} \
                    -a ${NAME_A} -b ${NAME_B} \
                    --abs ${ABS_TOL} --rel ${REL_TOL} \
                    ${TESTER_ARGS} ${AddTest_TESTER_ARGS}"
                    )
                endif()
            endforeach()
        else()
            message(
                FATAL_ERROR
                    "For vtkdiff tester the number of diff data arguments must be a multiple of six."
            )
        endif()
    elseif(AddTest_TESTER STREQUAL "xdmfdiff")
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLengthMod8 "${DiffDataLength} % 8")
        if(NOT ${DiffDataLengthMod8} EQUAL 0)
            message(
                FATAL_ERROR
                    "For xdmfdiff tester the number of diff data arguments must be a multiple of eight."
            )
        endif()
        if(${AddTest_ABSTOL} OR ${AddTest_RELTOL})
            message(
                FATAL_ERROR
                    "ABSTOL or RELTOL arguments must not be present."
            )
        endif()
        math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
        foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 8)
            list(GET AddTest_DIFF_DATA "${DiffDataIndex}"
                 REFERENCE_VTK_FILE
            )
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" VTK_FILE)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_A)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+3")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" NAME_B)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+4")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" ABS_TOL)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+5")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" REL_TOL)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+6")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" TIMESTEP_A)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+7")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" TIMESTEP_B)

            list(
                APPEND
                TESTER_COMMAND
                "${SELECTED_DIFF_TOOL_PATH} \
                ${AddTest_SOURCE_PATH}/${REFERENCE_VTK_FILE} \
                ${_binary_path}/${VTK_FILE} \
                -a ${NAME_A} -b ${NAME_B} \
                --abs ${ABS_TOL} --rel ${REL_TOL} \
                --timestep-a ${TIMESTEP_A} --timestep-b ${TIMESTEP_B}\
                ${TESTER_ARGS} ${AddTest_TESTER_ARGS}"
            )
        endforeach()
    elseif(AddTest_TESTER STREQUAL "vtkdiff-mesh")
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLengthMod3 "${DiffDataLength} % 3")
        if(${DiffDataLengthMod3} EQUAL 0)
            math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
            foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 3)
                list(GET AddTest_DIFF_DATA "${DiffDataIndex}"
                     REFERENCE_VTK_FILE
                )
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" VTK_FILE)
                math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
                list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" ABS_TOLERANCE)

                list(
                    APPEND
                    TESTER_COMMAND
                    "${SELECTED_DIFF_TOOL_PATH} -m \
                ${AddTest_SOURCE_PATH}/${REFERENCE_VTK_FILE} \
                ${_binary_path}/${VTK_FILE} \
                --abs ${ABS_TOLERANCE}"
                )
            endforeach()
        else()
            message(FATAL_ERROR "The number of diff data arguments must be a
            multiple of three: expected.vtu output.vtu absolute_tolerance."
            )
        endif()
    elseif(AddTest_TESTER STREQUAL "gmldiff")
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
        foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 3)
            list(GET AddTest_DIFF_DATA "${DiffDataIndex}" GML_FILE)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" ABS_TOL)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" REL_TOL)

            get_filename_component(FILE_EXPECTED ${GML_FILE} NAME)
            if(WIN32)
                file(TO_NATIVE_PATH "${Python_EXECUTABLE}" PY_EXE)
                # Dirty hack for Windows Python paths with spaces:
                string(REPLACE "Program Files" "\"Program Files\"" PY_EXE
                               ${PY_EXE}
                )
            else()
                set(PY_EXE ${Python_EXECUTABLE})
            endif()
            list(
                APPEND
                TESTER_COMMAND
                "${PY_EXE} ${PROJECT_SOURCE_DIR}/scripts/test/gmldiff.py \
                --abs ${ABS_TOL} --rel ${REL_TOL} \
                ${TESTER_ARGS} ${AddTest_TESTER_ARGS}\
                ${AddTest_SOURCE_PATH}/${FILE_EXPECTED} \
                ${_binary_path}/${GML_FILE}"
            )
        endforeach()
    elseif(AddTest_TESTER STREQUAL "memcheck")
        set(TESTER_COMMAND
            "! ${GREP_TOOL_PATH} definitely ${AddTest_SOURCE_PATH}/${AddTest_NAME}_memcheck.txt"
        )
    elseif(AddTest_TESTER STREQUAL "numdiff")
        list(LENGTH AddTest_DIFF_DATA DiffDataLength)
        math(EXPR DiffDataLastIndex "${DiffDataLength}-1")
        foreach(DiffDataIndex RANGE 0 ${DiffDataLastIndex} 3)
            list(GET AddTest_DIFF_DATA "${DiffDataIndex}" FILE)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+1")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" ABS_TOL)
            math(EXPR DiffDataAuxIndex "${DiffDataIndex}+2")
            list(GET AddTest_DIFF_DATA "${DiffDataAuxIndex}" REL_TOL)
            get_filename_component(FILE_EXPECTED ${FILE} NAME)
            list(
                APPEND
                TESTER_COMMAND
                "${SELECTED_DIFF_TOOL_PATH} -a ${ABS_TOL} -r ${REL_TOL}  \
        ${TESTER_ARGS} ${AddTest_TESTER_ARGS} ${AddTest_SOURCE_PATH}/${FILE_EXPECTED} \
        ${_binary_path}/${FILE}"
            )
        endforeach()
    endif()

    add_test(
        NAME ${TESTER_NAME}
        COMMAND
            ${CMAKE_COMMAND} -DSOURCE_PATH=${AddTest_SOURCE_PATH}
            -DSELECTED_DIFF_TOOL_PATH=${SELECTED_DIFF_TOOL_PATH}
            "-DTESTER_COMMAND=${TESTER_COMMAND}" -DBINARY_PATH=${_binary_path}
            -DGLOB_MODE=${GLOB_MODE}
            -DLOG_FILE_BASE=${PROJECT_BINARY_DIR}/logs/${TESTER_NAME} -P
            ${PROJECT_SOURCE_DIR}/scripts/cmake/test/AddTestTester.cmake
            --debug-output
        WORKING_DIRECTORY ${AddTest_SOURCE_PATH}
    )
    set(_tester_labels ${labels})
    if(NOT OGS_PETSC_NP4_INCLUDE_TESTER_VARIANTS)
        list(REMOVE_ITEM _tester_labels petsc_np1_source petsc_np4_variant)
    elseif("${TEST_NAME}" MATCHES "-np4")
        list(REMOVE_ITEM _tester_labels petsc_np1_source)
        list(APPEND _tester_labels petsc_np4_variant)
    endif()
    set_tests_properties(
        ${TESTER_NAME} PROPERTIES DEPENDS ${TEST_NAME} DISABLED
                                  ${AddTest_DISABLED} LABELS
                                  "tester;${_tester_labels}"
    )
endmacro()

# Checks if a test is expected to succeed based on the properties WILL_FAIL,
# PASS_REGULAR_EXPRESSION and FAIL_REGULAR_EXPRESSION. The function expects the
# test name (used only for debugging purposes) and the test properties as
# arguments. The test does not need to exist, yet. This function does not query
# any test case, but only uses the passed list of properties
function(isTestCommandExpectedToSucceed TEST_NAME)
    set(options WILL_FAIL)
    set(oneValueArgs PASS_REGULAR_EXPRESSION FAIL_REGULAR_EXPRESSION)
    set(multiValueArgs)
    cmake_parse_arguments(
        TEST_FAILURE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN}
    )

    message(DEBUG "failure properties for test ${TEST_NAME}:")
    list(APPEND CMAKE_MESSAGE_INDENT "  ")
    message(DEBUG "WILL_FAIL: ${TEST_FAILURE_WILL_FAIL}")
    message(DEBUG "PASS_RE: ${TEST_FAILURE_PASS_REGULAR_EXPRESSION}")
    message(DEBUG "FAIL_RE: ${TEST_FAILURE_FAIL_REGULAR_EXPRESSION}")
    list(POP_BACK CMAKE_MESSAGE_INDENT)

    if(${TEST_FAILURE_WILL_FAIL})
        if(DEFINED TEST_FAILURE_PASS_REGULAR_EXPRESSION)
            # Note: if the test property PASS_REGULAR_EXPRESSION is set, the
            # process return code will be ignored, see
            # https://cmake.org/cmake/help/latest/prop_test/PASS_REGULAR_EXPRESSION.html
            message(
                SEND_ERROR
                    "Error in test '${TEST_NAME}': Please do not use both WILL_FAIL and PASS_REGULAR_EXPRESSION in the same test. The logic will be unclear, then."
            )
        endif()
        if(DEFINED TEST_FAILURE_FAIL_REGULAR_EXPRESSION)
            message(
                SEND_ERROR
                    "Error in test '${TEST_NAME}': Please do not use both WILL_FAIL and FAIL_REGULAR_EXPRESSION in the same test. The logic will be unclear, then."
            )
        endif()

        set(TEST_COMMAND_IS_EXPECTED_TO_SUCCEED false)
    elseif(DEFINED TEST_FAILURE_PASS_REGULAR_EXPRESSION)
        if(DEFINED TEST_FAILURE_FAIL_REGULAR_EXPRESSION)
            message(
                SEND_ERROR
                    "Error in test '${TEST_NAME}': Please do not use both PASS_REGULAR_EXPRESSION and FAIL_REGULAR_EXPRESSION in the same test. The logic will be unclear, then."
            )
        endif()

        set(TEST_COMMAND_IS_EXPECTED_TO_SUCCEED false)
    else()
        set(TEST_COMMAND_IS_EXPECTED_TO_SUCCEED true)
    endif()

    set(TEST_COMMAND_IS_EXPECTED_TO_SUCCEED
        "${TEST_COMMAND_IS_EXPECTED_TO_SUCCEED}" PARENT_SCOPE
    )
endfunction()
