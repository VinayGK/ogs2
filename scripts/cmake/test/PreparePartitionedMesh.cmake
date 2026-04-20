# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

if(NOT DEFINED PROJECT_FILE)
    message(FATAL_ERROR "PROJECT_FILE was not provided.")
endif()
if(NOT DEFINED MESH_OUTPUT_DIR)
    message(FATAL_ERROR "MESH_OUTPUT_DIR was not provided.")
endif()
if(NOT DEFINED N_PARTITIONS)
    message(FATAL_ERROR "N_PARTITIONS was not provided.")
endif()

if(NOT EXISTS "${PROJECT_FILE}")
    message(FATAL_ERROR "Project file does not exist: ${PROJECT_FILE}")
endif()

file(MAKE_DIRECTORY "${MESH_OUTPUT_DIR}")

function(_resolve_project_file INPUT_PROJECT_FILE OUT_PROJECT_FILE)
    set(_resolved_project_file "${INPUT_PROJECT_FILE}")
    file(READ "${INPUT_PROJECT_FILE}" _project_xml LIMIT 4096)

    if(_project_xml MATCHES "<OpenGeoSysProjectDiff")
        get_filename_component(_project_dir "${INPUT_PROJECT_FILE}" DIRECTORY)
        string(REGEX MATCH "base_file=\"([^\"]+)\"" _base_file_match
                           "${_project_xml}"
        )
        set(_base_file "${CMAKE_MATCH_1}")
        if(NOT _base_file STREQUAL "")
            if(IS_ABSOLUTE "${_base_file}")
                set(_base_candidate "${_base_file}")
            else()
                set(_base_candidate "${_project_dir}/${_base_file}")
            endif()
            if(EXISTS "${_base_candidate}")
                set(_resolved_project_file "${_base_candidate}")
            endif()
        endif()

        if(_resolved_project_file STREQUAL "${INPUT_PROJECT_FILE}")
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

function(_collect_project_meshes PROJECT_FILE OUT_MESHES)
    file(STRINGS "${PROJECT_FILE}" _project_lines)
    set(_inside_meshes_block FALSE)
    set(_before_process_sections TRUE)
    set(_inside_mesh_tag FALSE)
    set(_pending_mesh "")
    set(_meshes "")
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
                    list(APPEND _meshes "${_mesh}")
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
                list(APPEND _meshes "${_mesh}")
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
                    list(APPEND _meshes "${_mesh}")
                endif()
            else()
                set(_inside_mesh_tag TRUE)
            endif()
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _meshes)
    set(${OUT_MESHES} "${_meshes}" PARENT_SCOPE)
endfunction()

function(_collect_projectdiff_mesh_replacements PROJECT_FILE OUT_MESHES)
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

function(_resolve_mesh_source PROJECT_DIR RESOLVED_PROJECT_DIR MESH_REFERENCE OUT_MESH_SOURCE)
    set(_mesh_source "")
    if(IS_ABSOLUTE "${MESH_REFERENCE}" AND EXISTS "${MESH_REFERENCE}")
        set(_mesh_source "${MESH_REFERENCE}")
    else()
        set(_candidate_dirs "${PROJECT_DIR}")
        if(NOT "${RESOLVED_PROJECT_DIR}" STREQUAL "${PROJECT_DIR}")
            list(APPEND _candidate_dirs "${RESOLVED_PROJECT_DIR}")
        endif()

        foreach(_dir IN LISTS _candidate_dirs)
            if(NOT _mesh_source STREQUAL "")
                break()
            endif()

            set(_candidate "${_dir}/${MESH_REFERENCE}")
            if(EXISTS "${_candidate}")
                set(_mesh_source "${_candidate}")
                break()
            endif()

            get_filename_component(_mesh_ext "${MESH_REFERENCE}" EXT)
            if(_mesh_ext STREQUAL "")
                foreach(_ext .vtu .msh .vtk .vtm .vti)
                    set(_candidate_with_ext "${_dir}/${MESH_REFERENCE}${_ext}")
                    if(EXISTS "${_candidate_with_ext}")
                        set(_mesh_source "${_candidate_with_ext}")
                        break()
                    endif()
                endforeach()
            endif()
        endforeach()
    endif()

    set(${OUT_MESH_SOURCE} "${_mesh_source}" PARENT_SCOPE)
endfunction()

function(
    _partition_meshes_for_cfg
    MAIN_MESH_SOURCE
    OTHER_MESH_SOURCES
    OUTPUT_DIR
    N_PARTITIONS
    PARTMESH_EXECUTABLE
    OUT_OK
    OUT_ERROR
)
    execute_process(
        COMMAND
            ${PARTMESH_EXECUTABLE} --ogs2metis -i "${MAIN_MESH_SOURCE}" -o
            "${OUTPUT_DIR}" ${OTHER_MESH_SOURCES}
        RESULT_VARIABLE _ogs2metis_exit
        OUTPUT_VARIABLE _ogs2metis_out
        ERROR_VARIABLE _ogs2metis_err
    )
    if(NOT _ogs2metis_exit EQUAL 0)
        set(${OUT_OK} FALSE PARENT_SCOPE)
        set(
            ${OUT_ERROR}
            "partmesh --ogs2metis failed (${_ogs2metis_exit}) for '${MAIN_MESH_SOURCE}'\n${_ogs2metis_out}\n${_ogs2metis_err}"
            PARENT_SCOPE
        )
        return()
    endif()

    execute_process(
        COMMAND
            ${PARTMESH_EXECUTABLE} --exe_metis -n ${N_PARTITIONS} -i
            "${MAIN_MESH_SOURCE}" -o "${OUTPUT_DIR}" ${OTHER_MESH_SOURCES}
        RESULT_VARIABLE _partition_exit
        OUTPUT_VARIABLE _partition_out
        ERROR_VARIABLE _partition_err
    )
    if(NOT _partition_exit EQUAL 0)
        set(${OUT_OK} FALSE PARENT_SCOPE)
        set(
            ${OUT_ERROR}
            "partmesh --exe_metis failed (${_partition_exit}) for '${MAIN_MESH_SOURCE}'\n${_partition_out}\n${_partition_err}"
            PARENT_SCOPE
        )
        return()
    endif()

    set(${OUT_OK} TRUE PARENT_SCOPE)
    set(${OUT_ERROR} "" PARENT_SCOPE)
endfunction()

_resolve_project_file("${PROJECT_FILE}" _resolved_project_file)
get_filename_component(_project_dir "${PROJECT_FILE}" DIRECTORY)
get_filename_component(_resolved_project_dir "${_resolved_project_file}" DIRECTORY)

set(_project_meshes "")
file(READ "${PROJECT_FILE}" _input_project_xml LIMIT 4096)
if(_input_project_xml MATCHES "<OpenGeoSysProjectDiff")
    _collect_projectdiff_mesh_replacements("${PROJECT_FILE}" _project_meshes)
endif()
if(NOT _project_meshes)
    _collect_project_meshes("${_resolved_project_file}" _project_meshes)
endif()

if(NOT _project_meshes)
    message(FATAL_ERROR "No mesh references found in '${PROJECT_FILE}'.")
endif()

set(_mesh_sources "")
set(_mesh_sources_to_partition "")
set(_expected_cfg_files "")

foreach(_mesh IN LISTS _project_meshes)
    _resolve_mesh_source("${_project_dir}" "${_resolved_project_dir}" "${_mesh}"
                         _mesh_source
    )
    if(_mesh_source STREQUAL "")
        message(
            FATAL_ERROR
                "Could not resolve mesh reference '${_mesh}' from '${PROJECT_FILE}'."
        )
    endif()

    get_filename_component(_mesh_source_dir "${_mesh_source}" DIRECTORY)
    get_filename_component(_mesh_name "${_mesh_source}" NAME)
    string(REGEX REPLACE "\\.[^.]+$" "" _mesh_base "${_mesh_name}")
    set(_source_cfg
        "${_mesh_source_dir}/${_mesh_base}_partitioned_msh_cfg${N_PARTITIONS}.bin"
    )
    set(_out_cfg
        "${MESH_OUTPUT_DIR}/${_mesh_base}_partitioned_msh_cfg${N_PARTITIONS}.bin"
    )
    list(APPEND _expected_cfg_files "${_out_cfg}")

    if(NOT EXISTS "${_out_cfg}")
        if(EXISTS "${_source_cfg}")
            file(
                GLOB _existing_partition_files
                "${_mesh_source_dir}/${_mesh_base}_partitioned_*${N_PARTITIONS}.bin"
            )
            if(_existing_partition_files)
                file(COPY ${_existing_partition_files} DESTINATION "${MESH_OUTPUT_DIR}")
            endif()
        endif()
    endif()

    if(NOT EXISTS "${_out_cfg}")
        list(APPEND _mesh_sources_to_partition "${_mesh_source}")
    endif()

    list(APPEND _mesh_sources "${_mesh_source}")
endforeach()

if(_mesh_sources_to_partition)
    if(NOT DEFINED PARTMESH_EXECUTABLE OR PARTMESH_EXECUTABLE STREQUAL "")
        message(
            FATAL_ERROR
                "Missing cfg${N_PARTITIONS} meshes and PARTMESH_EXECUTABLE is not available."
        )
    endif()
    if(NOT EXISTS "${PARTMESH_EXECUTABLE}")
        message(
            FATAL_ERROR
                "PARTMESH_EXECUTABLE does not exist: ${PARTMESH_EXECUTABLE}"
        )
    endif()

    message(
        STATUS
            "Preparing cfg${N_PARTITIONS} partitioned meshes for ${PROJECT_FILE} in ${MESH_OUTPUT_DIR}"
    )

    list(GET _mesh_sources 0 _main_mesh_source)
    set(_other_mesh_sources ${_mesh_sources})
    list(REMOVE_AT _other_mesh_sources 0)

    _partition_meshes_for_cfg(
        "${_main_mesh_source}" "${_other_mesh_sources}" "${MESH_OUTPUT_DIR}"
        "${N_PARTITIONS}" "${PARTMESH_EXECUTABLE}" _partition_ok _partition_error
    )

    if(NOT _partition_ok)
        message(
            FATAL_ERROR
                "Failed to prepare cfg${N_PARTITIONS} meshes:\n${_partition_error}"
        )
    endif()
endif()

set(_missing_cfg_outputs "")
foreach(_cfg_file IN LISTS _expected_cfg_files)
    if(NOT EXISTS "${_cfg_file}")
        list(APPEND _missing_cfg_outputs "${_cfg_file}")
    endif()
endforeach()

if(_missing_cfg_outputs)
    string(REPLACE ";" "\n" _missing_cfg_outputs_joined "${_missing_cfg_outputs}")
    message(
        FATAL_ERROR
            "Expected cfg${N_PARTITIONS} mesh files are still missing:\n${_missing_cfg_outputs_joined}"
    )
endif()
