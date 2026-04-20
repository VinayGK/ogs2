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

    set(${OUT_PROJECT_FILE} "${_resolved_project_file}" PARENT_SCOPE)
endfunction()

function(_collect_project_meshes PROJECT_FILE OUT_MESHES)
    file(STRINGS "${PROJECT_FILE}" _project_lines)
    set(_inside_meshes_block FALSE)
    set(_before_process_sections TRUE)
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

        if(NOT _line MATCHES "<mesh[^>]*>[^<]+</mesh>")
            continue()
        endif()

        if(NOT _inside_meshes_block AND NOT _before_process_sections)
            continue()
        endif()

        string(REGEX REPLACE ".*<mesh[^>]*>([^<]+)</mesh>.*" "\\1" _mesh "${_line}")
        string(STRIP "${_mesh}" _mesh)
        if(NOT _mesh STREQUAL "")
            list(APPEND _meshes "${_mesh}")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES _meshes)
    set(${OUT_MESHES} "${_meshes}" PARENT_SCOPE)
endfunction()

function(_resolve_mesh_source PROJECT_DIR MESH_REFERENCE OUT_MESH_SOURCE)
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

_resolve_project_file("${PROJECT_FILE}" _resolved_project_file)
if(NOT _resolved_project_file STREQUAL "${PROJECT_FILE}")
    message(
        STATUS
            "Resolved project-diff '${PROJECT_FILE}' to '${_resolved_project_file}' for mesh preparation."
    )
endif()

get_filename_component(_project_dir "${_resolved_project_file}" DIRECTORY)
_collect_project_meshes("${_resolved_project_file}" _project_meshes)
if(NOT _project_meshes)
    message(
        WARNING
            "No top-level <mesh> entries found in ${_resolved_project_file}; skipping mesh preparation."
    )
    return()
endif()

set(_mesh_sources "")
set(_missing_cfg FALSE)
set(_resolved_mesh_count 0)

foreach(_mesh IN LISTS _project_meshes)
    _resolve_mesh_source("${_project_dir}" "${_mesh}" _mesh_source)
    if(_mesh_source STREQUAL "")
        message(
            WARNING
                "Could not resolve top-level mesh reference '${_mesh}' in ${_resolved_project_file}; skipping this mesh."
        )
        continue()
    endif()
    math(EXPR _resolved_mesh_count "${_resolved_mesh_count} + 1")

    get_filename_component(_mesh_source_dir "${_mesh_source}" DIRECTORY)
    get_filename_component(_mesh_base "${_mesh_source}" NAME_WE)
    set(_source_cfg
        "${_mesh_source_dir}/${_mesh_base}_partitioned_msh_cfg${N_PARTITIONS}.bin"
    )
    set(_out_cfg
        "${MESH_OUTPUT_DIR}/${_mesh_base}_partitioned_msh_cfg${N_PARTITIONS}.bin"
    )

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
        set(_missing_cfg TRUE)
    endif()

    list(APPEND _mesh_sources "${_mesh_source}")
endforeach()

if(_resolved_mesh_count EQUAL 0)
    message(
        WARNING
            "No resolvable mesh files were found in ${_resolved_project_file}; skipping mesh preparation."
    )
    return()
endif()

list(GET _mesh_sources 0 _main_mesh)

if(_missing_cfg)
    if(NOT DEFINED PARTMESH_EXECUTABLE OR PARTMESH_EXECUTABLE STREQUAL "")
        message(
            WARNING
                "Missing cfg${N_PARTITIONS} partitioned meshes and PARTMESH_EXECUTABLE is not available; skipping generation."
        )
        return()
    endif()
    if(NOT EXISTS "${PARTMESH_EXECUTABLE}")
        message(
            WARNING
                "PARTMESH_EXECUTABLE does not exist: ${PARTMESH_EXECUTABLE}; skipping generation."
        )
        return()
    endif()

    list(REMOVE_AT _mesh_sources 0)

    message(
        STATUS
            "Preparing cfg${N_PARTITIONS} partitioned meshes for ${PROJECT_FILE} in ${MESH_OUTPUT_DIR}"
    )

    execute_process(
        COMMAND ${PARTMESH_EXECUTABLE} --ogs2metis -i "${_main_mesh}" -o
                "${MESH_OUTPUT_DIR}"
        RESULT_VARIABLE _partmesh_ogs2metis_exit
        OUTPUT_VARIABLE _partmesh_ogs2metis_out
        ERROR_VARIABLE _partmesh_ogs2metis_err
    )
    if(NOT _partmesh_ogs2metis_exit EQUAL 0)
        message(
            WARNING
                "partmesh --ogs2metis failed with exit code ${_partmesh_ogs2metis_exit}\n${_partmesh_ogs2metis_out}\n${_partmesh_ogs2metis_err}"
        )
        return()
    endif()

    set(_partmesh_cmd ${PARTMESH_EXECUTABLE} --exe_metis -n ${N_PARTITIONS}
                      -i "${_main_mesh}" -o "${MESH_OUTPUT_DIR}"
    )
    if(_mesh_sources)
        list(APPEND _partmesh_cmd -- ${_mesh_sources})
    endif()

    execute_process(
        COMMAND ${_partmesh_cmd}
        RESULT_VARIABLE _partmesh_partition_exit
        OUTPUT_VARIABLE _partmesh_partition_out
        ERROR_VARIABLE _partmesh_partition_err
    )
    if(NOT _partmesh_partition_exit EQUAL 0)
        message(
            WARNING
                "partmesh partitioning failed with exit code ${_partmesh_partition_exit}\n${_partmesh_partition_out}\n${_partmesh_partition_err}"
        )
        return()
    endif()
endif()

get_filename_component(_main_mesh_base "${_main_mesh}" NAME_WE)
set(_main_out_cfg
    "${MESH_OUTPUT_DIR}/${_main_mesh_base}_partitioned_msh_cfg${N_PARTITIONS}.bin"
)
if(NOT EXISTS "${_main_out_cfg}")
    message(
        WARNING
            "Partitioned cfg${N_PARTITIONS} main mesh was not created: ${_main_out_cfg}"
    )
endif()
