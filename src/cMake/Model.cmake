# based on: https://github.com/schrummy14/LIGGGHTS_Flexible_Fibers/blob/master/src/WINDOWS/CMake_patch.zip

INCLUDE(Macros)

# Source of the static contact-model whitelist (style_contact_model.h).
#  default : the tracked list src/contact_model_whitelist.txt, plus the local
#            src/style_contact_model_user.whitelist if it exists (same files as
#            the Make route, so both routes compile the same combinations)
#  -DLIGGGHTS_CONTACT_WHITELIST_FROM_OPTIONS=ON : legacy behaviour, the cross
#            product of the enabled ENABLE_MODEL_* options
# Combinations that are not compiled still run through the runtime fallback
# (with a one-time warning), unless LIGGGHTS_NO_CONTACT_MODEL_FALLBACK is defined.
SET(LIGGGHTS_CONTACT_WHITELIST "${CMAKE_CURRENT_SOURCE_DIR}/contact_model_whitelist.txt"
    CACHE FILEPATH "Tracked contact-model whitelist (GRAN_MODEL(...) lines) used to generate style_contact_model.h")
OPTION(LIGGGHTS_CONTACT_WHITELIST_FROM_OPTIONS
       "Generate the contact-model whitelist as the cross product of the ENABLE_MODEL_* options instead of reading LIGGGHTS_CONTACT_WHITELIST" OFF)

MACRO(WRITE_WHITELIST)
    SET(fileName ${CMAKE_CURRENT_SOURCE_DIR}/style_contact_model.h)
    SET(tmpFileName ${CMAKE_CURRENT_BINARY_DIR}/style_contact_model.h.tmp)
    # Create File
    GETDATETIME(NOW "%Y-%m-%d %H:%M:%S")
    FILE(WRITE ${tmpFileName} "/* created on ${NOW} */\n")

    SET(N 0)
    IF(LIGGGHTS_CONTACT_WHITELIST_FROM_OPTIONS)
        SET(CONTACT_WHITELIST_SOURCE "ENABLE_MODEL_* options")
        GET_PROPERTY(normal_models_local GLOBAL PROPERTY normal_models)
        GET_PROPERTY(tangential_models_local GLOBAL PROPERTY tangential_models)
        GET_PROPERTY(cohesion_models_local GLOBAL PROPERTY cohesion_models)
        GET_PROPERTY(rolling_models_local GLOBAL PROPERTY rolling_models)
        GET_PROPERTY(surface_models_local GLOBAL PROPERTY surface_models)

        FOREACH(nm ${normal_models_local})
            FOREACH(tm ${tangential_models_local})
                FOREACH(cm ${cohesion_models_local})
                    FOREACH(rm ${rolling_models_local})
                        FOREACH(sm ${surface_models_local})
                            MATH(EXPR N "${N}+1")
                            FILE(APPEND ${tmpFileName} "GRAN_MODEL(${nm}, ${tm}, ${cm}, ${rm}, ${sm})\n")
                        ENDFOREACH()
                    ENDFOREACH()
                ENDFOREACH()
            ENDFOREACH()
        ENDFOREACH()
    ELSE()
        IF(NOT EXISTS "${LIGGGHTS_CONTACT_WHITELIST}")
            MESSAGE(FATAL_ERROR "Contact-model whitelist ${LIGGGHTS_CONTACT_WHITELIST} not found. "
                                "Set -DLIGGGHTS_CONTACT_WHITELIST=<file> or -DLIGGGHTS_CONTACT_WHITELIST_FROM_OPTIONS=ON.")
        ENDIF()
        SET(CONTACT_WHITELIST_SOURCE "${LIGGGHTS_CONTACT_WHITELIST}")
        SET(whitelist_files "${LIGGGHTS_CONTACT_WHITELIST}")
        IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/style_contact_model_user.whitelist")
            LIST(APPEND whitelist_files "${CMAKE_CURRENT_SOURCE_DIR}/style_contact_model_user.whitelist")
            SET(CONTACT_WHITELIST_SOURCE "${CONTACT_WHITELIST_SOURCE} + style_contact_model_user.whitelist")
        ENDIF()
        # re-run configure when a whitelist file changes
        SET_PROPERTY(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${whitelist_files})
        SET(seen_entries "")
        FOREACH(wlf ${whitelist_files})
            FILE(STRINGS "${wlf}" wl_lines REGEX "^[ \t]*GRAN_MODEL[ \t]*\\(")
            FOREACH(line ${wl_lines})
                STRING(STRIP "${line}" line)
                # normalise whitespace so duplicates are detected
                STRING(REGEX REPLACE "[ \t]+" "" key "${line}")
                LIST(FIND seen_entries "${key}" idx)
                IF(idx EQUAL -1)
                    LIST(APPEND seen_entries "${key}")
                    MATH(EXPR N "${N}+1")
                    FILE(APPEND ${tmpFileName} "${line}\n")
                ENDIF()
            ENDFOREACH()
        ENDFOREACH()
    ENDIF()

    # only touch style_contact_model.h when its entries change, so that a
    # re-configure does not force the ~45 s recompile of lammps.cpp (PF-07)
    SET(write_new TRUE)
    IF(EXISTS ${fileName})
        FILE(STRINGS ${fileName} old_entries REGEX "^GRAN_MODEL")
        FILE(STRINGS ${tmpFileName} new_entries REGEX "^GRAN_MODEL")
        IF("${old_entries}" STREQUAL "${new_entries}")
            SET(write_new FALSE)
        ENDIF()
    ENDIF()
    IF(write_new)
        CONFIGURE_FILE(${tmpFileName} ${fileName} COPYONLY)
    ENDIF()
    FILE(REMOVE ${tmpFileName})

    SET(CONTACT_WHITELIST_COUNT ${N})
    MESSAGE(STATUS "There are ${N} contact model combinations (source: ${CONTACT_WHITELIST_SOURCE})")
ENDMACRO()

MACRO(ADD_CONTACT_MODEL normal_model tangential_model cohesion_model rolling_model surface_model)
    GET_PROPERTY(normal_models_local GLOBAL PROPERTY normal_models)
    GET_PROPERTY(tangential_models_local GLOBAL PROPERTY tangential_models)
    GET_PROPERTY(cohesion_models_local GLOBAL PROPERTY cohesion_models)
    GET_PROPERTY(rolling_models_local GLOBAL PROPERTY rolling_models)
    GET_PROPERTY(surface_models_local GLOBAL PROPERTY surface_models)

    SET_PROPERTY(GLOBAL PROPERTY normal_models "${normal_models_local}" "${normal_model}")
    SET_PROPERTY(GLOBAL PROPERTY tangential_models "${tangential_models_local}" "${tangential_model}")
    SET_PROPERTY(GLOBAL PROPERTY cohesion_models "${cohesion_models_local}" "${cohesion_model}")
    SET_PROPERTY(GLOBAL PROPERTY rolling_models "${rolling_models_local}" "${rolling_model}")
    SET_PROPERTY(GLOBAL PROPERTY surface_models "${surface_models_local}" "${surface_model}")
ENDMACRO()
