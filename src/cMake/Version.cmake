INCLUDE(Macros)

MACRO(GENERATE_VERSION_H)
    GETDATETIME(BUILD_TIME "%Y-%m-%d %H:%M:%S")

    EXECUTE_PROCESS(
      COMMAND git log -1 --format=%H
      WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
      OUTPUT_VARIABLE GIT_HASH
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # mark builds from a modified working tree (P0-10)
    EXECUTE_PROCESS(
      COMMAND git diff-index --quiet HEAD --
      WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
      RESULT_VARIABLE GIT_DIRTY
      OUTPUT_QUIET ERROR_QUIET
    )
    IF(NOT "${GIT_DIRTY}" STREQUAL "0")
      SET(GIT_HASH "${GIT_HASH}-dirty")
    ENDIF()

    # report the compiled contact-model whitelist, not only the ENABLE_MODEL_* options (P0-10)
    IF(LIGGGHTS_CONTACT_WHITELIST_FROM_OPTIONS)
      SET(CONTACT_MODELS_BANNER "models:${ENABLED_MODELS}, contact-model whitelist: ${CONTACT_WHITELIST_COUNT} combinations from ENABLE_MODEL_* options")
    ELSE()
      GET_FILENAME_COMPONENT(CONTACT_WHITELIST_NAME "${LIGGGHTS_CONTACT_WHITELIST}" NAME)
      SET(CONTACT_MODELS_BANNER "contact-model whitelist: ${CONTACT_WHITELIST_COUNT} combinations from ${CONTACT_WHITELIST_NAME}")
    ENDIF()

    EXECUTE_PROCESS(
      COMMAND whoami
      WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
      OUTPUT_VARIABLE USER
      OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    STRING(STRIP ${USER} USER)
    STRING(REPLACE "\\" "\\\\" USER ${USER})

    FILE(WRITE version_liggghts.h
      "#define LIGGGHTS_VERSION \"LIGGGHTS-PUBLIC ${LIGGGHTS_VERSION}, "
      "compiled ${BUILD_TIME} by ${USER}, git commit ${GIT_HASH}, "
      "build configuration:${ENABLED_OPTIONS}, "
      "${CONTACT_MODELS_BANNER} (+ runtime fallback)\""
    )
ENDMACRO()