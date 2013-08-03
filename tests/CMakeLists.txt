# Author: petter.strandmark@gmail.com (Petter Strandmark)

MACRO (CREATE_TEST NAME)
	ADD_EXECUTABLE(test_${NAME} 
	               test_${NAME}.cpp
	               ${SPII_HEADERS})

	TARGET_LINK_LIBRARIES(test_${NAME} easyip)

	ADD_TEST(NAME ${NAME}
	         COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test_${NAME})
ENDMACRO (CREATE_TEST)

CREATE_TEST(ip)
CREATE_TEST(logical_expression)
CREATE_TEST(sum)
