# Author: petter.strandmark@gmail.com (Petter Strandmark)

MACRO (CREATE_EXAMPLE NAME)
  ADD_EXECUTABLE(${NAME}
                 ${NAME}.cpp
                 ${EASYIP_HEADERS})
  TARGET_LINK_LIBRARIES(${NAME} easyip)
ENDMACRO (CREATE_EXAMPLE)

CREATE_EXAMPLE(knapsack)
CREATE_EXAMPLE(parking)
CREATE_EXAMPLE(sat)
CREATE_EXAMPLE(simple_example)
CREATE_EXAMPLE(sudoku)
