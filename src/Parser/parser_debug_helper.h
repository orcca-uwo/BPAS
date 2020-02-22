#ifndef PARSER_DEBUG_HELPER
#define PARSER_DEBUG_HELPER

#include <stdio.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define printf_red(str){\
	printf(ANSI_COLOR_RED "%s\t\t", str);\
	printf(ANSI_COLOR_RESET);\
}
#define printf_green(str){\
	printf(ANSI_COLOR_GREEN "%s\t\t", str);\
	printf(ANSI_COLOR_RESET);\
}
#define printf_blue(str){\
	printf(ANSI_COLOR_BLUE "%s\t\t", str);\
	printf(ANSI_COLOR_RESET);\
}
#define printf_yellow(str){\
	printf(ANSI_COLOR_YELLOW "%s\t\t", str);\
	printf(ANSI_COLOR_RESET);\
}
#define printf_magenta(str){\
	printf(ANSI_COLOR_MAGENTA "%s\t\t", str);\
	printf(ANSI_COLOR_RESET);\
}
#define printf_cyan(str){\
	printf(ANSI_COLOR_CYAN "%s\t\t", str);\
	printf(ANSI_COLOR_RESET);\
}

#endif 
