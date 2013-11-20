#include "CommandLine.h"

command_line::command_line(int argc, char* argv[])
{
	args.reserve(argc);
	while (argc--)
	{
		args.push_back(*(argv++));
	}
}

command_line::~command_line()
{
}
