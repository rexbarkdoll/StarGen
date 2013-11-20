#include "UnitTest++.h"
#include "CommandLine.h"

TEST(OneArgument)
{
	char* cargv[] = { "test.exe", "--a" };
	int argc = 2;
	char** argv = cargv;
	command_line cmd(--argc, ++argv);
	CHECK_EQUAL(argc, cmd.count());
}

TEST(MultipleArguments)
{
	char* cargv[] = { "test.exe", "--a", "--b", "--c" };
	int argc = 4;
	char** argv = cargv;
	command_line cmd(--argc, ++argv);
	CHECK_EQUAL(argc, cmd.count());
}
