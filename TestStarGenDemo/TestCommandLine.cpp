#include "UnitTest++.h"
#include "CommandLine.h"
#include <limits>

TEST(NoArgumentsReturnsZeroCount)
{
	char* cargv[] = { "test.exe" };
	int argc = 1;
	char** argv = cargv;
	command_line cmd(--argc, ++argv);
	CHECK(cmd.empty());
}

TEST(GetOneIntValue)
{
	char* cargv[] = { "test.exe", "--s", "123" };
	int argc = 3;
	char** argv = cargv;
	int value;
	command_line cmd(--argc, ++argv);
	CHECK_EQUAL(argc, cmd.count());
	CHECK(cmd.get("--s", value));
	CHECK_EQUAL(123, value);
}

TEST(MissingParameterIsFailure)
{
	char* cargv[] = { "test.exe", "--s" };
	int argc = 2;
	char** argv = cargv;
	int value;
	command_line cmd(--argc, ++argv);
	CHECK_EQUAL(argc, cmd.count());
	CHECK(!cmd.get("--s", value));
}

TEST(MalformedParameterIsFailure)
{
	char* cargv[] = { "test.exe", "--s", "--t" };
	int argc = 3;
	char** argv = cargv;
	int value;
	command_line cmd(--argc, ++argv);
	CHECK_EQUAL(argc, cmd.count());
	CHECK(!cmd.get("--s", value));
}

TEST(GetOneDoubleValue)
{
	char* cargv[] = { "test.exe", "--s", "123.456" };
	int argc = 3;
	char** argv = cargv;
	double value;
	command_line cmd(--argc, ++argv);
	CHECK_EQUAL(argc, cmd.count());
	CHECK(cmd.get("--s", value));
	CHECK_CLOSE(123.456, value, std::numeric_limits<double>::epsilon());
}

TEST(GetBooleanResult)
{
	char* cargv[] = { "test.exe", "--b" };
	int argc = 2;
	char** argv = cargv;
	command_line cmd(--argc, ++argv);
	CHECK_EQUAL(argc, cmd.count());
	CHECK(cmd.get("--b"));
}

TEST(GetMissingBooleanIsFailure)
{
	char* cargv[] = { "test.exe", "--a" };
	int argc = 2;
	char** argv = cargv;
	command_line cmd(--argc, ++argv);
	CHECK_EQUAL(argc, cmd.count());
	CHECK(!cmd.get("--b"));
}
