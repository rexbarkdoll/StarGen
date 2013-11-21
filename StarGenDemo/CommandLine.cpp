#include "CommandLine.h"
#include <algorithm>
#include <sstream>

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

bool command_line::get(const std::string& target) const
{
	auto it = std::find_if(args.begin(), args.end(), [&](const std::string& arg)->bool{ return (target.compare(arg) == 0); });
	if (it != args.end())
	{
		return true;
	}
	return false;
}

bool command_line::get(const std::string& target, int& value) const
{
	auto it = std::find_if(args.begin(), args.end(), [&](const std::string& arg)->bool{ return (target.compare(arg) == 0); });
	if (it != args.end())
	{
		++it;
		if (it != args.end())
		{
			return (std::stringstream(*it) >> value) ? true : false;
		}
	}
	return false;
}

bool command_line::get(const std::string& target, double& value) const
{
	auto it = std::find_if(args.begin(), args.end(), [&](const std::string& arg)->bool{ return (target.compare(arg) == 0); });
	if (it != args.end())
	{
		++it;
		if (it != args.end())
		{
			return (std::stringstream(*it) >> value) ? true : false;
		}
	}
	return false;
}
