#pragma once
#include <string>
#include <vector>
class command_line
{
public:
	command_line(int argc, char* argv[]);
	~command_line();

	bool empty() const { return args.empty(); }

	size_t count() const { return args.size(); }

private:
	std::vector<std::string> args;
};

