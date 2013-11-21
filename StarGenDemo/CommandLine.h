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

	bool get(const std::string& arg) const;

	bool get(const std::string& arg, int& value) const;

	bool get(const std::string& arg, double& value) const;

private:
	std::vector<std::string> args;
};

