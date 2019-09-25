#include "Condition.h"

Condition::Condition ()
{
	type = 0;
	function = 0;
}

Condition::Condition (const Condition & condition)
{
	type = condition.type;
	function = condition.function;
}

Condition::~Condition ()
{
}

void Condition::set_Type (int Type)
{
	type = Type;
}

void Condition::set_function (int Function)
{
	function = Function;
}

int Condition::Type ()
{
	return type;
}

int Condition::Function ()
{
	return function;
}

void Condition::Copy (const Condition & condition)
{
	type = condition.type;
	function = condition.function;
}
