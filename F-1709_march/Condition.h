#pragma once

class Condition
{
private:
	int type;
	int function;
public:
	Condition ();
	Condition (const Condition & condition);
	~Condition ();

	void set_Type (int Type);
	void set_function (int Function);

	int Type ();
	int Function ();

	void Copy (const Condition & condition);
};