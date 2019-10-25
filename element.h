#pragma once

#include <string>

struct ElementData
{
	std::string n;
	int q;
	int s;

	//new stuff
	std::string symbol;
	std::string name;
	int z;
	double mass;
	// stuff that needs revisting
	std::string r_cov;
	std::string e_affinity;
	std::string p_eig;
	std::string s_eig;
	std::string abundance;
	std::string el_neg;
	std::string ion_pot;
	std::string dipol;
};

class Element
{
	ElementData _self;
public:

	Element();
	Element(const ElementData& d);
	Element(std::string s, double i, int e = 0);

	const ElementData& operator()() const;
};
